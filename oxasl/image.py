"""
OXASL -  Classes for representing ASL data and constructing instances from command line parameters

Copyright (c) 2008 University of Nottingham
"""
import collections
import json
import os
import sys
import warnings

import six
import numpy as np
import scipy.linalg

from fsl.data.image import Image

from .options import OptionCategory, OptionGroup

try:
    import oxasl_bids
except ImportError:
    oxasl_bids = None

class Options(OptionCategory):
    """
    OptionGroup which contains options for describing an ASL image
    """

    def __init__(self, fname_opt="-i"):
        OptionCategory.__init__(self, "image")
        self.fname_opt = fname_opt

    def groups(self, parser):
        group = OptionGroup(parser, "ASL data")
        group.add_option("--asldata", self.fname_opt, help="ASL data file")
        group.add_option("--iaf", help="input ASl format: diff=differenced,tc=tag-control,ct=control-tag,mp=multiphase,ve=vessel-encoded,vediff=pairwise subtracted VE data, hadamard=Hadamard time encoded data")
        group.add_option("--order", help="Data order as sequence of 2 or 3 characters: t=TIs/PLDs, r=repeats, l=labelling (tag/control/phases etc). First character is fastest varying")
        group.add_option("--tis", help="TIs (s) as comma-separated list")
        group.add_option("--tes", help="TSs (s) as comma-separated list for multi-TE ASL data")
        group.add_option("--plds", help="PLDs (s) as comma-separated list - alternative to --tis")
        group.add_option("--ntis", help="Number of TIs (for use when processing does not require actual values)", type="int")
        group.add_option("--nplds", help="Equivalent to --ntis", type="int")
        group.add_option("--rpts", help="Variable repeats as comma-separated list, one per TI/PLD", default=None)
        group.add_option("--nphases", help="For --iaf=mp, number of phases (assumed to be evenly spaced)", default=None, type="int")
        group.add_option("--nenc", help="For --iaf=ve/vediff, number of encoding cycles in original acquisition", type="int", default=8)
        group.add_option("--hadamard-size", help="Size of Hadamard encoding matrix, must be power of 2 or 12", type="int", default=8)
        group.add_option("--casl", help="Acquisition was pseudo cASL (pcASL) rather than pASL", action="store_true", default=False)
        group.add_option("--tau", "--taus", "--bolus", "--lds", help="Labelling duration (s). Can be single value or comma separated list, one per TI/PLD or per sub-bolus for Hadamard data")
        group.add_option("--slicedt", help="Timing difference between slices (s) for 2D readout", type=float, default=0.0)
        group.add_option("--sliceband", help="Number of slices per pand in multi-band setup", type=int)
        group.add_option("--artsupp", help="Arterial suppression (vascular crushing) was used", action="store_true", default=False)
        group.add_option("--ibf", help="input block format - alternative to --order for compatibility. rpt=Blocks of repeats (i.e. repeats are slowest varying), tis=Blocsk of TIs/PLDs")
        group.add_option("--calib-first-vol", help="First volume of the data is an M0 image for calibration - remainder is ASL data", action="store_true", default=False)
        return [group, ]

def summary(img, log=sys.stdout):
    """
    Write a summary of an Image to a stream

    For an AslImage the ``summary`` method is used which displays all the ASL data in human
    readable form. For generic Image objects just basic information is displayed

    :param img: fsl.data.image.Image object which may or may not be an AslImage
    :param log: Stream-like object - default is ``sys.stdout``
    """
    log.write("%s:\n" % (img.name.ljust(30)))
    if hasattr(img, "summary"):
        img.summary(log=log)

def data_order(iaf, ibf, order, multite=False):
    """
    Determine the data format and ordering from ``iaf`` and ``ibf`` options

    If ``iaf`` is not specified (None or empty string), TC pairs are specified if the
    ``order`` parameter is given and contains the ``l`` character. Otherwise differenced
    data is assumed.

    If neither ``ibf`` nor ``order`` are specified, ``ibf=rpt`` is guessed.

    If ``iaf`` is not ``diff`` and ``order`` does not include the ``l` character it is assumed
    that the encoded images (e.g. TC pairs) are the fastest varying.

    If any parameters had to be guessed (rather than inferred from other information) a
    warning is output.

    :return: Tuple of: IAF (ASL format, 'tc', 'ct', 'diff', 'quant', 've', 'vediff', 'hadamard' or 'mp'), ordering (sequence of
             2 or three chars, 'l'=labelling images, 'r'=repeats, 't'=TIs/PLDs, 'e'=TEs. The characters
             are in order from fastest to slowest varying), Boolean indicating whether the block
             format needed to be guessed
    """
    if order is not None:
        order = order.lower()

    if not iaf:
        if order is not None and "l" in order:
            # Order specified labelling images so guess TC pairs as this is most common
            warnings.warn("Data format was not specified - assuming TC pairs")
            iaf = "tc"
        elif not order:
            # No order specified  - guess differenced for compatibility with Oxford ASL
            warnings.warn("Data format was not specified - assuming differenced")
            iaf = "diff"
        else:
            # Order specified and did not include labelling images so we are entitled
            # to assume differenced data without a warning
            iaf = "diff"
    elif iaf not in ("diff", "tc", "ct", 'quant', "mp", "ve", "vediff", "hadamard"):
        raise ValueError("Unrecognized data format: iaf=%s" % iaf)

    ibf_guessed = False
    if not order:
        if not ibf:
            # Guess but defer warning until we have extracted TIs as it doesn't matter
            # for single TI data
            ibf = "rpt"
            ibf_guessed = True

        order_map = {
            "rpt" : "tr",
            "tis" : "rt",
        }
        order = order_map.get(ibf.lower(), None)
        if not order:
            raise ValueError("Unrecognized data block format: ibf=%s" % ibf)

    if iaf not in ("diff", 'quant') and "l" not in order:
        order = "l" + order

    if multite and "e" not in order:
        # TE ordering is optional since most data is single TE. Insert it at the
        # start as this seems to make most sense in terms of multi-TE acquisition
        order = "e" + order

    for char in order:
        if char not in ('l', 'r', 't', 'e'):
            raise ValueError("Unrecognized character in data ordering: '%s'" % char)

    return iaf, order, ibf_guessed

class AslImage(Image):
    """
    Subclass of fsl.data.image.Image which adds ASL structure information

    An AslImage contains information about the structure of the data enabling it to perform
    operations such as reordering and label/control differencing.

    As a *minimum* you must provide a means of determining the number of TIs/PLDs in the data.
    Specifying the data format and ordering explicitly is recommended, but a default ordering will be
    used (with a warning) if you do not.

    Ordering can be defined in two ways:

    1. Setting the ``order`` parameters to a sequence of characters (case insensitive):
      - ``l`` - Labelling images (e.g. label/control pairs, sequence of multi-phases, vessel encoding cycles)
      - ``t`` - TIs/PLDs
      - ``r`` - Repeats
      - ``e`` - TEs (optional)

      The sequence is in order from fastest varying (innermost grouping) to slowest varying (outermost
      grouping). If ``l`` is not included this describes data which is already differenced.

    2. Specifying the ``ibf`` option
      - ``ibf`` - ``rpt`` Blocked by repeats, i.e. first repeat of all TIs, followed by second repeat of all TIs...
                  ``tis`` Blocked by TIs/PLDs, i.e. all repeats of first TI, followed by all repeats of second TI...
                  When using --ibf, the labelling images (e.g. label/control pairs) are always adjacent

    The data format is defined using the ``iaf`` parameter:

      - ``iaf`` - ``tc`` = tag then control, ``ct`` = control then tag, ``mp`` = multiphase,
        ``ve`` = vessel encoded, ``vediff`` = Pairwise subtracted vessel encoded, ``diff`` = already differenced,
        ``hadamard` = Hadamard encoded, ``quant`` = already quantified

    If oxasl_bids is installed, some/all of the image metadata can be extracted from a .json sidecar if one exists

    :ivar nvols:  Number of volumes in data
    :ivar iaf:    Data format - see above
    :ivar order:  Data ordering string - see above
    :ivar ntc:    Number of labelling images in data (e.g. 2 for TC pairs, 1 for differenced data,
                  for vessel encoded or multiphase data the number of encoding cycles/phases, for Hadamard data the
                  size of the Hadamard matrix, e.g. 8)
    :ivar casl:   If True, acquisition was CASL/pCASL rather than PASL
    :ivar phases: List of phases for multiphase data (``iaf='mp'``)
    :ivar ntis:   Number of TIs/PLDs
    :ivar tis:    Optional list of TIs
    :ivar plds:    Optional list of PLDs
    :ivar have_plds: True if PLDs were provided rather than TIs
    :ivar tau:    Bolus durations - one per TI/PLD. If ``have_plds`` is True, tis are derived by adding the bolus duration to the PLDs
    :ivar rpts:   Repeats, one value per TI (may be given as a constant but always stored as list)
    :ivar slicedt: Increase in TI/PLD per slice in z dimension for 2D acquisitions (default: 0)
    :ivar sliceband: Number of slices per band for multiband acquisitions (default: None)
    :ivar artsupp: If True, data was acquired with arterial suppression
    """

    DIFFERENCED = 0
    TC_PAIRS = 1
    MULTIPHASE = 2

    def __init__(self, image, name=None, **kwargs):
        if image is None:
            raise ValueError("No image data (filename, Nibabel object or Numpy array)")

        # This is sort-of a bug in nibabel or fslpy - it passes the kwargs to the
        # nibabel.load function which does not expect extra keyword arguments. So we will extract
        # those that fslpy/nibabel might expect and only pass these
        img_kwargs = ("header", "xform", "loadData", "calcRange", "indexed", "threaded", "dataSource")
        img_args = dict([(k, v) for k, v in kwargs.items() if k in img_kwargs])

        if kwargs.pop("calib_first_vol", False):
            # First volume is an M0 calibration image - stip it off and initialize the image with the
            # remaining data
            temp_img = Image(image, name="temp", **img_args)
            img_args.pop("header", None)
            self.setMeta("calib", Image(temp_img.data[..., 0], name="calib", header=temp_img.header, **img_args))
            Image.__init__(self, temp_img.data[..., 1:], name=name, header=temp_img.header, **img_args)
        else:
            self.setMeta("calib", kwargs.pop("calib", None))
            Image.__init__(self, image, name=name, **img_args)

        # If we are initializing from a filename and we have oxasl_bids we can get metadata
        # from the BIDS sidecar
        if oxasl_bids is not None and kwargs.pop("use_bids", False) and isinstance(image, str):
            if ".nii" in image:
                sidecar_fname = image[:image.index(".nii")] + ".json"
            else:
                sidecar_fname = image + ".json"
            if os.path.exists(sidecar_fname):
                with open(sidecar_fname, "r") as f:
                    metadata = json.load(f)
                    metadata["img_shape"] = self.shape
                kwargs.update(oxasl_bids.oxasl_config_from_metadata(metadata, "asl"))

        order = kwargs.pop("order", None)
        iaf = kwargs.pop("iaf", None)
        ibf = kwargs.pop("ibf", None)
        ntis = kwargs.pop("ntis", None)
        nplds = kwargs.pop("nplds", None)
        tis = kwargs.pop("tis", None)
        tes = kwargs.pop("tes", None)
        plds = kwargs.pop("plds", None)
        rpts = kwargs.pop("rpts", kwargs.pop("nrpts", None))
        phases = kwargs.pop("phases", None)
        nphases = kwargs.pop("nphases", None)
        nenc = kwargs.pop("nenc", kwargs.get("ntc", None))
        had_size = kwargs.pop("hadamard_size", kwargs.pop("nenc", kwargs.pop("ntc", None)))

        if self.ndim not in (3, 4):
            raise ValueError("3D or 4D data expected")

        # Check for multi-TE data
        #
        # Note that for the moment we assume that multi-TEs are applied at every TI/PLD.
        # Variable repeats/taus are only per-TI. Default is one TE with value 0
        #
        # Sets the attributes tes, ntes
        if tes is not None:
            if isinstance(tes, six.string_types):
                tes = [float(te) for te in tes.split(",")]
            self.setMeta("tes", tes)
            self.setMeta("ntes", len(tes))
        else:
            self.setMeta("tes", [0,])
            self.setMeta("ntes", 1)

        # Determine the data ordering
        #
        # Sets the metadata: iaf (str), order (str)
        iaf, order, ibf_guessed = data_order(iaf, ibf, order, self.ntes > 1)
        self.setMeta("order", order)
        self.setMeta("iaf", iaf.lower())

        # Determine the number of labelling images present. This may be label/control
        # pairs, a set of multiple phases, vessel encoding cycles or already differenced data
        #
        # Sets the metadata: ntc (int), phases (list or None)
        if self.iaf in ("tc", "ct"):
            ntc = 2
        elif self.iaf == "mp":
            if phases is None and nphases is None:
                raise ValueError("Multiphase data specified but number of phases not given")
            elif phases is not None:
                if nphases is not None and nphases != len(phases):
                    raise ValueError("Number of phases is not consistent with length of phases list")
            else:
                phases = [pidx * 360 / nphases for pidx in range(nphases)]

            if isinstance(phases, six.string_types):
                phases = [float(ph) for ph in phases.split(",")]
            ntc = len(phases)
        elif self.iaf in ("ve", "vediff"):
            if nenc is None:
                raise ValueError("Vessel encoded data specified but number of encoding cycles not given")
            elif self.iaf == "ve":
                ntc = nenc
            elif nenc % 2 == 0:
                ntc = int(nenc / 2)
            else:
                raise ValueError("Subtracted vessel encoded data must have had an even number of encoding cycles")
        elif self.iaf == "hadamard":
            if had_size is None:
                raise ValueError("Hadamard encoded data specified but hadamard size not given")
            elif had_size not in (2, 4, 8, 12, 16, 32, 64):
                raise ValueError("Hadamard encoded data specified but hadamard size (%s) not a supported value" % had_size)
            ntc = had_size
        else:
            ntc = 1
        self.setMeta("ntc", ntc)
        self.setMeta("phases", phases)

        # Determine the number and type of delay images present. These may be given as TIs or PLDs.
        #
        # Internally we always have TIs, and only have PLDs as well if they were provided. If PLDs
        # were provided, the TIs are derived by adding on the bolus duration
        #
        # Sets the attributes tis (list), ntis (int), have_plds (bool), plds (list)
        if tis is not None and plds is not None:
            raise ValueError("Cannot specify PLDs and TIs at the same time")

        have_plds = False
        if nplds is not None:
            ntis = nplds
            have_plds = True

        if plds is not None:
            tis = plds
            have_plds = True

        if ntis is None and tis is None and self.iaf != "quant":
            raise ValueError("Number of TIs/PLDs not specified")
        elif tis is not None:
            print(tis, type(tis))
            if isinstance(tis, six.string_types):
                tis = [float(ti) for ti in tis.split(",")]
            if ntis is not None and len(tis) != ntis:
                raise ValueError("Number of TIs/PLDs specified as: %i, but a list of %i TIs/PLDs was given" % (ntis, len(tis)))
            if len(tis) == self.nvols and self.iaf in ("tc", "ct"):
                # Special case: we have label/control pairs with a TI/PLD for each volume.
                # However label/control TIs/PLDs must be equal in ASL and we expect a single value
                # for each pair. Check this is the case and reduce the number (note that
                # this scenario is common with BIDS datasets or JSON sidecars)
                tis_even = tis[::2]
                tis_odd = tis[1::2]
                if not np.allclose(tis_even, tis_odd):
                    raise ValueError(f"Inconsistent TIs/PLDs for label/control pairs: {tis_even} vs {tis_odd}")
                tis = tis_even
            ntis = len(tis)

        if self.iaf == "quant":
            if (ntis is not None and ntis > 1) or (tis is not None and tis != [0]):
                raise ValueError("Cannot specify TIs/PLDs for quantified data")
            ntis, tis = 1, [0]

        self.setMeta("ntis", int(ntis))
        self.setMeta("have_plds", have_plds)
        if have_plds:
            self.setMeta("plds", tis)
        else:
            self.setMeta("tis", tis)

        if ibf_guessed and len(self.tis) > 1:
            warnings.warn("Data order was not specified for multi-TI data - assuming blocks of repeats")

        # Determine the number of repeats (fixed or variable)
        #
        # Sets the attribute rpts (list, one per TI/PLD)
        if self.iaf == "quant" and self.nvols > 1:
            raise ValueError("Cannot have multi-volume data if iaf=quant")

        nvols_norpts = self.ntc * self.ntis * self.ntes
        if rpts is None:
            # Calculate fixed number of repeats
            if self.nvols % nvols_norpts != 0:
                raise ValueError("Data contains %i volumes, inconsistent with %i TIs and %i labelling images at %i TEs" % (self.nvols, self.ntis, self.ntc, self.ntes))
            rpts = [int(self.nvols / nvols_norpts)] * self.ntis
        else:
            if isinstance(rpts, six.string_types):
                rpts = [int(rpt) for rpt in rpts.split(",")]
            elif isinstance(rpts, (int, np.integer)):
                rpts = [rpts,]

            if len(rpts) == 1:
                rpts *= self.ntis
            elif len(rpts) != self.ntis:
                raise ValueError("%i TIs specified, inconsistent with %i variable repeats" % (self.ntis, len(rpts)))

            if sum(rpts) * self.ntc * self.ntes != self.nvols:
                raise ValueError("Data contains %i volumes, inconsistent with %i labelling images at %i TEs and total of %i repeats" % (self.nvols, self.ntc, self.ntes, sum(rpts)))
        self.setMeta("rpts", rpts)

        # Bolus durations should be a sequence same length as TIs/PLDs
        #
        # Sets the attributes taus (list, one per TI/PLD)
        taus = kwargs.pop("bolus", kwargs.pop("taus", kwargs.pop("tau", None)))
        if taus is None:
            taus = 1.8
        if isinstance(taus, six.string_types):
            taus = [float(tau) for tau in taus.split(",")]
        elif isinstance(taus, (float, int, np.integer)):
            taus = [float(taus),] * self.ntis

        if len(taus) == 1:
            taus = taus * ntis
        if self.iaf == "hadamard":
            # With Hadamard encoding we might have one tau per PLD as usual, but we might also have hybrid
            # protocols with 1 tau per sub-bolus. We want to expand this out to 1 per sub-bolus per PLD
            if len(taus) == self.ntis and len(taus) == self.ntc-1:
                warnings.warn("Number of LDs matches both number of PLDs and sub-boli. This is ambiguous, we will assume they are matched to PLDs")
            if len(taus) == self.ntis:
                # one per PLD - duplicate for each sub-bolus
                taus = sum([[t] * (self.ntc-1) for t in taus], [])
            elif len(taus) == self.ntc - 1:
                # One per sub-bolus - duplicate for each PLD
                taus = list(taus) * self.ntis
            elif len(taus) == self.ntis * (self.ntc-1):
                # One per PLD/sub-bolus - no need to change
                pass
            else:
                raise ValueError("%i bolus durations specified, inconsistent with %i TIs/PLDs and %i Hadamard sub-boli" % (len(taus), self.ntis, self.ntc-1))
        elif len(taus) != self.ntis:
            raise ValueError("%i bolus durations specified, inconsistent with %i TIs/PLDs" % (len(taus), self.ntis))
        self.setMeta("taus", taus)

        # Labelling type. CASL/pCASL normally associated with PLDs but can pass TIs instead.
        # However we would not expect PLDs for a non-CASL aquisition so this generates a warning
        #
        # If PLDs were provided, TIs are derived by adding the bolus duration to the PLDs
        #
        # Sets the attributes casl (bool), updates tis (list)
        casl = kwargs.pop("casl", None)
        if casl is None:
            casl = self.have_plds
        if self.have_plds:
            if not casl:
                warnings.warn("PLDs specified but aquisition was not CASL/pCASL - will treat these as TIs")
        self.setMeta("casl", casl)

        # Other acquisition parameters
        self.setMeta("slicedt", kwargs.pop("slicedt", 0))
        self.setMeta("sliceband", kwargs.pop("sliceband", None))
        self.setMeta("artsupp", kwargs.pop("artsupp", False))

    def __getattr__(self, name):
        return self.getMeta(name, None)

    @property
    def nvols(self):
        if self.ndim == 4:
            return self.shape[3]
        else:
            return 1

    @property
    def sliceband(self):
        sb = self.getMeta("sliceband", None)
        if sb is not None:
            return int(sb)
        else:
            return None

    @property
    def tis(self):
        if self.iaf == "hadmard":
            raise RuntimeError("TIs not defined for Hadamard encoded data")
        elif self.casl and self.have_plds:
            return [pld + tau for pld, tau in zip(self.plds, self.taus)]
        else:
            return self.getMeta("tis", self.getMeta("plds"))

    @property
    def nplds(self):
        return self.ntis

    @property
    def nphases(self):
        return self.getMeta("ntc", None)

    @property
    def nenc(self):
        return self.getMeta("ntc", None)

    @property
    def hadamard_size(self):
        return self.getMeta("ntc", None)

    def is_var_repeats(self):
        """
        :return: True if this data set has repeats which vary by time point
        """
        return min(self.rpts) != max(self.rpts)

    def get_vol_index(self, label_idx, ti_idx, rpt_idx, te_idx=0, order=None):
        """
        Get the volume index for a specified label, TI and repeat index

        :param label_idx: Label index starting from 0, e.g. for ``iaf=ct`` 0 would be the control image,
                          for ``iaf=mp`` 3 would be the 4th phase encoded image
        :param ti_idx: TI/PLD index, starting from 0
        :param rpt_idx: Repeat index, starting from 0
        :param te_idx: TE index for multi-TE data, starting from 0
        :param order: If specified use custom data ordering string (does not change ordering
                      within this AslImage - use ``reorder`` for that)
        """
        if order is None:
            order = self.order
        else:
            order = order.lower()
        
        if "l" not in order:
            # Need labelling image ordering even if differenced - harmless if ntc == 1
            order = "l" + order
        
        if "e" not in order:
            # Need ordering for TEs even if data is single-TE - this is harmless if ntes == 1
            order = "e" + order

        if ti_idx >= self.ntis:
            raise ValueError("Requested TI index %i but only %i TIs present" % (ti_idx, self.ntis))
        if te_idx >= self.ntes:
            raise ValueError("Requested TE index %i but only %i TEs present" % (te_idx, self.ntes))

        # Not yet found a simple way to do this which works for variable repeats
        # without crude iteration!
        ti_pos, rpt_pos, label_pos, te_pos = order.index("t"), order.index("r"), order.index("l"), order.index("e")
        its = [0, 0, 0, 0]
        vol_idx = 0
        while 1:
            ti, rpt, label, te = its[ti_pos], its[rpt_pos], its[label_pos], its[te_pos]
            if (ti, rpt, label, te) == (ti_idx, rpt_idx, label_idx, te_idx):
                return vol_idx
            its[0] += 1
            if its[0] == self._get_ncomp(order[0], ti):
                its[0] = 0
                its[1] += 1
            if its[1] == self._get_ncomp(order[1], ti):
                its[1] = 0
                its[2] += 1
            if its[2] == self._get_ncomp(order[2], ti):
                its[2] = 0
                its[3] += 1
            if its[3] == self._get_ncomp(order[3], ti):
                # Normally this is the end but we might have run out of
                # repeats but not run out of TIs
                pass
            else:
                vol_idx += 1
            if vol_idx > self.nvols:
                break

        raise ValueError("No volume for supplied TI, TE, label and repeat")

    def _get_ncomp(self, comp_id, ti):
        ret = {"t": self.ntis, "r" : self.rpts[ti], "l" : self.ntc, "e" : self.ntes}
        if comp_id in ret:
            return ret[comp_id]
        else:
            raise RuntimeError("Unknown ordering character: %s" % comp_id)

    def reorder(self, out_order=None, iaf=None, name=None):
        """
        Re-order ASL data

        The order is defined by a string in which
        r=repeats, l=labelling images, and t=tis/plds.
        The first character is the fastest varying

        So for a tag-control data set with 3 TIs and 2 repeats an order of "ltr" would be:
        TC (TI1), TC (TI2), TC (TI3), TC(TI1, repeat 2), TC(TI2 repeat 2), etc.
        """
        if out_order is None:
            out_order = self.order

        if iaf is None:
            iaf = self.iaf

        if self.iaf == "diff" and "l" in out_order:
            raise ValueError("Data is differenced but output order is not")
        elif "l" not in out_order and self.iaf != "diff":
            raise ValueError("Data is not differenced but output_order does not contain labelling")
        elif iaf != self.iaf and (iaf not in ("tc", "ct") or self.iaf not in ("tc", "ct")):
            raise ValueError("Can't change data format from '%s' to '%s'" % (self.iaf, iaf))
        
        if self.ntes > 1 and "e" not in out_order:
            out_order = "e" + out_order

        input_data = self.data
        if input_data.ndim == 3:
            input_data = input_data[..., np.newaxis]
        output_data = np.zeros(input_data.shape, dtype=input_data.dtype)
        tags = range(self.ntc)
        for te in range(self.ntes):
            for ti in range(self.ntis):
                for rpt in range(self.rpts[ti]):
                    for tag in tags:
                        if iaf != self.iaf:
                            # Change from TC to CT or vice versa
                            out_tag = 1-tag
                        else:
                            out_tag = tag
                        in_idx = self.get_vol_index(tag, ti, rpt, te)
                        out_idx = self.get_vol_index(out_tag, ti, rpt, te, order=out_order)
                        output_data[:, :, :, out_idx] = input_data[:, :, :, in_idx]

        if not name:
            name = self.name + "_reorder"
        return self.derived(image=output_data, name=name, iaf=iaf, order=out_order)

    def single_ti(self, ti_idx, order=None, name=None):
        """
        Extract the subset of data for a single TI/PLD

        FIXME will not correctly set have_plds flag in output if input has PLDs

        :param ti_idx: Index of the TI/PLD required starting from 0
        :param order: If specified the ordering string for the returned data
        :param name: Optional name for returned image. Defaults to original name with
                     suffix ``_ti<n>`` where ``<n>`` is the TI/PLD index
        :return: AslImage instance containing data only for the selected TI/PLD
        """
        if order is None:
            if self.iaf == "diff":
                order = "r"
            else:
                order = "lr"
        elif "t" in order:
            order = order.remove("t")
        order = order + "t"

        # Re-order so that TIs are together
        reordered = self.reorder(order)

        # Find the start index for this TI and the number of times it was repeated
        start = 0
        for idx in range(ti_idx):
            start += self.rpts[idx]*self.ntc*self.ntes
        nrpts = self.rpts[ti_idx]
        nvols = nrpts*self.ntc*self.ntes
        output_data = reordered.data[:, :, :, start:start+nvols]
        tis, plds = None, None
        if self.plds is not None:
            plds = [self.plds[ti_idx],]
        if self.tis is not None:
            tis = [self.tis[ti_idx],]

        if self.taus is not None:
            if self.iaf == "hadamard":
                # With Hadamard we have 1 tau per sub-bolus
                taus = self.taus[ti_idx:ti_idx+self.ntc-1]
            else:
                taus = self.taus[ti_idx]
        else:
            taus = None

        if not name:
            name = self.name + "_ti%i" % ti_idx
        return self.derived(image=output_data, name=name, order=order, tis=tis, plds=plds, taus=taus, ntis=1, nplds=1, rpts=nrpts)

    def diff(self, name=None):
        """
        Perform tag-control subtraction.

        Data will be reordered so the tag/control pairs are together

        Note that currently differencing is not supported for multiphase or vessel encoded data.
        For Hadamard data, the output will be the decoded sub-boli images.

        :param name: Optional name for returned image. Defaults to original name with suffix ``_diff``
        :return: AslImage instance containing differenced data
        """
        extra_kwargs = {}
        if self.iaf in ("diff", "quant"):
            # Already differenced
            return self
        elif self.iaf in ("tc", "ct"):
            output_data = np.zeros(list(self.shape[:3]) + [int(self.nvols/self.ntc)])

            # Re-order so that TC pairs are together with the tag first
            out_order = self.order.replace("l", "")
            differenced_order = out_order
            reordered = self.reorder("l" + out_order, iaf="tc").data

            for vol in range(int(self.nvols / 2)):
                tag = 2*vol
                ctrl = tag+1
                output_data[..., vol] = reordered[..., ctrl] - reordered[..., tag]
        elif self.iaf == "hadamard":
            # Hadamard decoding
            had_matrix = scipy.linalg.hadamard(self.ntc)
            output_data = np.zeros(list(self.shape[:3]) + [(self.ntc-1)*int(self.nvols/self.ntc)])

            # Re-order so that Hadamard encoding images are together and PLDs are next varying.
            # We will re-order again at the end to make sure final order is consistent with
            # the input order.
            out_order = self.order.replace("l", "")
            differenced_order = "tr"
            reordered = self.reorder("ltr").data

            # Decode each block of N Hadamard images into N-1 perfusion weighted images
            # Assumes that the Hadamard images are in adjacent volumes but we don't need to
            # worry at this point about what PLD/repeat they correspond to (see below)
            # Each column of the Hadamard matrix (excluding the first) contains +1/-1 values
            # which form a linear combination of the images in the sequence to extract the
            # decoded sub-bolus for that column number. Note that the convention is that +1
            # represents a label condition and thus should be subtracted
            for vol_idx in range(0, self.nvols):
                rpt_idx = int(vol_idx / self.ntc)
                had_idx = vol_idx % self.ntc
                for sub_bolus in range(self.ntc-1):
                    output_data[..., rpt_idx*(self.ntc-1)+sub_bolus] -= had_matrix[had_idx, sub_bolus+1]*reordered[..., vol_idx]
            output_data /= (self.ntc/2)

            # Now calculate the PLDs for the data decoded about which vary by sub-bolus
            # because for each sub-bolus, all of the sub-boli *after* it will contribute
            # to the PLD. Note that we have by now expanded the taus list so it is one
            # value per sub-bolus if it wasn't originally (not per PLD)
            #
            # For this to work with the decoded data we must have had all the original PLDs
            # together, hence why we reordered to 'ltr' ordering
            decoded_plds, decoded_taus, decoded_rpts = [], [], []
            for pld_idx in range(0, self.nplds):
                encoding_taus = self.taus[pld_idx * (self.ntc-1): (pld_idx + 1) * (self.ntc-1)]
                for subbolus_idx in range(self.ntc-1):
                    total_pld = self.plds[pld_idx]
                    for subsequent_subbolus_idx in range(subbolus_idx+1, self.ntc-1):
                        total_pld += encoding_taus[subsequent_subbolus_idx]

                    decoded_plds.append(total_pld)
                    decoded_taus.append(encoding_taus[subbolus_idx])
                    decoded_rpts.append(self.rpts[pld_idx])

            extra_kwargs["plds"] = decoded_plds
            extra_kwargs["rpts"] = decoded_rpts
            extra_kwargs["taus"] = decoded_taus
            extra_kwargs["nplds"] = self.nplds * (self.ntc-1)
        else:
            raise ValueError("Data is not tag-control pairs - cannot difference")

        if not name:
            name = self.name + "_diff"
        ret = self.derived(image=output_data, name=name, iaf="diff", order=differenced_order, **extra_kwargs)

        # It is possible that to do decoding we might have changed the data order - so put this right
        # if necessary (will do nothing if the ordering is already correct)
        return ret.reorder(out_order)

    def mean_across_repeats(self, name=None, diff=True):
        """
        Calculate the mean ASL signal across all repeats

        :param name: Optional name for returned image. Defaults to original name with suffix ``_mean``
        :param diff: If False do not difference the data before taking the mean (default: True)
        :return: Label-control subtracted AslImage with one volume per TI/PLD
        """
        if diff and self.ntc > 1:
            # Have tag-control pairs - need to subtract
            data = self.diff()
            out_order = "ret"
        elif self.ntc > 1:
            data = self
            out_order = "relt"
        else:
            data = self
            out_order = "ret"

        # Reorder so repeats are together saving original order. Note that
        # rt and tr are equivalent in the output but we want to preserve
        # whatever it was beforehand
        orig_order = data.order
        data = data.reorder(out_order)
        input_data = data.data
        if input_data.ndim == 3:
            input_data = input_data[..., np.newaxis]

        # Create output data - one repeat per ti. Note ordering of nested loops follows
        # out_order defined above
        output_data = np.zeros(list(data.shape[:3]) + [data.ntis * data.ntc * data.ntes])
        start = 0
        for ti, nrp in enumerate(data.rpts):
            for label in range(data.ntc):
                for te in range(data.ntes):
                    repeat_data = input_data[..., start:start+nrp]
                    output_data[..., te+data.ntes*(label+data.ntc*ti)] = np.mean(repeat_data, 3)
                    start += nrp

        if not name:
            name = self.name + "_mean"
        return data.derived(image=output_data, name=name, iaf=data.iaf, order=orig_order, rpts=1)

    def mean(self, name=None):
        """
        Take the mean across all volumes

        This takes a naive mean without differencing or grouping by TIs

        :param name: Optional name for returned image. Defaults to original name with suffix ``_mean``
        :return: 3D fsl.data.image.Image. Not an AslImage as timing information has been lost
        """
        meandata = self.data
        if meandata.ndim > 3:
            meandata = np.mean(meandata, axis=-1)
        if not name:
            name = self.name + "_mean"
        return Image(image=meandata, name=name, header=self.header)

    def perf_weighted(self, name=None):
        """
        Generate a perfusion weighted image by taking the mean over repeats and then
        the mean over TIs

        :param name: Optional name for returned image. Defaults to original name with suffix ``_pwi``
        :return: 3D fsl.data.image.Image. Not an AslImage as timing information lost
        """
        
        if not name:
            name = self.name + "_pwi"

        if (self.iaf == "mp") or (self.iaf == "ve"):
            # Special case for multiphase and vessel-encoded data - we cannot
            # difference all time points but we can take the difference between
            # maximal in/out phase using FFT for multiphase data or SD across
            # encoding cycles for vessel-encoded data. We do not make this part
            # of diff() because it is not a 'proper' differencing and is just a
            # workaround for generating a PWI, e.g. for registration

            # First take the mean across repeats, leaving TIs (PLDs) as the last dimension
            meandata = self.mean_across_repeats(diff=False).reorder(out_order="lrt")
            
            # If there are multiple TIs (PLDs), then average across these first
            # NB. Uses a similar approach to the init-loc functionality in oxasl_ve
            if self.ntis > 1:                
                mean_acrosstis = np.zeros(list(meandata.data.shape[:3]) + [self.ntc], dtype=np.float32)
                for idx in range(self.ntis):
                    mean_acrosstis += meandata[..., idx*self.ntc:(idx+1)*self.ntc]
                mean_acrosstis /= self.ntis
            else:
                mean_acrosstis = meandata.data

            # Now run specific processing for multiphase or vessel-encoded data
            if self.iaf == "mp":
                # Take the FFT across the multiphase dimension
                fft = np.fft.fft(mean_acrosstis, axis=-1)

                # The first (non-DC) component represents the sinusoidal
                # variation across cycles which gives us an approximation to the
                # perfusion signal
                mean_diffdata = np.abs(fft[..., 1])

            elif self.iaf == "ve":
                # Special case for vessel-encoded data. In this case, without
                # knowing the specifics of the vessel-encoding procedure, we
                # make the (reasonable) assumption that each vessel is
                # reasonably efficiently encoded across the vessel-encoding
                # cycles, in which case the standard deviation across cycles
                # will give a reasonable representation of the perfusion signal,
                # irrespective of the order in which the vessels have been
                # labelled/controlled.
                mean_diffdata = np.std(mean_acrosstis,axis=-1)                

        else: # Standard differencing
            mean_diffdata = self.diff().mean_across_repeats().data
            if mean_diffdata.ndim > 3:
                mean_diffdata = np.mean(mean_diffdata, axis=-1)


        ret = Image(image=mean_diffdata, name=name, header=self.header)
        return ret

    def split_epochs(self, epoch_size, overlap=0, time_order=None):
        """
        Split ASL data into 'epochs' of a specified size, with optional overlap
        """
        asldata = self.diff()
        if time_order is not None:
            asldata = asldata.reorder(time_order)

        epoch = 0
        epoch_start = 0
        input_data = asldata.data
        ret = []
        while 1:
            epoch_end = min(epoch_start + epoch_size, asldata.nvols)
            tis = []
            rpts = []
            for te in range(asldata.ntes):
                for ti in range(asldata.ntis):
                    for rpt in range(asldata.rpts[ti]):
                        vol_idx = asldata.get_vol_index(0, ti, rpt, te)
                        if vol_idx >= epoch_start and vol_idx < epoch_end:
                            ti_val = self.tis[ti]
                            if ti_val not in tis:
                                tis.append(ti_val)
                                rpts.append(1)
                            else:
                                idx = tis.index(ti_val)
                                rpts[idx] += 1
            epoch_data = input_data[..., epoch_start:epoch_end]
            epoch_img = AslImage(image=epoch_data,
                                 name=self.name + "_epoch%i" % epoch,
                                 iaf=asldata.iaf, order=asldata.order,
                                 tis=tis, rpts=rpts, tes=self.tes,
                                 header=self.header).mean_across_repeats()

            ret.append(epoch_img)
            epoch += 1
            epoch_start += epoch_size - overlap
            # Finish if start of next epoch is out of range or if current epoch
            # ended at the last volume
            if epoch_start >= asldata.nvols or epoch_end == asldata.nvols:
                break

        return ret

    def summary(self, log=sys.stdout):
        """
        Write a summary of the data to a file stream

        :param log: Stream-like object to write the summary to. Defaults to ``sys.stdout``
        """
        for key, value in self.metadata_summary().items():
            log.write("%s: %s\n" % (key.ljust(30), str(value)))

    def metadata_summary(self):
        """
        Generate a human-readable dictionary of metadata

        :return: Dictionary mapping human readable metadata name to value (e.g. 'Label type'
                 might map to 'Label-control pairs'). The keys and values are not guaranteed
                 to be fixed and should not be parsed - use the instance attributes instead.
        """
        md = collections.OrderedDict()
        md["Data shape"] = str(self.shape)

        if self.iaf == "tc":
            label_type = "Label-control pairs"
        elif self.iaf == "ct":
            label_type = "Control-Label pairs"
        elif self.iaf == "mp":
            label_type = "Multiphase"
        elif self.iaf == "quant":
            label_type = "Already quantified"
        elif self.iaf == "ve":
            label_type = "Vessel encoded"
        elif self.iaf == "vediff":
            label_type = "Vessel encoded (pairwise subtracted)"
        elif self.iaf == "hadamard":
            label_type = "Hadamard time-encoded"
        else:
            label_type = "Already differenced"
        md["Label type"] = label_type

        if self.iaf == "quant":
            # Other parameters basically irrelevant if data is already quantified
            return md

        if self.iaf == "mp":
            md["Phases"] = str(self.phases)
        elif self.iaf in ("ve", "vediff"):
            md["Encoding cycles"] = self.nenc

        md["Labelling"] = "CASL/pCASL" if self.casl else "PASL"
        if self.have_plds:
            md["PLDs (s)"] = str(self.plds)
        else:
            md["TIs (s)"] = str(self.tis)

        md["Repeats at each TI"] = str(self.rpts)
        if self.taus:
            md["Bolus durations (s)"] = str(self.taus)
        if self.slicedt:
            md["Time per slice (s)"] = self.slicedt
        if self.sliceband:
            md["Slices per band"] = self.sliceband
        if self.artsupp:
            md["Arterial suppression"] = "Was used"
        
        if len(self.tes) > 0:
            md["TEs (s)"] = str(self.tes)
            
        return md

    def derived(self, image, name=None, suffix="", **kwargs):
        """
        Create a derived ASL image based on this one, but with different data

        This is only possible if the number of volumes match, otherwise we cannot
        use the existing information about TIs, repeats etc. If the number of volumes
        do not match a generic Image is returned instead

        Any further keyword parameters are passed to the AslImage constructor, overriding existing
        attributes, so this can be used to create a derived image with different numbers of
        repeats, etc, provided the data is consistent with this.

        Note that the image space does not have to match, however in this case the 'header'
        kwarg should be passed to create the new image header

        :param data: Numpy data for derived image
        :param name: Name for new image (can be simple name or full filename)
        :param suffix: If name not specified, construct by adding suffix to original image name
        :return: derived AslImage. However if the AslImage constructor fails, a basic
                 fsl.data.image.Image is returned and a warning is output.
        """
        if name is None:
            name = self.name
        name = name + suffix

        derived_attrs = ["iaf", "order", "rpts", "taus", "phases", "tes",
                         "casl", "sliceband", "slicedt", "artsupp", "calib"]
        if self.have_plds:
            derived_attrs.append("plds")
            derived_attrs.append("nplds")
        else:
            derived_attrs.append("tis")
            derived_attrs.append("ntis")

        derived_kwargs = {}
        for attr in derived_attrs:
            derived_kwargs[attr] = kwargs.get(attr, getattr(self, attr, None))
        if derived_kwargs["iaf"] in ("ve", "vediff") and "nenc" not in derived_kwargs:
            derived_kwargs["nenc"] = self.nenc
        if derived_kwargs["iaf"] == "hadamard" and "hadamard_size" not in derived_kwargs:
            derived_kwargs["hadamard_size"] = self.hadamard_size

        try:
            return AslImage(image=image, name=name, header=kwargs.get("header", self.header), **derived_kwargs)
        except ValueError as exc:
            if kwargs.get("raise_on_fail", False):
                raise
            else:
                warnings.warn("AslImage.derived failed (%s) - returning fsl.data.image.Image" % str(exc))
                return Image(image=image, name=name, header=kwargs.get("header", self.header), **kwargs)
