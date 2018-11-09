"""
Subclass of fsl.data.image.Image which represents ASL data
"""

import sys
import warnings

import numpy as np

from fsl.data.image import Image

from .options import OptionCategory, IgnorableOptionGroup

class AslImageOptions(OptionCategory):
    """
    OptionGroup which contains options for describing an ASL image
    """

    def __init__(self, title="Input ASL image", fname_opt="-i", **kwargs):
        OptionCategory.__init__(self, "image", **kwargs)
        self.title = title
        self.fname_opt = fname_opt

    def groups(self, parser):
        group = IgnorableOptionGroup(parser, self.title, ignore=self.ignore)
        group.add_option("--asldata", self.fname_opt, help="ASL data file")
        group.add_option("--iaf", help="input ASl format: diff=differenced,tc=tag-control,ct=control-tag,mp=multiphase,ve=vessel-encoded")
        group.add_option("--order", help="Data order as sequence of 2 or 3 characters: t=TIs/PLDs, r=repeats, l=labelling (tag/control/phases etc). First character is fastest varying")
        group.add_option("--tis", help="TIs as comma-separated list")
        group.add_option("--plds", help="PLDs as comma-separated list - alternative to --tis")
        group.add_option("--ntis", help="Number of TIs (for use when processing does not require actual values)", type="int")
        group.add_option("--nplds", help="Equivalent to --ntis", type="int")
        group.add_option("--rpts", help="Variable repeats as comma-separated list, one per TI/PLD", default=None)
        group.add_option("--nphases", help="For --iaf=mp, number of phases (assumed to be evenly spaced)", default=None, type="int")
        group.add_option("--nenc", help="For --iaf=ve, number of encoding cycles", type="int", default=8)
        group.add_option("--casl", help="Acquisition was pseudo cASL (pcASL) rather than pASL", action="store_true", default=False)
        group.add_option("--tau", "--taus", "--bolus", help="Bolus duration. Can be single value or comma separated list, one per TI/PLD")
        group.add_option("--slicedt", help="Timing difference between slices (ms) for 2D readout", type=float, default=0.0)
        group.add_option("--sliceband", help="Number of slices per pand in multi-band setup", type=int)
        group.add_option("--artsupp", help="Arterial suppression (vascular crushing) was used", action="store_true", default=False)
        group.add_option("--ibf", help="input block format - alternative to --order for compatibility. rpt=Blocks of repeats (i.e. repeats are slowest varying), tis=Blocsk of TIs/PLDs")
        return [group, ]

def summary(img, log=sys.stdout):
    """
    Write a summary of the data to a file stream
    """
    log.write("%s:\n" % (img.name.ljust(30)))
    if hasattr(img, "summary"):
        img.summary(log=log)

class AslImage(Image):
    """
    Subclass of fsl.data.image.Image which adds ASL structure information

    An AslImage contains information about the structure of the data enabling it to perform
    operations such as reordering and label/control differencing.

    As a minimum you must provide a means of determining the number of TIs/PLDs in the data. 

    Specifying the data format and ordering explicitly is recommended, but a default ordering will be
    used (with a warning) if you do not.

    Ordering can be defined in two ways: 
    
    1. Setting the ``order`` parameters to a sequence of characters (case insensitive):
      - ``l`` - Labelling images (e.g. label/control pairs, sequence of multi-phases, vessel encoding cycles)
      - ``t`` - TIs/PLDs
      - ``r`` - Repeats

      The sequence is in order from fastest varying (innermost grouping) to slowest varying (outermost
      grouping). If ``p/P`` is not included this describes data which is already differenced.

    2. Specifying the ``ibf`` option
      - ``ibf`` - ``rpt`` Blocked by repeats, i.e. first repeat of all TIs, followed by second repeat of all TIs...
                  ``tis`` Blocked by TIs/PLDs, i.e. all repeats of first TI, followed by all repeats of second TI...
                  When using --ibf, the labelling images (e.g. label/control pairs) are always adjacent
    
    The data format is defined using the ``iaf`` parameter:

      - ``iaf`` - ``tc`` = tag then control, ``ct`` = control then tag, ``mp`` = multiphase, ``ve`` = vessel encoded

    Attributes:

      ``nvols`` - Number of volumes in data
      ``iaf`` - Data format - see above
      ``order`` - Data ordering string - see above
      ``ntc`` - Number of labelling images in data (e.g. 2 for TC pairs, 1 for differenced data)
      ``phases`` - List of phases for multiphase data (``iaf='mp'``)
      ``ntis`` - Number of TIs/PLDs
      ``tis`` - Optional list of TIs
      ``plds`` - Optional list of PLDs
      ``have_plds`` - True if PLDs were provided
      ``tau`` - Bolus durations - one per TI/PLD. If ``have_plds`` is True, tis are derived by adding the bolus duration to the PLDs
      ``rpts`` - Repeats, one value per TI (may be given as a constant but always stored as list)
    """
  
    DIFFERENCED = 0
    TC_PAIRS = 1
    MULTIPHASE = 2

    def __init__(self, image, name=None, **kwargs):
        if image is None:
            raise ValueError("No image data (filename, Nibabel object or Numpy array)")

        # This is sort-of a bug in nibabel or fslpy - it passes the kwargs to the
        # nibabel.load function which does not expect extra keyword arguments
        img_kwargs = ("header", "xform", "loadData", "calcRange", "indexed", "threaded", "dataSource")
        img_args = dict([(k, v) for k, v in kwargs.items() if k in img_kwargs])
        Image.__init__(self, image, name=name, **img_args)
        
        order = kwargs.pop("order", None)
        iaf = kwargs.pop("iaf", None)
        ibf = kwargs.pop("ibf", None)
        ntis = kwargs.pop("ntis", None)
        nplds = kwargs.pop("nplds", None)
        tis = kwargs.pop("tis", None)
        plds = kwargs.pop("plds", None)
        rpts = kwargs.pop("rpts", kwargs.pop("nrpts", None))
        phases = kwargs.pop("phases", None)
        nphases = kwargs.pop("nphases", None)
        nenc = kwargs.pop("nenc", None)
        
        if self.ndim == 4:
            self.nvols = self.shape[3]
        elif self.ndim == 3:
            self.nvols = 1
        else:
            raise ValueError("3D or 4D data expected")

        # Determine the data format and ordering
        #
        # Sets the attributes: iaf (str), order (str)
        if not iaf:
            if order is not None and "l" in order:
                warnings.warn("Data format was not specified - assuming TC pairs")
                iaf = "tc"
            elif not order:
                warnings.warn("Data format was not specified - assuming differenced")
                iaf = "diff"
            else:
                # Order specified and did not include labelling images so we are entitled
                # to assume differenced data without a warning
                iaf = "diff"
        elif iaf not in ("diff", "tc", "ct", "mp", "ve"):
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

        if iaf != "diff" and "l" not in order:
            order = "l" + order
        for char in order.lower():
            if char not in ('l', 'r', 't'):
                raise ValueError("Unrecognized character in data ordering: '%s'" % char)

        self.order = order.lower()
        self.iaf = iaf.lower()

        # Determine the number of labelling images present. This may be label/control
        # pairs, a set of multiple phases, vessel encoding cycles or already differenced data
        #
        # Sets the attributes: ntc (int), phases (list or None)
        self.phases = None
        if self.iaf in ("tc", "ct"):
            self.ntc = 2
        elif self.iaf == "mp":
            if phases is None and nphases is None:
                raise ValueError("Multiphase data specified but number of phases not given")
            elif phases is not None:
                if nphases is not None and nphases != len(phases):
                    raise ValueError("Number of phases is not consistent with length of phases list")
            else:
                phases = [pidx * 360 / nphases for pidx in range(nphases)]

            if isinstance(phases, str): phases = [float(ph) for ph in phases.split(",")]
            self.phases = phases
            self.ntc = len(phases)
        elif self.iaf == "ve":
            if nenc is None:
                raise ValueError("Vessel encoded data specified but number of encoding cycles not given")
            self.ntc = nenc
        else:
            self.ntc = 1

        # Determine the number and type of delay images present. These may be given as TIs or PLDs.
        #
        # Internally we always have TIs, and only have PLDs as well if they were provided. If PLDs
        # were provided, the TIs are derived by adding on the bolus duration
        #
        # Sets the attributes tis (list), ntis (int), have_plds (bool), plds (list)
        if (tis is not None and plds is not None) or (ntis is not None and nplds is not None) or \
           (tis is not None and nplds is not None) or (plds is not None and ntis is not None):
            raise ValueError("Cannot specify PLDs and TIs at the same time")

        self.have_plds = False

        if nplds is not None:
            ntis = nplds
            self.have_plds = True
        
        if plds is not None:
            tis = plds
            self.have_plds = True

        if ntis is None and tis is None:
            raise ValueError("Number of TIs/PLDs not specified")
        elif tis is not None:
            if isinstance(tis, str): tis = [float(ti) for ti in tis.split(",")]
            ntis = len(tis)
            if ntis is not None and len(tis) != ntis:
                raise ValueError("Number of TIs/PLDs specified as: %i, but a list of %i TIs/PLDs was given" % (ntis, len(tis)))
        self.tis = tis
        self.ntis = int(ntis)
        if self.have_plds:
            self.plds = tis

        if ibf_guessed and len(self.tis) > 1:
            warnings.warn("Data order was not specified for multi-TI data - assuming blocks of repeats")
                
        # Determine the number of repeats (fixed or variable)
        #
        # Sets the attribute rpts (list, one per TI/PLD)
        if rpts is None:
            # Calculate fixed number of repeats 
            if self.nvols % (self.ntc * self.ntis) != 0:
                raise ValueError("Data contains %i volumes, inconsistent with %i TIs and %i labelling images" % (self.nvols, self.ntis, self.ntc))        
            rpts = [int(self.nvols / (self.ntc * self.ntis))] * self.ntis
        else:
            if isinstance(rpts, str): rpts = [int(rpt) for rpt in rpts.split(",")]
            elif isinstance(rpts, int): rpts = [rpts,]
            if len(rpts) == 1:
                rpts *= self.ntis
            elif len(rpts) != self.ntis:
                raise ValueError("%i TIs specified, inconsistent with %i variable repeats" % (self.ntis, len(rpts)))        
            elif sum(rpts) * self.ntc != self.nvols:
                raise ValueError("Data contains %i volumes, inconsistent with %i labelling images and total of %i repeats" % (self.nvols, self.ntc, sum(rpts)))        
        self.rpts = rpts

        # Bolus durations should be a sequence same length as TIs/PLDs
        #
        # Sets the attributes taus (list)
        self.taus = kwargs.pop("bolus", kwargs.pop("taus", kwargs.pop("tau", None)))
        if self.taus is None:
            self.taus = 1.8
        if isinstance(self.taus, str): self.taus = [float(tau) for tau in self.taus.split(",")]
        elif isinstance(self.taus, (float, int)): self.taus = [float(self.taus),] * self.ntis

        if len(self.taus) == 1:
            self.taus = self.taus * self.ntis
        if len(self.taus) != self.ntis:
            raise ValueError("%i bolus durations specified, inconsistent with %i TIs/PLDs" % (len(self.taus), self.ntis))
           
        # Labelling type. CASL/pCASL normally associated with PLDs but can pass TIs instead. 
        # However we would not expect PLDs for a non-CASL aquisition so this generates a warning
        #
        # If PLDs were provided, TIs are derived by adding the bolus duration to the PLDs
        #
        # Sets the attributes casl (bool), updates tis (list)
        self.casl = kwargs.pop("casl", None)
        if self.casl is None: 
            self.casl = self.have_plds
        if self.have_plds:
            if not self.casl:
                warnings.warn("PLDs specified but aquisition was not CASL/pCASL - will treat these as TIs")
            else:
                self.tis = [pld + tau for pld, tau in zip(self.plds, self.taus)]

        # Other acquisition parameters
        self.slicedt = kwargs.pop("slicedt", 0)
        self.sliceband = kwargs.pop("sliceband", None)
        self.artsupp = kwargs.pop("artsupp", False)

    def get_vol_index(self, label_idx, ti_idx, rpt_idx, order=None):
        """
        Get the volume index for a specified label, TI and repeat index

        :param label_idx: Label index starting from 0, e.g. for ``iaf=ct`` 0 would be the control image,
                          for ``iaf=mp`` 3 would be the 4th phase encoded image
        :param ti_idx: TI/PLD index, starting from 0
        :param rpt_idx: Repeat index, starting from 0
        :param order: If specified use custom data ordering string (does not change ordering
                      within this AslImage - use ``reorder`` for that)
        """
        if order is None:
            order = self.order
        if len(order) < 3: order = "l" + order

        # Not yet found a simple way to do this which works for variable repeats
        # without crude iteration!
        ti_pos, rpt_pos, label_pos = order.index("t"), order.index("r"), order.index("l")
        its = [0, 0, 0]
        vol_idx = 0
        while 1:
            ti, rpt, label = its[ti_pos], its[rpt_pos], its[label_pos]
            if (ti, rpt, label) == (ti_idx, rpt_idx, label_idx):
                return vol_idx
            its[0] += 1
            if its[0] == self._get_ncomp(order[0], ti):
                its[0] = 0
                its[1] += 1
            if its[1] == self._get_ncomp(order[1], ti):
                its[1] = 0
                its[2] += 1
            if its[2] == self._get_ncomp(order[2], ti):
                # Normally this is the end but we might have run out of 
                # repeats but not run out of TIs
                pass
            else:
                vol_idx += 1
            if vol_idx > self.nvols:
                break

        raise ValueError("No volume for supplied TI, label and repeat")

    def _get_ncomp(self, comp_id, ti):
        ret = {"t": self.ntis, "r" : self.rpts[ti], "l" : self.ntc}
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
        
        input_data = self.data
        if input_data.ndim == 3:
            input_data = input_data[..., np.newaxis]
        output_data = np.zeros(input_data.shape, dtype=input_data.dtype)
        tags = range(self.ntc)
        for ti in range(self.ntis):
            for rpt in range(self.rpts[ti]):    
                for tag in tags:
                    if iaf != self.iaf:
                        # Change from TC to CT or vice versa
                        out_tag = 1-tag
                    else:
                        out_tag = tag
                    in_idx = self.get_vol_index(tag, ti, rpt)
                    out_idx = self.get_vol_index(out_tag, ti, rpt, order=out_order)
                    output_data[:, :, :, out_idx] = input_data[:, :, :, in_idx]

        if not name:
            name = self.name + "_reorder"
        return self.derived(image=output_data, name=name, iaf=iaf, order=out_order)
        
    def single_ti(self, ti_idx, order=None, name=None):
        """
        Extract the subset of data for a single TI/PLD

        FIXME will not correctly set have_plds flag in output if input has PLDs
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
            start += self.rpts[idx]*self.ntc
        nrpts = self.rpts[ti_idx]
        nvols = nrpts * self.ntc
        output_data = reordered.data[:, :, :, start:start+nvols]
        tis, plds = None, None
        if self.have_plds and self.plds is not None:
            plds = [self.plds[ti_idx],]
        elif not self.have_plds and self.tis is not None:
            tis = [self.tis[ti_idx],]

        if self.taus is not None:
            taus = self.taus[ti_idx]
        else:
            taus = None
        
        if not name:
            name = self.name + "_ti%i" % ti_idx
        return self.derived(image=output_data, name=name, order=order, tis=tis, plds=plds, taus=taus, ntis=1, rpts=nrpts)
        
    def diff(self, name=None):
        """
        Perform tag-control subtraction. 
        
        Data will be reordered so the tag/control pairs are together
        """
        if self.iaf == "diff":
            # Already differenced
            return self
        elif self.iaf not in ("tc", "ct"):
            raise ValueError("Data is not tag-control pairs - cannot difference")
        else:
            output_data = np.zeros(list(self.shape[:3]) + [int(self.nvols/2)])

            # Re-order so that TC pairs are together with the tag first
            out_order = self.order.replace("l", "")
            reordered = self.reorder("l" + out_order, iaf="tc").data
            
            for vol in range(int(self.nvols / 2)):
                tag = 2*vol
                ctrl = tag+1
                output_data[..., vol] = reordered[..., ctrl] - reordered[..., tag]
        
        out_order = self.order.replace("l", "")
        
        if not name:
            name = self.name + "_diff"
        return self.derived(image=output_data, name=name, iaf="diff", order=out_order)

    def mean_across_repeats(self, name=None, diff=True):
        """
        Calculate the mean ASL signal across all repeats

        :return: Label-control subtracted AslImage with one volume per TI/PLD
        """
        if diff and self.ntc > 1:
            # Have tag-control pairs - need to subtract
            data = self.diff()
            out_order = "rt"
        elif self.ntc > 1:
            data = self
            out_order = "rlt"
        else:
            data = self
            out_order = "rt"

        # Reorder so repeats are together saving original order. Note that
        # rt and tr are equivalent in the output but we want to preserve
        # whatever it was beforehand
        orig_order = data.order
        data = data.reorder(out_order)
        input_data = data.data
        if input_data.ndim == 3:
            input_data = input_data[..., np.newaxis]

        # Create output data - one repeat per ti
        output_data = np.zeros(list(self.shape[:3]) + [self.ntis * data.ntc])
        start = 0
        for ti, nrp in enumerate(self.rpts):
            for label in range(data.ntc):
                repeat_data = input_data[..., start:start+nrp]
                output_data[..., label+data.ntc*ti] = np.mean(repeat_data, 3)
                start += nrp
        
        if not name:
            name = self.name + "_mean"
        return self.derived(image=output_data, name=name, iaf=data.iaf, order=orig_order, rpts=1)

    def mean(self, name=None):
        """
        Take the mean across all volumes

        This takes a naive mean without differencing or grouping by TIs

        :return: 3D Image. Not an AslImage as timing information lost
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

        :return: 3D Image. Not an AslImage as timing information lost
        """
        meandata = self.diff().mean_across_repeats().data
        if meandata.ndim > 3:
            meandata = np.mean(meandata, axis=-1)
        if not name:
            name = self.name + "_pwi"
        return Image(image=meandata, name=name, header=self.header)
            
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
            for ti in range(asldata.ntis):
                for rpt in range(asldata.rpts[ti]):
                    vol_idx = asldata.get_vol_index(0, ti, rpt)
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
                                 tis=tis, rpts=rpts,
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
        """
        ti_str = "TIs "
        if self.have_plds: ti_str = "PLDs"
        log.write("Data shape                    : %s\n" % str(self.shape))
        log.write("%s                          : %s\n" % (ti_str, str(self.tis)))
        log.write("Number of repeats at each TI  : %s\n" % str(self.rpts))
        log.write("Label-Control                 : ")
        if self.ntc == 2:
            if self.iaf == "tc": log.write("Label-control pairs\n")
            else: log.write("Control-Label pairs\n")
        elif self.iaf == "mp":
            log.write("Multiple phases (%s)" % str(self.phases))
        elif self.iaf == "ve":
            log.write("Vessel encoded (%i encoding cycles)" % self.ntc)
        else:
            log.write("Already differenced\n")

    def derived(self, image, name=None, suffix="", **kwargs):
        """
        Create a derived ASL image based on this one, but with different data

        This is only possible if the number of volumes match, otherwise we cannot
        use the existing information about TIs, repeats etc. If the number of volumes
        do not match a generic Image is returned instead
        
        :param data: Numpy data for derived image
        :param name: Name for new image (can be simple name or full filename)
        :param suffix: If name not specified, construct by adding suffix to original image name

        Any further keyword parameters are passed to the AslImage constructor, overriding existing
        attributes, so this can be used to create a derived image with different numbers of 
        repeats, etc, provided the data is consistent with this.

        If the AslImage constructor fails, a basic fsl.data.image.Image is returned.
        """
        if name is None:
            name = self.name
        name = name + suffix

        DERIVED_ATTRS = ["iaf", "order", "rpts", "taus", "phases",
                         "casl", "sliceband", "slicedt", "artsupp"]
        if self.have_plds:
            DERIVED_ATTRS.append("plds")
            DERIVED_ATTRS.append("nplds")
        else:
            DERIVED_ATTRS.append("tis")
            DERIVED_ATTRS.append("ntis")

        derived_kwargs = {}
        for attr in DERIVED_ATTRS:
            derived_kwargs[attr] = kwargs.get(attr, getattr(self, attr, None))
        if self.iaf == "ve":
            derived_kwargs["nenc"] = self.ntc

        try:
            return AslImage(image=image, name=name, header=self.header, **derived_kwargs)
        except ValueError as exc:
            warnings.warn("AslImage.derived failed (%s) - returning fsl.data.image.Image" % str(exc))
            return Image(image=image, name=name, header=self.header, **kwargs)
            