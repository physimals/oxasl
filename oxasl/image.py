"""
Basic classes for ASL data processing
"""

import sys
import math
from optparse import OptionGroup

import numpy as np
import scipy

from fsl.data.image import Image
from fsl.wrappers import mcflirt, fslmaths, LOAD

from .fslwrap import Workspace

class AslOptionGroup(OptionGroup):
    """
    OptionGroup with support for ignoring certain options
    """
    def __init__(self, *args, **kwargs):
        self._ignore = kwargs.pop("ignore", [])
        OptionGroup.__init__(self, *args, **kwargs)

    def add_option(self, name=None, *args, **kwargs):
        if name not in self._ignore and name.lstrip("-") not in self._ignore and ("dest" not in kwargs or kwargs["dest"] not in self._ignore):
            OptionGroup.add_option(self, name, *args, **kwargs)

def add_data_options(parser, fname_opt="-i", output_type="directory", **kwargs):
    parser.add_option("-o", dest="output", help="Output %s" % output_type, default=None)
    parser.add_option("--debug", dest="debug", help="Debug mode", action="store_true", default=False)

    g = AslOptionGroup(parser, "Input data", **kwargs)
    g.add_option(fname_opt, dest="asldata", help="ASL data file")
    g.add_option("--order", dest="order", help="Data order as sequence of 2 or 3 characters: t=TIs/PLDs, r=repeats, p/P=TC/CT pairs. First character is fastest varying")
    g.add_option("--tis", dest="tis", help="TIs as comma-separated list")
    g.add_option("--plds", dest="plds", help="PLDs as comma-separated list")
    g.add_option("--nrpts", dest="nrpts", help="Fixed number of repeats per TI", default=None)
    g.add_option("--rpts", dest="rpts", help="Variable repeats as comma-separated list, one per TI", default=None)
    g.add_option("--iaf", dest="iaf", help="input ASl format: diff,tc,ct")
    g.add_option("--ibf", dest="ibf", help="input block format (for multi-TI): rpt,tis")
    parser.add_option_group(g)

def summary(img, log=sys.stdout):
    """
    Write a summary of the data to a file stream
    """
    log.write("%s:\n" % (img.name.ljust(30)))
    if hasattr(img, "summary"):
        img.summary(log=log)

class AslImage(Image):

    DIFFERENCED = 0
    TC_PAIRS = 1
    MULTIPHASE = 2

    """
    Subclass of fslwrap.Image which adds ASL structure information

    An AslImage contains information about the structure of the data enabling it to perform
    operations such as reordering and tag/control differencing.

    As a minimum you must provide a data order and a means of determining the number of TIs/PLDs 
    in the data. 

    Ordering is defined by a sequence of characters:
      - ``p`` - Tag/Control pairs
      - ``P`` - Control/Tag pairs
      - ``t`` - TIs/PLDs
      - ``r`` - Repeats
      - ``m`` - Multiple phases

    The sequence is in order from fastest varying (innermost grouping) to slowest varying (outermost
    grouping). If ``p/P`` is not included this describes data which is already differenced.

    Attributes:

      ``nvols`` - Number of volumes in data
      ``order`` - Data ordering string
      ``ntc`` - Number of tag/control images in data
      ``tagfirst`` - True if tag/control pairs have tag first
      ``multiphase`` - True if tag/control images are multiphase
      ``phases`` - List of phases for multiphase data
      ``ntis`` - Number of TIs/PLDs
      ``tis`` - Optional list of TIs
      ``have_plds`` - True if TIs are actually PLDs
      ``rpts`` - Repeats, one value per TI (may be constant but always stored as list)
    """
  
    def __init__(self, image, name=None, **kwargs):
        if image is None:
            raise ValueError("No image data (filename, Nibabel object or Numpy array)")

        img_args = dict([(k, v) for k, v in kwargs.items() if k in ("header")])
        Image.__init__(self, image, name=name, **img_args)
        
        order = kwargs.pop("order", None)
        iaf = kwargs.pop("iaf", None)
        ibf = kwargs.pop("ibf", None)
        ntis = kwargs.pop("ntis", None)
        nplds = kwargs.pop("nplds", None)
        tis = kwargs.pop("tis", None)
        plds = kwargs.pop("plds", None)
        nrpts = kwargs.pop("nrpts", None)
        rpts = kwargs.pop("rpts", None)
        phases = kwargs.pop("phases", None)
        nphases = kwargs.pop("nphases", None)
        
        if self.ndim == 4:
            self.nvols = self.shape[3]
        elif self.ndim == 3:
            self.nvols = 1
        else:
            raise RuntimeError("3D or 4D data expected")

        if iaf is not None:
            if order:
                raise ValueError("Can't specifiy IAF and order parameters together")
            raise RuntimeError("iaf is not implemented yet")
        elif order is None:
            #warnings.warn("Data order was not specified - assuming TC pairs in blocks of repeats")
            raise ValueError("Data order must be specified")
        self.order = order

        # Determine the number and type of tag/control images present. This may be tag/control
        # pairs (in either order), a set of multiple phases, or already differenced data
        #
        # Sets the attributes: ntc (int), tagfirst (bool), multiphase (bool), phases (list)
        if "p" in order or "P" in order:
            self.ntc = 2
            self.tagfirst = "p" in self.order
            self.multiphase = False
            self.phases = []
        elif "m" in order:
            if phases is None and nphases is None:
                raise RuntimeError("Multiphase data specified but number of phases not given")
            elif phases is not None:
                if nphases is not None and nphases != len(phases):
                    raise RuntimeError("Number of phases is not consistent with length of phases list")
            else:
                phases = [pidx * 360 / nphases for pidx in range(nphases)]

            self.phases = phases
            self.ntc = len(phases)
            self.tagfirst = False
            self.multiphase = True
        else:
            self.ntc = 1
            self.tagfirst = False
            self.multiphase = False
            self.phases = []

        # Determine the number and type of delay images present. These may be given as TIs or PLDs.
        # Internally we always refer to these as TIs with the plds attribute telling us if 
        # they are really PLDs.
        #
        # Sets the attributes tis (list), ntis (int), have_plds (bool)
        if (tis is not None and plds is not None) or (ntis is not None and nplds is not None) or \
           (tis is not None and nplds is not None) or (plds is not None and ntis is not None):
            raise RuntimeError("Cannot specify PLDs and TIs at the same time")

        self.have_plds = False

        # ntis/nplds and tis/plds are synonyms internally but we flag which we have
        if nplds is not None:
            self.have_plds = True
            ntis = nplds
        
        if plds is not None:
            tis = plds
            self.have_plds = True
        
        if ntis is None and tis is None:
            raise RuntimeError("Number of TIs/PLDs not specified")
        elif tis is not None:
            if isinstance(tis, str): tis = [float(ti) for ti in tis.split(",")]
            ntis = len(tis)
            if ntis is not None and len(tis) != ntis:
                raise RuntimeError("Number of TIs/PLDs specified as: %i, but a list of %i TIs/PLDs was given" % (ntis, len(tis)))
        self.tis = tis
        self.ntis = int(ntis)
        
        # Determine the number of repeats (fixed or variable)
        #
        # Sets the attribute rpts (list, one per TI/PLD)
        if nrpts is not None and rpts is not None:
            raise RuntimeError("Cannot specify both fixed and variable numbers of repeats")        
        elif nrpts is None and rpts is None:
            # Calculate fixed number of repeats 
            if self.nvols % (self.ntc * self.ntis) != 0:
                raise RuntimeError("Data contains %i volumes, inconsistent with %i TIs and %i tag/control images" % (self.nvols, self.ntis, self.ntc))        
            rpts = [int(self.nvols / (self.ntc * self.ntis))] * self.ntis
        elif nrpts is not None:
            nrpts = int(nrpts)
            if nrpts * self.ntis * self.ntc != self.nvols:
                raise RuntimeError("Data contains %i volumes, inconsistent with %i TIs, %i tag/control images and %i repeats" % (self.nvols, self.ntis, self.ntc, nrpts))
            rpts = [nrpts] * self.ntis
        else:
            if isinstance(rpts, str): rpts = [int(rpt) for rpt in rpts.split(",")]
            if len(rpts) != self.ntis:
                raise RuntimeError("%i TIs specified, inconsistent with %i variable repeats" % (self.ntis, len(rpts)))        
            elif sum(rpts) * self.ntc != self.nvols:
                raise RuntimeError("Data contains %i volumes, inconsistent with %i tag/control images and total of %i variable repeats" % (self.nvols, self.ntc, sum(rpts)))        
        self.rpts = rpts
        
    def _get_order_idx(self, order, tag, ti, rpt):
        idx = 0
        first = True
        for comp in order[::-1]:
            #print("comp: %s" % comp)
            if not first:
                idx *= self._get_ncomp(comp, ti)
                #print("Multiplied by %i" % self._get_ncomp(comp, ti))
            idx += self._get_comp(comp, tag, ti, rpt)
            #print("Added %i" % self._get_comp(comp, tag, ti, rpt))
            first = False
        return idx

    def _get_comp(self, comp_id, tag, ti, rpt):
        ret = {"t": ti, "r" : rpt, "p" : tag, "P" : 1-tag, "m" : tag}
        if comp_id in ret: 
            return ret[comp_id]
        else:
            raise RuntimeError("Unknown ordering character: %s" % comp_id)

    def _get_ncomp(self, comp_id, ti):
        ret = {"t": self.ntis, "r" : self.rpts[ti], "p" : 2, "P" : 2, "m" : self.ntc}
        if comp_id in ret: 
            return ret[comp_id]
        else:
            raise RuntimeError("Unknown ordering character: %s" % comp_id)

    def reorder(self, out_order):
        """
        Re-order ASL data 

        The order is defined by a string in which
        r=repeats, p=tag-control pairs, P=control-tag pairs and t=tis/plds.
        The first character is the fastest varying

        So for a data set with 3 TIs and 2 repeats an order of "ptr" would be:
        TC (TI1), TC (TI2), TC (TI3), TC(TI1, repeat 2), TC(TI2 repeat 2), etc.
        """
        if self.ntc == 1 and ("p" in out_order or "P" in out_order):
            raise RuntimeError("Data contains TC pairs but output order does not")
        elif ("p" in self.order or "P" in self.order) and ("p" not in out_order and "P" not in out_order):
            raise RuntimeError("Output order contains TC pairs but input data  does not")
        elif "m" in self.order and "m" not in out_order:
            raise RuntimeError("Data is multiphase but output order is not")
        elif "m" in out_order and "m" not in self.order:
            raise RuntimeError("Output order contains multiphases but data does not")

        #print("reordering from %s to %s" % (self.order, out_order))
        input_data = self.nibImage.get_data()
        output_data = np.zeros(self.shape, dtype=input_data.dtype)
        if input_data.ndim == 3:
            input_data = input_data[..., np.newaxis]
        tags = range(self.ntc)
        for ti in range(self.ntis):
            for rpt in range(self.rpts[ti]):
                for tag in tags:
                    #print("ti=%i, rpt=%i, tag=%i" % (ti, rpt, tag))
                    in_idx = self._get_order_idx(self.order, tag, ti, rpt)
                    #print("Input (%s) index %i" % (self.order, in_idx))
                    out_idx = self._get_order_idx(out_order, tag, ti, rpt)
                    #print("Output (%s) index %i" % (out_order, out_idx))
                    output_data[:, :, :, out_idx] = input_data[:, :, :, in_idx]
                    #print("")
        return AslImage(image=output_data, name=self.name + "_reorder", header=self.header,
                        order=out_order, tis=self.tis, ntis=self.ntis, rpts=self.rpts, phases=self.phases,
                        base=self)

    def single_ti(self, ti_idx, order=None):
        """
        Extract the subset of data for a single TI/PLD

        FIXME will not correctly set have_plds flag in output if input has PLDs
        """
        if order is None:
            if "p" in self.order or "P" in self.order: 
                order = "pr"
            elif "m" in self.order: 
                order = "mr"
            else: order = "r"
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
        output_data = reordered.nibImage.get_data()[:, :, :, start:start+nvols]
        if self.tis is not None:
            tis = [self.tis[ti_idx],]
        else:
            tis = None
        return AslImage(image=output_data, name=self.name + "_ti%i" % ti_idx, 
                        order=order, tis=tis, ntis=1, nrpts=nrpts, phases=self.phases, base=self)

    def diff(self):
        """
        Perform tag-control differencing. 
        
        Data will be reordered so the tag/control pairs are together
        """
        if "m" in self.order:
            raise RuntimeError("Cannot difference multiphase data")
        elif "p" not in self.order and "P" not in self.order:
            # Already differenced
            output_data = self.nibImage.get_data()
        elif self.nvols % 2 != 0:
            raise RuntimeError("Invalid number of volumes for TC data: %i" % self.nvols)
        else:
            output_data = np.zeros(list(self.shape[:3]) + [int(self.nvols/2)])

            # Re-order so that TC pairs are together with the tag first
            out_order = self.order.replace("p", "").replace("P", "")
            reordered = self.reorder("p" + out_order).nibImage.get_data()
            
            for t in range(int(self.nvols / 2)):
                tag = 2*t
                ctrl = tag+1
                output_data[..., t] = reordered[..., ctrl] - reordered[..., tag]
        
        out_order = self.order.replace("p", "").replace("P", "")
        return AslImage(image=output_data, name=self.name + "_diff", 
                        order=out_order, tis=self.tis, ntis=self.ntis, rpts=self.rpts, base=self)

    def mean_across_repeats(self):
        """
        Calculate the mean signal across all repeats
        """
        if self.ntc > 1:
            # Have tag-control pairs - need to diff
            diff = self.diff()
        else:
            diff = self
        
        # Reorder so repeats are together saving original order. Note that
        # rt and tr are equivalent in the output but we want to preserve
        # whatever it was beforehand
        orig_order = diff.order
        diff = diff.reorder("rt")
        input_data = diff.nibImage.get_data()

        # Create output data - one volume per ti
        output_data = np.zeros(list(self.shape[:3]) + [self.ntis])
        start = 0
        for ti, nrp in enumerate(self.rpts):
            repeat_data = input_data[..., start:start+nrp]
            output_data[..., ti] = np.mean(repeat_data, 3)
            start += nrp
        
        return AslImage(image=output_data, name=self.name + "_mean", 
                        order=orig_order, tis=self.tis, ntis=self.ntis, nrpts=1,
                        base=self)

    def perf_weighted(self):
        """
        Generate a perfusion weighted image by taking the mean over repeats and then
        the mean over TIs
        """
        meandata = self.diff().mean_across_repeats().nibImage.get_data()
        if meandata.ndim > 3:
            meandata = np.mean(meandata, axis=-1)
        return Image(image=meandata, name=self.name + "_perf_weighted", base=self)
            
    def split_epochs(self, epoch_size, overlap=0, time_order=None):
        asldata = self.diff()
        if time_order is not None:
            asldata = asldata.reorder(time_order)
        
        epoch = 0
        epoch_start = 0
        input_data = asldata.nibImage.get_data()
        ret = []
        while 1:
            epoch_end = min(epoch_start + epoch_size, asldata.nvols)
            tis = []
            rpts = []
            for ti in range(asldata.ntis):
                for rpt in range(asldata.rpts[ti]):
                    vol_idx = asldata._get_order_idx(asldata.order, 0, ti, rpt)
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
                                 order=asldata.order,
                                 tis=tis, rpts=rpts,
                                 base=asldata).mean_across_repeats()

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
        #fsl.Image.summary(self, log)
        log.write("Data shape                    : %s\n" % str(self.shape))
        log.write("%s                          : %s\n" % (ti_str, str(self.tis)))
        log.write("Number of repeats at each TI  : %s\n" % str(self.rpts))
        log.write("Label-Control                 : ")
        if self.ntc == 2:
            if self.tagfirst: log.write("Label-control pairs\n")
            else: log.write("Control-Label pairs\n")
        elif self.multiphase:
            log.write("Multiple phases (%s)" % str(self.phases))
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

        Any further keyword parameters are passed to the Image constructor
        """
        if name is None:
            name = self.name
        name = name + suffix
        if image.ndim != self.ndim or (image.ndim == 4 and image.shape[3] != self.shape[3]):
            print(kwargs)
            return Image(image=image, name=name, **kwargs)
        else:
            
            return AslImage(image=image, name=name, base=self,
                            order=self.order, ntis=self.ntis, tis=self.tis, rpts=self.rpts, phases=self.phases, **kwargs)

class AslWorkspace(Workspace):
    """
    Adds some functionality to an fslwrap.Workspace which is useful to ASL images
    """

    def smooth(self, img, fwhm, output_name=None):
        if output_name is None:
            output_name = img.name + "_smooth"
        sigma = round(fwhm/2.355, 2)

        self.log.write("Spatial smoothing with FWHM: %f (sigma=%f)\n" % (fwhm, sigma))
        smoothed = scipy.ndimage.gaussian_filter(img.nibImage.get_data(), sigma=sigma)
        return img.derived(smoothed, suffix="_smooth")

    def preprocess(self, asldata, diff=False, reorder=None, mc=False, smooth=False, fwhm=None, ref=None, **kwargs):
        self.log.write("ASL preprocessing...\n")

        # Keep original AslImage with info about TIs, repeats, etc
        orig = asldata

        if diff: 
            self.log.write("  - Tag-control subtraction\n")
            asldata = asldata.diff()
            
        if reorder:
            self.log.write("  - Re-ordering to %s\n" % reorder)
            if "p" in reorder.lower() and diff:
                reorder = reorder.replace("p", "").replace("P", "") 
            asldata = asldata.reorder(reorder)

        if mc: 
            self.log.write("  - Motion correction\n")
            output = mcflirt(asldata, cost="mutualinfo", out=LOAD)
            asldata = asldata.derived(output["out"].nibImage.get_data(), suffix="_mc")

        if smooth:
            asldata = self.smooth(asldata, fwhm=fwhm)

        self.log.write("DONE\n\n")
        return asldata

    def reg(self, wsp, ref, reg_targets, options, ref_str="asl"):
        """ 
        FIXME not functional yet
        """
        self.log.write("Segmentation and co-registration...\n")

        # Brain-extract ref image
        ref_bet = wsp.bet(ref)

        # This is done to avoid the contrast enhanced rim resulting from low intensity ref image
        d = ref_bet.nibImage.get_data()
        thr_ref = np.percentile(d[d != 0], 10.0)
        d[d < thr_ref] = 0
        raw_bet = ref_bet.derived(d, save=True)

        for imgs in reg_targets:
            reg = imgs[0]
            reg_bet = wsp.bet(reg, args="-B -f 0.3")
            name = "%s_2%s" % (reg.name, ref_str)
            name_inv = "%s_2%s" % (ref_str, reg.name)
            postreg, mat, invmat = wsp.flirt(ref_bet, args="-dof 7", 
                                             output_name=name,
                                             output_mat=name + ".mat",
                                             output_invmat=name_inv + ".mat")
            wsp.write_file(matrix_to_text(invmat), name_inv)
            for coreg in imgs[1:]:
                wsp.apply_xfm(coreg, ref_bet, name_inv, args="-interp nearestneighbour", output_name="%s_fast_seg_2asl" % t1.name)
    
        self.log.write("DONE\n\n")
