"""
Structural data module for ASL

Copyright (c) 2008-2018 University of Oxford
"""
import os

import fsl.wrappers as fsl
from fsl.data.image import Image

from .options import OptionCategory, IgnorableOptionGroup

class StructuralImageOptions(OptionCategory):
    """
    OptionGroup which contains options for describing a structural image
    """

    def __init__(self, title="Structural image", **kwargs):
        OptionCategory.__init__(self, "struc", **kwargs)
        self.title = title

    def groups(self, parser):
        group = IgnorableOptionGroup(parser, self.title, ignore=self.ignore)
        group.add_option("-s", "--struc", dest="struc", help="Structural image", type="image", default=None)
        group.add_option("--sbet", "--struc-brain", "--struc-bet", dest="struc_brain", type="image", help="Structural image (brain extracted)", default=None)
        group.add_option("--struc2asl", help="Structural->ASL transformation matrix", default=None)
        group.add_option("--asl2struc", help="ASL->Structural transformation matrix", default=None)
        group.add_option("--wm-seg", help="White matter segmentation of structural image", type="image", default=None)
        group.add_option("--gm-seg", help="Grey matter segmentation of structural image", type="image", default=None)
        group.add_option("--csf-seg", help="CSF segmentation of structural image", type="image", default=None)
        group.add_option("--fslanat", help="FSL_ANAT output directory for structural information", default=None)
        group.add_option("--fastsrc", dest="fastsrc", help="Images from a FAST segmentation - if not set FAST will be run on structural image")
        group.add_option("--senscorr", dest="senscorr", help="Use bias field (from segmentation) for sensitivity correction", action="store_true", default=False)
        group.add_option("--struc2std", help="Structural to MNI152 linear registration (.mat)")
        group.add_option("--struc2std-warp", help="Structural to MNI152 non-linear registration (warp)")

        return [group, ]

def preproc_struc(wsp):
    """
    Do preprocessing on supplied structural data - copy relevant image and do brain extraction
    """
    wsp.log.write("\nPre-processing structural data\n")
    if wsp.fslanat:
        wsp.log.write(" - Using FSL_ANAT output directory for structural data: %s\n" % wsp.fslanat)
        biascorr = os.path.join(wsp.fslanat, "T1_biascorr")
        biascorr_brain = os.path.join(wsp.fslanat, "T1_biascorr_brain")
        if os.path.isfile(biascorr) and os.path.isfile(biascorr_brain):
            wsp.log.write(" - Using bias-corrected structural images")
            wsp.struc = Image(biascorr)
            wsp.struc_brain = Image(biascorr_brain)
        else:
            wsp.log.write(" - Using non bias-corrected structural images")
            wsp.struc = Image(os.path.join(wsp.fslanat, "T1"))
            wsp.struc_brain = Image(os.path.join(wsp.fslanat, "T1_brain"))
            
        warp = os.path.join(wsp.fslanat, "T1_to_MNI_nonlin_coeff")
        mat = os.path.join(wsp.fslanat, "T1_to_MNI_lin.mat")
        if os.path.isfile(warp):
            wsp.struc2std_warp = warp
        elif os.path.isfile(mat):
            wsp.struc2std_mat = mat

    elif wsp.struc:
        wsp.log.write(" - Using structural image provided: %s\n" % wsp.struc.name)
    #elif wsp.struc_lores
    #    wsp.log.write("Low-resolution tructural image: %s\n" % wsp.struc_lores.name)
    else:
        wsp.log.write(" - No structural data supplied - output will be native space only\n")

    if wsp.struc is not None and wsp.struc_brain is None:
        wsp.log.write(" - Brain-extracting structural image\n")
        bet_result = fsl.bet(wsp.struc, output=fsl.LOAD, seg=True, mask=True, log=wsp.fsllog)
        wsp.struc_brain = bet_result["output"]
        wsp.struc_brain_mask = bet_result["output_mask"]
    elif wsp.struc_brain is not None and wsp.struc_brain_mask is None:
        wsp.struc_brain_mask = fsl.fslmaths(wsp.struc_brain).bin().run()
    # FIXME
    # wsp.preproc_struc = wsp._done

def segment(wsp):
    """
    Segment the structural image
    """
    if None in (wsp.wm_seg, wsp.gm_seg, wsp.csf_seg):
        preproc_struc(wsp)
        wsp.log.write("\nGetting structural segmentation\n")
        if wsp.fslanat:
            wsp.log.write(" - Using FSL_ANAT output\n")
            wsp.csf_pv_struc = Image(os.path.join(wsp.fslanat, "T1_fast_pve_0"))
            wsp.gm_pv_struc = Image(os.path.join(wsp.fslanat, "T1_fast_pve_1"))
            wsp.wm_pv_struc = Image(os.path.join(wsp.fslanat, "T1_fast_pve_2"))
        
            try:
                wsp.bias_struc = Image(os.path.join(wsp.fslanat, "T1_fast_bias"))
                wsp.log.write(" - Bias field extracted sucessfully")
            except:
                wsp.log.write(" - No bias field found")
        elif wsp.fastdir:
            raise NotImplementedError("Specifying FAST output directory")
        elif wsp.struc:
            wsp.log.write(" - Running FAST\n")
            fast_result = fsl.fast(wsp.struc_brain, out=fsl.LOAD, log=wsp.fsllog)
            print(fast_result)
            wsp.csf_pv_struc = fast_result["out_pve_0"]
            wsp.gm_pv_struc = fast_result["out_pve_1"]
            wsp.wm_pv_struc = fast_result["out_pve_2"]
            #wsp.bias_struc = fast_result["fast_bias"]
        else:
            raise ValueError("No structural data provided - cannot segment")

        wsp.csf_seg_struc = wsp.csf_pv_struc.data > 0.5
        wsp.gm_seg_struc = wsp.gm_pv_struc.data > 0.5
        wsp.wm_seg_struc = wsp.wm_pv_struc.data > 0.5
