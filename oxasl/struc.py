"""
Structural data module for ASL

Copyright (c) 2008-2018 University of Oxford
"""
import os
import glob

import numpy as np

import fsl.wrappers as fsl
from fsl.data.image import Image
from fsl.utils.path import PathError

from oxasl.options import OptionCategory, IgnorableOptionGroup
from oxasl.reporting import LightboxImage

class StructuralImageOptions(OptionCategory):
    """
    OptionGroup which contains options for describing a structural image
    """

    def __init__(self, title="Structural image", **kwargs):
        OptionCategory.__init__(self, "struc", **kwargs)
        self.title = title

    def groups(self, parser):
        group = IgnorableOptionGroup(parser, self.title, ignore=self.ignore)
        group.add_option("--struc", "-s", help="Structural image", type="image", default=None)
        group.add_option("--struc-brain", "--sbet", "--struc-bet", type="image", help="Structural image (brain extracted)", default=None)
        group.add_option("--struc2asl", help="Structural->ASL transformation matrix", default=None)
        group.add_option("--asl2struc", help="ASL->Structural transformation matrix", default=None)
        group.add_option("--wm-seg", help="White matter segmentation of structural image", type="image", default=None)
        group.add_option("--gm-seg", help="Grey matter segmentation of structural image", type="image", default=None)
        group.add_option("--csf-seg", help="CSF segmentation of structural image", type="image", default=None)
        group.add_option("--fslanat", help="FSL_ANAT output directory for structural information", default=None)
        group.add_option("--fastsrc", help="Images from a FAST segmentation - if not set FAST will be run on structural image")
        group.add_option("--struc2std", help="Structural to MNI152 linear registration (.mat)")
        group.add_option("--struc2std-warp", help="Structural to MNI152 non-linear registration (warp)")

        return [group, ]

def init(wsp):
    """
    Do initialization on supplied structural data - copy relevant image and do brain extraction

    FIXME copy across all supplied structural data
    """
    if wsp.structural is not None:
        return

    wsp.log.write("\nInitialising structural data\n")
    wsp.sub("structural")

    if wsp.fslanat:
        wsp.log.write(" - Using FSL_ANAT output directory for structural data: %s\n" % wsp.fslanat)
        biascorr = os.path.join(wsp.fslanat, "T1_biascorr")
        biascorr_brain = os.path.join(wsp.fslanat, "T1_biascorr_brain")
        if glob.glob(biascorr + ".*") and glob.glob(biascorr_brain + ".*"):
            wsp.log.write(" - Using bias-corrected structural images\n")
            wsp.structural.struc = Image(biascorr)
            wsp.structural.brain = Image(biascorr_brain)
        else:
            wsp.log.write(" - Using non bias-corrected structural images\n")
            wsp.structural.struc = Image(os.path.join(wsp.fslanat, "T1"))
            wsp.structural.brain = Image(os.path.join(wsp.fslanat, "T1_brain"))
    elif wsp.struc:
        wsp.log.write(" - Using structural image provided by user: %s\n" % wsp.struc.name)
        wsp.structural.struc = wsp.struc
        if wsp.wm_seg is not None:
            wsp.structural.wm_seg = wsp.wm_seg # FIXME user override segmentation

    #elif wsp.structural.struc_lores
    #    wsp.log.write("Low-resolution tructural image: %s\n" % wsp.structural.struc_lores.name)
    else:
        wsp.log.write(" - No structural data supplied - output will be native space only\n")

    if wsp.structural.struc is not None and wsp.structural.brain is None:
        wsp.log.write(" - Brain-extracting structural image\n")
        bet_result = fsl.bet(wsp.structural.struc, output=fsl.LOAD, seg=True, mask=True, log=wsp.fsllog)
        wsp.structural.brain = bet_result["output"]
        #wsp.structural.brain_mask = bet_result["output_mask"]

    if wsp.structural.brain is not None and wsp.structural.brain_mask is None:
        # FIXME - for now get the mask by binarising the brain image but gives slightly
        # different results compared to using the mask returned by BET
        wsp.structural.brain_mask = Image((wsp.structural.brain.data != 0).astype(np.int), header=wsp.structural.struc.header)

    if wsp.structural.struc is not None:
        segment(wsp)

def segment(wsp):
    """
    Segment the structural image
    """
    init(wsp)
    if None in (wsp.structural.wm_seg, wsp.structural.gm_seg, wsp.structural.csf_seg):
        init(wsp)
        page = wsp.report.page("seg")
        page.heading("Segmentation of structural image")

        wsp.log.write("\nGetting structural segmentation\n")
        if wsp.fslanat:
            wsp.log.write(" - Using FSL_ANAT output\n")
            page.text("Segmentation taken from FSL_ANAT output at ``%s``" % wsp.fslanat)
            wsp.structural.csf_pv = Image(os.path.join(wsp.fslanat, "T1_fast_pve_0"))
            wsp.structural.gm_pv = Image(os.path.join(wsp.fslanat, "T1_fast_pve_1"))
            wsp.structural.wm_pv = Image(os.path.join(wsp.fslanat, "T1_fast_pve_2"))

            try:
                wsp.structural.bias = Image(os.path.join(wsp.fslanat, "T1_fast_bias"))
                wsp.log.write(" - Bias field extracted sucessfully\n")
            except PathError:
                wsp.log.write(" - No bias field found")
        elif wsp.fastsrc:
            # FIXME should be possible to implement this
            raise NotImplementedError("Specifying FAST output directory")
            #img = os.path.split(wsp.fastsrc)[1]
            #wsp.structural.csf_pv = os.path.join(wsp.fastsrc, "%s_pve_0" % img)
            #wsp.structural.gm_pv = os.path.join(wsp.fastsrc, "%s_pve_1" % img)
            #wsp.structural.wm_pv = os.path.join(wsp.fastsrc, "%s_pve_2" % img)
        elif wsp.structural.struc:
            wsp.log.write(" - Running FAST\n")
            page.text("FAST run to segment structural image")
            fast_result = fsl.fast(wsp.structural.brain, out=fsl.LOAD, log=wsp.fsllog)
            wsp.structural.csf_pv = fast_result["out_pve_0"]
            wsp.structural.gm_pv = fast_result["out_pve_1"]
            wsp.structural.wm_pv = fast_result["out_pve_2"]
            #wsp.bias_struc = fast_result["fast_bias"]
        else:
            raise ValueError("No structural data provided - cannot segment")

        wsp.structural.csf_seg = Image((wsp.structural.csf_pv.data > 0.5).astype(np.int), header=wsp.structural.struc.header)
        wsp.structural.gm_seg = Image((wsp.structural.gm_pv.data > 0.5).astype(np.int), header=wsp.structural.struc.header)
        wsp.structural.wm_seg = Image((wsp.structural.wm_pv.data > 0.5).astype(np.int), header=wsp.structural.struc.header)

        page.heading("Segmentation image", level=1)
        page.text("CSF partial volume")
        page.image("csf_pv", LightboxImage(wsp.structural.csf_pv, bgimage=wsp.structural.brain))
        page.text("Grey matter partial volume")
        page.image("gm_pv", LightboxImage(wsp.structural.gm_pv, bgimage=wsp.structural.brain))
        page.text("White matter partial volume")
        page.image("wm_pv", LightboxImage(wsp.structural.wm_pv, bgimage=wsp.structural.brain))
