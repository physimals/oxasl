"""
OXASL - Quantification using Fabber

Includes partial volume correction

Copyright (c) 2008 University of Nottingham
"""

import numpy as np

from fsl.data.image import Image

from oxasl import mask, reg
from . import multistep_fit


def run(wsp):
    # Basic non-PVC run
    multistep_fit.run(wsp.sub("basil"))
    wsp.quantify_wsps.append("basil")

    # Re-do registration using PWI as reference
    reg.run(wsp, redo=True, struc_bbr=True, struc_flirt=False, use_quantification_wsp=wsp.basil)

    # Quantification in alternate spaces
    for quantify_space in ("struc", "std", "custom"):
        if wsp.ifnone("quantify_%s" % quantify_space, False):
            quantify_name = "basil_%s" % quantify_space
            quantify_wsp = wsp.sub(quantify_name)
            quantify_wsp.image_space = quantify_space
            multistep_fit.run(quantify_wsp) 
            wsp.quantify_wsps.append(quantify_name)

    # If the user has provided manual PV maps (pvgm and pvgm) then do PVEc, even if they
    # have not explicitly given the --pvcorr option 
    wsp.user_pv_flag = ((wsp.pvwm is not None) and (wsp.pvgm is not None))
    if wsp.pvcorr or wsp.user_pv_flag:
        # Partial volume correction is very sensitive to the mask, so recreate it
        # if it came from the structural image as this requires accurate ASL->Struc registration
        if wsp.rois.mask_src == "struc":
            wsp.rois.mask_orig = wsp.rois.mask
            wsp.rois.mask = None
            mask.run(wsp)

        if wsp.pvcorr or wsp.user_pv_flag:
            _default_pvcorr(wsp)

def _default_pvcorr(wsp):
    # Do partial volume correction fitting
    #
    # FIXME: We could at this point re-apply all corrections derived from structural space?
    # But would need to make sure corrections module re-transforms things like sensitivity map
    
    # Prepare GM and WM partial volume maps from FAST segmentation
    if wsp.user_pv_flag:
        wsp.log.write("\nUsing user-supplied PV estimates\n")
        wsp.structural.wm_pv_asl = wsp.pvwm
        wsp.structural.gm_pv_asl = wsp.pvgm
    elif wsp.structural.struc is not None:
        wsp.structural.wm_pv_asl = reg.change_space(wsp, wsp.structural.wm_pv, "asl")
        wsp.structural.gm_pv_asl = reg.change_space(wsp, wsp.structural.gm_pv, "asl")
    else:
        raise RuntimeError("Can't do partial volume correction without either user PV estimates or structural image")

    wsp.basil_options = wsp.ifnone("basil_options", {})
    wsp.basil_options.update({"pwm" : wsp.structural.wm_pv_asl, 
                                "pgm" : wsp.structural.gm_pv_asl})
    multistep_fit.run(wsp.sub("basil_pvcorr"), prefit=False)
    wsp.quantify_wsps.append("basil_pvcorr")

