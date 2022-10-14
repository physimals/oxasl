"""
OXASL - Partial volume correction module

Copyright (c) 2008-2020 Univerisity of Oxford
"""

import numpy as np

from fsl.data.image import Image

from oxasl import basil, mask, reg
from oxasl.options import OptionCategory, OptionGroup

try:
    import oxasl_surfpvc
except ImportError:
    oxasl_surfpvc = None

class Options(OptionCategory):
    """
    Options for corrections of the input data
    """

    def __init__(self, **kwargs):
        OptionCategory.__init__(self, "corrections")

    def groups(self, parser):
        ret = []

        g = OptionGroup(parser, "Partial volume correction (PVEc)")
        g.add_option("--pvcorr", help="Apply PVEc using FAST estimates taken from --fslanat dir", action="store_true", default=False)
        g.add_option("--surf-pvcorr", help="Apply PVEc using surface PV estimates taken from --fslanat dir w/ surfaces (not mutually exclusive with --pvcorr)", action="store_true", default=False)
        g.add_option('--cores', help="Number of processor cores to use for --surf-pvcorr", type=int)
        g.add_option("--pvgm", help="GM PV estimates in ASL space (apply PVEc only, don't estimate PVs)", type="image", default=None)
        g.add_option("--pvwm", help="As above, WM PV estimates in ASL space", type="image", default=None)
        g.add_option("--pvcsf", help="As above, CSF PV estimates in ASL space", type="image", default=None)
        ret.append(g)

        return ret

def run(wsp):
    # If the user has provided manual PV maps (pvgm and pvgm) then do PVEc, even if they
    # have not explicitly given the --pvcorr option 
    user_pv_flag = ((wsp.pvwm is not None) and (wsp.pvgm is not None))
    if wsp.pvcorr or wsp.surf_pvcorr or user_pv_flag:
        # Partial volume correction is very sensitive to the mask, so recreate it
        # if it came from the structural image as this requires accurate ASL->Struc registration
        if wsp.rois.mask_src == "struc":
            wsp.rois.mask_orig = wsp.rois.mask
            wsp.rois.mask = None
            mask.run(wsp)

        if wsp.pvcorr or user_pv_flag:
            # Do partial volume correction fitting
            #
            # FIXME: We could at this point re-apply all corrections derived from structural space?
            # But would need to make sure corrections module re-transforms things like sensitivity map
            
            # Prepare GM and WM partial volume maps from FAST segmentation
            if user_pv_flag:
                wsp.log.write("\nUsing user-supplied PV estimates\n")
                wsp.structural.wm_pv_asl = wsp.pvwm
                wsp.structural.gm_pv_asl = wsp.pvgm
            else:
                wsp.structural.wm_pv_asl = reg.change_space(wsp, wsp.structural.wm_pv, "asl")
                wsp.structural.gm_pv_asl = reg.change_space(wsp, wsp.structural.gm_pv, "asl")

            wsp.basil_options = wsp.ifnone("basil_options", {})
            wsp.basil_options.update({"pwm" : wsp.structural.wm_pv_asl, 
                                      "pgm" : wsp.structural.gm_pv_asl})
            basil.run(wsp.sub("basil_pvcorr"), prefit=False)
            wsp.basildirs.append("pvcorr")

        if wsp.surf_pvcorr:
            if oxasl_surfpvc is None:
                raise RuntimeError("Surface-based PVC requested but oxasl_surfpvc is not installed")
            if user_pv_flag:
                wsp.log.write(" - WARNING: Performing surface based PVC ignores user-specified PV maps\n")
            # Prepare GM and WM partial volume maps from surface using Toblerone plugin
            # Then reform the ASL ROI mask - Toblerone does not handle the cerebellum so need
            # to mask it out
            oxasl_surfpvc.prepare_surf_pvs(wsp)
            wsp.rois.mask_pvcorr = wsp.rois.mask
            min_pv = 0.01
            new_roi = (wsp.basil_options["pwm"].data > min_pv) | (wsp.basil_options["pgm"].data > min_pv)
            wsp.rois.mask = Image(new_roi.astype(np.int8), header=wsp.rois.mask_pvcorr.header)
        
            basil.run(wsp.sub("basil_surf_pvcorr"), prefit=False)
            wsp.basildirs.append("surf_pvcorr")
