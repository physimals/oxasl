"""
OXASL - Module to generate additional ROIs for analysis

Copyright (c) 2008-2020 Univerisity of Oxford
"""
import numpy as np

from fsl.data.image import Image
from fsl.data.atlases import AtlasRegistry

from oxasl import reg

def run(wsp):
    """
    Generate additional ROIs used for reporting and region analysis

    wsp.rois.gm_asl: GM mask in ASL space (>0.5 probability)
    wsp.rois.wm_asl: WM mask in ASL space (>0.5 probability)
    wsp.rois.pure_gm_asl: 'pure' GM mask in ASL space (default > 0.8 probability, but uses pure_gm_thresh option)
    wsp.rois.pure_wm_asl: 'pure' WM mask in ASL space (default > 0.9 probability, but uses pure_gm_thresh option)
    wsp.rois.cortial_gm_asl: Cortical GM mask in ASL space, same as pure_gm_asl but masked to remove subcortical structures
    wsp.rois.cerebral_wm_asl: Cererbral WM mask in ASL space, same as pure_wm_asl but masked to remove subcortical structures
    """
    if wsp.structural.struc is not None:
        # Get the cortex 
        atlases = AtlasRegistry()
        atlases.rescanAtlases()
        atlas = atlases.loadAtlas("harvardoxford-subcortical", loadSummary=False, resolution=2)
        wsp.rois.cortex_std = Image(np.mean(atlas.data[..., 0:2] + atlas.data[..., 11:13], axis=-1) * 2, header=atlas.header)
        cortex_asl = reg.change_space(wsp, wsp.rois.cortex_std, "asl")
        wsp.rois.cortex_asl = Image((cortex_asl.data > 50).astype(int), header=cortex_asl.header)

        if wsp.structural.gm_pv_asl is None:
            wsp.structural.gm_pv_asl = reg.change_space(wsp, wsp.structural.gm_pv, "asl")
            wsp.structural.wm_pv_asl = reg.change_space(wsp, wsp.structural.wm_pv, "asl")

        gm = np.asarray(wsp.structural.gm_pv_asl.data)
        wm = np.asarray(wsp.structural.wm_pv_asl.data)

        wsp.rois.pure_gm_thresh = wsp.ifnone("pure_gm_thresh", 0.8)
        wsp.rois.pure_wm_thresh = wsp.ifnone("pure_wm_thresh", 0.0)
        some_gm, some_wm = gm > 0.5, wm > 0.5
        pure_gm, pure_wm = gm > wsp.rois.pure_gm_thresh, wm > wsp.rois.pure_wm_thresh
        wsp.rois.gm_asl = Image(some_gm.astype(int), header=wsp.structural.gm_pv_asl.header)
        wsp.rois.pure_gm_asl = Image(pure_gm.astype(int), header=wsp.structural.gm_pv_asl.header)
        wsp.rois.wm_asl = Image(some_wm.astype(int), header=wsp.structural.wm_pv_asl.header)
        wsp.rois.pure_wm_asl = Image(pure_wm.astype(int), header=wsp.structural.wm_pv_asl.header)
        wsp.rois.cortical_gm_asl = Image(np.logical_and(pure_gm, cortex_asl.data).astype(int), header=wsp.structural.gm_pv_asl.header)
        wsp.rois.cerebral_wm_asl = Image(np.logical_and(pure_wm, cortex_asl.data).astype(int), header=wsp.structural.wm_pv_asl.header)
