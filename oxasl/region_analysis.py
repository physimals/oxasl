#!/bin/env python
"""
OXASL - Module which generates perfusion stats within various ROIs
"""
import os
import sys
import argparse
import csv
import glob

import numpy as np
import scipy.stats
import pandas as pd

from fsl.data import atlases
from fsl.data.image import Image
import fsl.wrappers as fsl

from oxasl import reg
from oxasl.options import OptionCategory, IgnorableOptionGroup

# This is the probability threshold below which we do not 
# consider a voxel relevant to GM/WM averages
PVE_THRESHOLD_BASE = 0.1

class RegionAnalysisOptions(OptionCategory):
    """
    OptionGroup which contains options for region analysis
    """

    def __init__(self, title="Region analysis", **kwargs):
        OptionCategory.__init__(self, "struc", **kwargs)
        self.title = title

    def groups(self, parser):
        group = IgnorableOptionGroup(parser, self.title, ignore=self.ignore)
        group.add_option("--region-analysis", help="Perform region analysis", action="store_true", default=False)
        group.add_option("--roi-min-nvoxels", default=10, type=int,
                          help="Minimum number of relevant voxels required to report statistics")
        group.add_option("--gm-thresh", default=0.8, type=float,
                          help="Probability threshold for 'pure' grey matter")
        group.add_option("--wm-thresh", default=0.9, type=float,
                          help="Probability threshold for 'pure' white matter")
        #group.add_option("--roi-native", nargs="*", default=[],
        #                  help="Additional ROI as binarised mask in ASL space. The name of the ROI will be the stripped filename. May be specified multiple times")
        #group.add_option("--roi-struct", nargs="*", default=[],
        #                  help="Additional ROI as binarised mask in structural space. The name of the ROI will be the stripped filename. May be specified multiple times")
        #group.add_option("--roi-mni", nargs="*", default=[],
        #                  help="Additional ROI as binarised mask in MNI space. The name of the ROI will be the stripped filename. May be specified multiple times")
        group.add_option("--add-atlas-rois", action="store_true", default=False,
                          help="Add ROIs from Harvard-Oxford cortical/subcortical atlases")
        group.add_option("--save-mni-rois", action="store_true", default=False,
                          help="Save ROIs in MNI space")
        group.add_option("--save-native-rois", action="store_true", default=False,
                          help="Save ROIs in native (ASL) space")
        group.add_option("--save-native-masks", action="store_true", default=False,
                          help="Save binary masks in native (ASL) space")

        return [group, ]

def _addvar(f):
    """ Add an unused 'variance' parameter to a function which doesn't use it """
    def _f_with_var(val, var):
        return f(val)
    return _f_with_var

def mean_invvarweighted(val, var):
    """ Inverse variance weighted mean (i.e. precision weighted mean) """
    prec = 1 / var
    prec[~np.isfinite(prec)] = 0
    return np.sum(val * prec) / np.sum(prec)

def i2(val, var):
    """ I^2 Measure of heterogenaity """
    prec = 1 / var
    prec[~np.isfinite(prec)] = 0
    n = len(val)
    mu_bar = mean_invvarweighted(val, var)

    #out = []
    Q = np.sum(prec * (val - mu_bar)**2)
    #H = np.sqrt(Q/(n - 1))
    i2 = (Q-(n-1))/Q

    # Negative values map to 0 (see https://www.ncbi.nlm.nih.gov/pmc/articles/PMC192859/)
    i2 = max(i2, 0)
    #tau2_DL = (Q - n + 1)/(sum(w) - sum([x**2 for x in w])/sum(w))
    #tau2_DL = max(0, tau2_DL)
    #out.extend([mu_bar,Q,H,I2,tau2_DL])

    #out = pd.DataFrame(out)
    #out['index']=['Weighted_Mean', 'Q', 'H', 'I2', 'tau2_DL']
    #out = out.set_index('index', drop=True)

    # Return I2 as an integer percentage
    return int(100*i2+0.5)

STATS_FNS = {
    "Mean" : _addvar(np.mean),
    "Std" : _addvar(np.std),
    "Median" : _addvar(np.median),
    "IQR" : _addvar(scipy.stats.iqr),
    "Precision-weighted mean" : mean_invvarweighted,
    "I2" : i2,
}

def get_stats(stats, img, var_img, roi, suffix="", ignore_nan=True, ignore_inf=True, min_nvoxels=10, mask=None):
    """
    Get a set of statistics for a 3D image within an roi

    :param img: 3D Numpy array
    :param roi: 3D Numpy array with same dimensions as img and boolean data type
    :param ignore_nan: Voxels with care NaN in img are ignored
    :param ignore_inf: Voxels which are infinite in img are ignored
    :param min_nvoxels: If the number of voxels in the ROI is less than this number
                       (after removing Nan and infinte values) no value will be returned

    :return: Mapping from name of statistic to value. This may be NaN or infinite depending
             on the input arguments. If the number of eligible voxels is less than min_nvoxels,
             None is returned (not NaN or zero).
    """
    if list(img.shape) != list(roi.shape):
        raise ValueError("Image must have same dimensions as ROI")
    if list(var_img.shape) != list(roi.shape):
        raise ValueError("Variance image must have same dimensions as ROI")
    if mask is not None and list(mask.shape) != list(roi.shape):
        raise ValueError("Mask must have same dimensions as ROI")

    effective_roi = roi
    if ignore_nan:
        effective_roi = np.logical_and(effective_roi, ~np.isnan(img))
    if ignore_inf:
        effective_roi = np.logical_and(effective_roi, np.isfinite(img))
    if mask is not None:
        effective_roi = np.logical_and(effective_roi, mask)

    sample_data = img[effective_roi]
    sample_var = var_img[effective_roi]
    nvoxels = len(sample_data)
    stats["Nvoxels" + suffix] = nvoxels
    for stat, fn in STATS_FNS.items():
        if nvoxels < min_nvoxels:
            stats[stat + suffix] = None
        else:
            stats[stat + suffix] = fn(sample_data, sample_var)

def oxasl_add_roi(wsp, rois, name, roi_native, threshold, roi_struct=None, roi_mni=None):
    wsp.log.write(" - %s..." % name)
    rois.append({
        "name" : name,
        "roi_native" : roi_native,
        "mask_native" : Image((roi_native.data > threshold).astype(np.int32), header=roi_native.header),
        "roi_struct" : roi_struct,
        "roi_mni" : roi_mni
    })
    wsp.log.write("DONE\n")

def oxasl_add_atlas(wsp, rois, atlas_name, resolution=2, threshold=0.5):
    """
    Get ROIs from an FSL atlas
    
    :param rois: Mapping from name to ROI array which will be updated
    :param mni2struc_warp: Warp image containing MNI->structural space warp
    :param struct2asl_mat: Matrix for struct->ASL transformation
    :param atlas_name: Name of the FSL atlas
    :param resolution: Resolution in mm
    :param threshold: Threshold for probabilistic atlases
    """
    wsp.log.write("\nAdding ROIs from standard atlas: %s (resolution=%imm, thresholding at %.2f)\n" % (atlas_name, resolution, threshold))
    reg.reg_struc2std(wsp, fnirt=True)
    registry = atlases.registry
    registry.rescanAtlases()
    desc = registry.getAtlasDescription(atlas_name)
    atlas = registry.loadAtlas(desc.atlasID, resolution=2)
    for label in desc.labels[:3]:
        roi_mni = atlas.get(label=label)
        roi_native = reg.std2asl(wsp, roi_mni)
        oxasl_add_roi(wsp, rois, label.name, roi_native, threshold=50, roi_mni=roi_mni)

def oxasl_perfusion_data(wsp):
    perfusion_data = [
        {
            "suffix" : "", 
            "f" : wsp.output.native.perfusion_calib,
            "var" :  wsp.output.native.perfusion_var_calib,
            "mask" : wsp.rois.mask.data,
        },
    ]
    if wsp.pvcorr:
        wsp.log.write(" - Found partial volume corrected results - will mask ROIs using 'base' GM/WM masks (threshold: %.2f)\n" % PVE_THRESHOLD_BASE)
        perfusion_data.extend([
            {
                "suffix" : "_gm", 
                "f" : wsp.output_pvcorr.native.perfusion_calib,
                "var" : wsp.output_pvcorr.native.perfusion_var_calib,
                "mask" : np.logical_and(wsp.rois.mask.data, wsp.structural.gm_pv_asl.data > PVE_THRESHOLD_BASE),
            },
            {
                "suffix" : "_wm", 
                "f" : wsp.output_pvcorr.native.pvcorr.perfusion_wm_calib,
                "var" : wsp.output_pvcorr.native.perfusion_wm_var_calib,
                "mask" : np.logical_and(wsp.rois.mask.data, wsp.structural.wm_pv_asl.data > PVE_THRESHOLD_BASE),
            },
        ])
    else:
        wsp.log.write(" - No partial volume corrected results - will mask ROIs using 'pure' GM/WM masks (PVE thresholds: %.2f / %.2f)\n" % (wsp.gm_thresh, wsp.wm_thresh))
        perfusion_data.extend([
            {
                "suffix" : "_gm",
                "f" : wsp.output.native.perfusion_calib,
                "var" :  wsp.output.native.perfusion_var_calib,
                "mask" : np.logical_and(wsp.rois.mask.data, wsp.structural.gm_pv_asl.data > wsp.gm_thresh),
            },
            {
                "suffix" : "_wm",
                "f" : wsp.output.native.perfusion_calib,
                "var" :  wsp.output.native.perfusion_var_calib,
                "mask" : np.logical_and(wsp.rois.mask.data, wsp.structural.wm_pv_asl.data > wsp.wm_thresh),
            },
        ])
    return perfusion_data

def run(wsp):
    """ Entry point for OXASL """
    if wsp.pvwm is not None:
        wsp.structural.wm_pv_asl = wsp.pvwm
    else:
        wsp.structural.wm_pv_asl = reg.struc2asl(wsp, wsp.structural.wm_pv)

    if wsp.pvwm is not None:
        wsp.structural.gm_pv_asl = wsp.pvgm
    else:
        wsp.structural.gm_pv_asl = reg.struc2asl(wsp, wsp.structural.gm_pv)
    
    wsp.gm_thresh, wsp.wm_thresh = wsp.ifnone("gm_thresh", 0.8), wsp.ifnone("wm_thresh", 0.9)

    wsp.log.write("\nRegionwise analysis\n")

    wsp.log.write("\nLoading perfusion images\n")
    perfusion_data = oxasl_perfusion_data(wsp)
    
    rois = []
    wsp.log.write("\nLoading generic ROIs\n")
    oxasl_add_roi(wsp, rois, "10%+GM", wsp.structural.gm_pv_asl, 0.1)
    oxasl_add_roi(wsp, rois, "10%+WM", wsp.structural.wm_pv_asl, 0.1)
    oxasl_add_roi(wsp, rois, "%i%%+GM" % (wsp.gm_thresh*100), wsp.structural.gm_pv_asl, wsp.gm_thresh)
    oxasl_add_roi(wsp, rois, "%i%%+WM" % (wsp.wm_thresh*100), wsp.structural.wm_pv_asl, wsp.wm_thresh)
    
    # Add ROIs from standard atlases
    oxasl_add_atlas(wsp, rois, "harvardoxford-cortical", threshold=0.5)
    oxasl_add_atlas(wsp, rois, "harvardoxford-subcortical", threshold=0.5)

    # Get stats in each ROI. Add name to stats dict to make TSV output easier
    wsp.log.write("\nGetting stats - minimum of %i voxels to report in region\n" % wsp.roi_min_nvoxels)
    for item in perfusion_data:
        stats = []
        for roi in rois:
            roi_stats = {"name" : roi["name"]}
            get_stats(roi_stats, item["f"].data, item["var"].data, roi["mask_native"].data, mask=item["mask"], min_nvoxels=wsp.roi_min_nvoxels)
            stats.append(roi_stats)
            columns = list(roi_stats.keys())

        setattr(wsp, "roi_stats%s" % item["suffix"], pd.DataFrame(stats, columns=columns))

    # Save output masks/PVE maps
    for roi in rois:
        fname = roi["name"].replace(" ", "_").replace(",", "").lower()
        if wsp.save_native_rois and "roi_native" in roi:
            setattr(wsp, fname, roi["roi_native"])
        if wsp.save_native_masks and "mask_native" in roi:
            setattr(wsp, fname, roi["mask_native"])
        if wsp.save_mni_rois and "roi_mni" in roi:
            setattr(wsp, fname, roi["roi_mni"])

    wsp.log.write("\nDONE\n")
