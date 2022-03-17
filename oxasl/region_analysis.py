#!/bin/env python
"""
OXASL - Module which generates perfusion stats within various ROIs

Copyright (c) 2008-2020 Univerisity of Oxford
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
from oxasl.options import OptionCategory, OptionGroup

class Options(OptionCategory):
    """
    Options for region analysis
    """

    def __init__(self):
        OptionCategory.__init__(self, "region_analysis")

    def groups(self, parser):
        group = OptionGroup(parser, "Region analysis")
        group.add_option("--region-analysis", help="Perform region analysis", action="store_true", default=False)
        group.add_option("--roi-min-nvoxels", default=10, type=int,
                          help="Minimum number of relevant voxels required to report statistics")
        group.add_option("--gm-thresh", default=0.8, type=float,
                          help="Probability threshold for 'pure' grey matter")
        group.add_option("--wm-thresh", default=0.9, type=float,
                          help="Probability threshold for 'pure' white matter")
        group.add_option("--min-gm-thresh", default=0.1, type=float,
                          help="Probability threshold for a voxel to be included in GM stats")
        group.add_option("--min-wm-thresh", default=0.1, type=float,
                          help="Probability threshold for a voxel to be included in WM stats")
        group.add_option("--add-roi",
                          help="Additional ROI as binarised mask in ASL, structural, or MNI space. The name of the ROI will be the stripped filename. May be specified multiple times")
        group.add_option("--add-mni-atlas",
                          help="Additional regions as labelled image in MNI space. If a single region is contained, the name of the ROI will be the stripped filename, otherwise use --add-mni-atlas-names. May be specified multiple times (comma separated)")
        group.add_option("--add-mni-atlas-labels",
                          help="Comma separated filenames containing names of regions in atlas given in --add-mni-atlas")
        #group.add_option("--add-standard-atlases", action="store_true", default=False,
        #                  help="Add ROIs from Harvard-Oxford cortical/subcortical atlases")
        group.add_option("--save-mni-rois", action="store_true", default=False,
                          help="Save ROIs in MNI space")
        group.add_option("--save-struct-rois", action="store_true", default=False,
                          help="Save ROIs in structural space")
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
    if Q == 0:
        i2 = 0
    else:
        i2 = (Q-(n-1))/Q
    #H = np.sqrt(Q/(n - 1))

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
    # Variance should not be zero but sometimes is - maybe masking?
    sample_var[sample_var == 0] = 1e-6
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
        "mask_native" : Image((roi_native.data > threshold).astype(np.int), header=roi_native.header),
        "roi_struct" : roi_struct,
        "roi_mni" : roi_mni
    })
    wsp.log.write("DONE\n")

def oxasl_add_fsl_atlas(wsp, rois, atlas_name, resolution=2, threshold=0.5):
    """
    Get ROIs from an FSL atlas
    
    :param rois: Mapping from name to ROI array which will be updated
    :param atlas_name: Name of the FSL atlas
    :param resolution: Resolution in mm
    :param threshold: Threshold for probabilistic atlases
    """
    wsp.log.write("\nAdding ROIs from standard atlas: %s (resolution=%imm, thresholding at %.2f)\n" % (atlas_name, resolution, threshold))
    registry = atlases.registry
    registry.rescanAtlases()
    desc = registry.getAtlasDescription(atlas_name)
    atlas = registry.loadAtlas(desc.atlasID, resolution=2)
    for label in desc.labels:
        roi_mni = atlas.get(label=label)
        roi_native = reg.change_space(wsp, roi_mni, "native")
        oxasl_add_roi(wsp, rois, label.name, roi_native, threshold=50, roi_mni=roi_mni)

def oxasl_add_custom_atlas(wsp, rois, atlas_img, region_names):
    """
    Get ROIs from an atlas described by a label image in MNI space
 
    :param rois: Mapping from name to ROI array which will be updated
    :param atlas_img: Atlas label image
    :param region_names: Atlas label image region names
    """
    wsp.log.write("\nAdding ROIs from MNI atlas image: %s\n" % (atlas_img.name))
    labels = [idx for idx in np.unique(atlas_img.data) if idx != 0]
    if len(labels) != len(region_names):
        region_names = ["Region %i" % label for label in labels]
    for idx, label in enumerate(labels):
        roi_mni = atlas_img.data.copy()
        roi_mni[roi_mni != label] = 0
        name = region_names[idx]
        roi_native = reg.change_space(wsp, roi_mni, "native")
        oxasl_add_roi(wsp, rois, name, roi_native, threshold=0.5, roi_mni=roi_mni)

def oxasl_perfusion_data(wsp):
    perfusion_data = [
        {
            "suffix" : "", 
            "f" : wsp.perfusion,
            "var" :  wsp.perfusion_var,
            "mask" : wsp.rois.mask.data,
        },
    ]
    if wsp.perfusion_wm is not None:
        wsp.log.write(" - Found partial volume corrected results - will mask ROIs using 'base' GM/WM masks (PVE thresholds: %.2f / %.2f)\n" % (wsp.min_gm_thres, wsp.min_wm_thresh))
        perfusion_data.extend([
            {
                "suffix" : "_gm", 
                "f" : wsp.perfusion,
                "var" : wsp.perfusion_var,
                "mask" : np.logical_and(wsp.rois.mask.data, wsp.structural.gm_pv_asl.data > wsp.min_gm_thresh),
            },
            {
                "suffix" : "_wm", 
                "f" : wsp.perfusion_wm,
                "var" : wsp.perfusion_wm_var,
                "mask" : np.logical_and(wsp.rois.mask.data, wsp.structural.wm_pv_asl.data > wsp.min_wm_thresh),
            },
        ])
    else:
        wsp.log.write(" - No partial volume corrected results - will mask ROIs using 'pure' GM/WM masks (PVE thresholds: %.2f / %.2f)\n" % (wsp.gm_thresh, wsp.wm_thresh))
        perfusion_data.extend([
            {
                "suffix" : "_gm",
                "f" : wsp.perfusion,
                "var" :  wsp.perfusion_var,
                "mask" : np.logical_and(wsp.rois.mask.data, wsp.structural.gm_pv_asl.data > wsp.gm_thresh),
            },
            {
                "suffix" : "_wm",
                "f" : wsp.perfusion,
                "var" :  wsp.perfusion_var,
                "mask" : np.logical_and(wsp.rois.mask.data, wsp.structural.wm_pv_asl.data > wsp.wm_thresh),
            },
        ])
    return perfusion_data

def run(wsp):
    """
    Entry point for OXASL

    wsp: Output data workspace
    """
    if not wsp.region_analysis:
        return

    wsp.log.write("\nRegionwise analysis\n")

    if wsp.pvwm is not None:
        wsp.structural.wm_pv_asl = wsp.pvwm
    else:
        wsp.structural.wm_pv_asl = reg.change_space(wsp, wsp.structural.wm_pv, "native")

    if wsp.pvwm is not None:
        wsp.structural.gm_pv_asl = wsp.pvgm
    else:
        wsp.structural.gm_pv_asl = reg.change_space(wsp, wsp.structural.gm_pv, "native")
    
    wsp.gm_thresh, wsp.wm_thresh = wsp.ifnone("gm_thresh", 0.8), wsp.ifnone("wm_thresh", 0.9)
    wsp.min_gm_thresh, wsp.min_wm_thresh = wsp.ifnone("min_gm_thresh", 0.1), wsp.ifnone("min_wm_thresh", 0.1)

    rois = []
    wsp.log.write("\nLoading generic ROIs\n")
    oxasl_add_roi(wsp, rois, "%i%%+GM" % (wsp.min_gm_thresh*100), wsp.structural.gm_pv_asl, wsp.min_gm_thresh)
    oxasl_add_roi(wsp, rois, "%i%%+WM" % (wsp.min_wm_thresh*100), wsp.structural.wm_pv_asl, wsp.min_wm_thresh)
    oxasl_add_roi(wsp, rois, "%i%%+GM" % (wsp.gm_thresh*100), wsp.structural.gm_pv_asl, wsp.gm_thresh)
    oxasl_add_roi(wsp, rois, "%i%%+WM" % (wsp.wm_thresh*100), wsp.structural.wm_pv_asl, wsp.wm_thresh)

    # Add ROIs from command line
    print("\nLoading user-specified ROIs")
    add_roi = [l.strip() for l in wsp.ifnone("add_roi", "").split(",") if l.strip() != ""]

    for fname in add_roi:
        roi_native = reg.change_space(wsp, Image(fname), "native")
        oxasl_add_roi(wsp, rois, os.path.basename(fname).split(".")[0], roi_native, 0.5)
    
    add_atlas = [l.strip() for l in wsp.ifnone("add_mni_atlas", "").split(",") if l.strip() != ""]
    atlas_labels = [l.strip() for l in wsp.ifnone("add_mni_atlas_labels", "").split(",") if l.strip() != ""]
    for idx, fname in enumerate(add_atlas):
        if idx < len(atlas_labels):
            with open(atlas_labels[idx]) as f:
                names = [l.strip() for l in f.readlines()]
        else:
            names = [os.path.basename(fname).split(".")[0],]
        oxasl_add_custom_atlas(wsp, rois, Image(fname), names)

    # Add ROIs from standard atlases
    oxasl_add_fsl_atlas(wsp, rois, "harvardoxford-cortical", threshold=0.5)
    oxasl_add_fsl_atlas(wsp, rois, "harvardoxford-subcortical", threshold=0.5)

    # Save output masks/PVE maps
    for roi in rois:
        fname = roi["name"].replace(" ", "_").replace(",", "").lower()
        if wsp.save_native_rois and "roi_native" in roi:
            setattr(wsp, fname, roi["roi_native"])
        if wsp.save_native_masks and "mask_native" in roi:
            setattr(wsp, fname, roi["mask_native"])
        if wsp.save_struct_rois and "roi_struct" in roi:
            setattr(wsp, fname, roi["roi_struct"])
        if wsp.save_mni_rois and "roi_mni" in roi:
            setattr(wsp, fname, roi["roi_mni"])

    for calib_method in wsp.calibration.calib_method:
        wsp.log.write("\nCalibration method: %s\n" % calib_method)
        calib_wsp = getattr(wsp.native, "calib_%s" % calib_method)

        wsp.log.write("\nLoading perfusion images\n")
        perfusion_data = oxasl_perfusion_data(calib_wsp)

        # Get stats in each ROI. Add name to stats dict to make TSV output easier
        wsp.log.write("\nGetting stats - minimum of %i voxels to report in region\n" % wsp.roi_min_nvoxels)
        for item in perfusion_data:
            stats = []
            for roi in rois:
                roi_stats = {"name" : roi["name"]}
                get_stats(roi_stats, item["f"].data, item["var"].data, roi["mask_native"].data, mask=item["mask"], min_nvoxels=wsp.roi_min_nvoxels)
                stats.append(roi_stats)
                columns = list(roi_stats.keys())

            setattr(calib_wsp, "roi_stats%s" % item["suffix"], pd.DataFrame(stats, columns=columns))

    wsp.log.write("\nDONE\n")
