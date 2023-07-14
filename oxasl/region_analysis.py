#!/bin/env python
"""
OXASL - Module which generates perfusion stats within various ROIs

Copyright (c) 2008-2020 Univerisity of Oxford
"""
import os

import numpy as np
import scipy.stats
import pandas as pd
from scipy.fftpack import fft, ifft

from fsl.data import atlases
from fsl.data.image import Image

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
        group.add_option("--pure-gm-thresh", "--gm-thresh", default=0.8, type=float,
                          help="Probability threshold for 'pure' grey matter")
        group.add_option("--pure-wm-thresh", "--wm-thresh", default=0.9, type=float,
                          help="Probability threshold for 'pure' white matter")
        group.add_option("--min-gm-thresh", default=0.1, type=float,
                          help="Probability threshold for a voxel to be included in GM stats")
        group.add_option("--min-wm-thresh", default=0.1, type=float,
                          help="Probability threshold for a voxel to be included in WM stats")
        group.add_option("--att-stats", help="Generate statistics for arterial transit time", action="store_true", default=False)
        group.add_option("--fuzzy-sets", help="ROI sets (groups of ROIs with total PV=1) will be modelled as fuzzy rather than binarized", action="store_true", default=False)
        group.add_option("--psf", help="Point-spread function for ROI sets. If specified, PSF will be applied to ROI sets (in ASL space) and 'fuzzy' mean used")
        group.add_option("--add-roi",
                          help="Additional ROI as binarised mask in ASL, structural, or MNI space. The name of the ROI will be the stripped filename. May be specified multiple times")
        group.add_option("--add-atlas", "--add-std-atlas", "--add-mni-atlas",
                          help="Additional regions as labelled image. If a single region is contained, the name of the ROI will be the stripped filename, otherwise use --add-std-atlas-names. May be specified multiple times (comma separated)")
        group.add_option("--add-atlas-labels", "--add-std-atlas-labels", "--add-mni-atlas-labels",
                          help="Comma separated filenames containing names of regions in atlas given in --add-atlas")
        #group.add_option("--add-standard-atlases", action="store_true", default=False,
        #                  help="Add ROIs from Harvard-Oxford cortical/subcortical atlases")
        group.add_option("--save-std-rois", "--save-mni-rois", action="store_true", default=False,
                          help="Save ROIs in standard (MNI) space")
        group.add_option("--save-struct-rois", action="store_true", default=False,
                          help="Save ROIs in structural space")
        group.add_option("--save-asl-rois", "--save-native-rois", action="store_true", default=False,
                          help="Save ROIs in ASL space")
        group.add_option("--save-asl-masks", "--save-native-masks", action="store_true", default=False,
                          help="Save binary masks in ASL space")

        return [group, ]

def _addvar(f):
    """ Add an unused 'variance' parameter to a function which doesn't use it """
    def _f_with_var(val, var):
        return f(val)
    return _f_with_var

def mean_invvarweighted(val, var):
    """ Inverse variance weighted mean (i.e. precision weighted mean) """
    if var is None:
        return None
    prec = 1 / var
    prec[~np.isfinite(prec)] = 0
    return np.sum(val * prec) / np.sum(prec)

def i2(val, var):
    """ I^2 Measure of heterogenaity """
    if var is None:
        return None
    prec = 1 / var
    prec[~np.isfinite(prec)] = 0
    n = len(val)
    mu_bar = mean_invvarweighted(val, var)

    Q = np.sum(prec * (val - mu_bar)**2)
    if Q == 0:
        i2 = 0
    else:
        i2 = (Q-(n-1))/Q

    # Negative values map to 0 (see https://www.ncbi.nlm.nih.gov/pmc/articles/PMC192859/)
    i2 = max(i2, 0)

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

def get_stats_binary(wsp, stats, img, var_img, roi, suffix="",
                     ignore_nan=True, ignore_inf=True, ignore_zerovar=True,
                     min_nvoxels=10, mask=None):
    """
    Get a set of statistics for a 3D image within a binary roi

    :param img: 3D Numpy array
    :param roi: 3D Numpy array with same dimensions as img and boolean data type
    :param ignore_nan: Voxels with care NaN in img are ignored
    :param ignore_inf: Voxels which are infinite in img are ignored
    :param ignore_zerovar: Ignore voxels with zero variance
    :param min_nvoxels: If the number of voxels in the ROI is less than this number
                       (after removing Nan and infinte values) no value will be returned
    :param mask: Additional binary mask to apply to ROI

    :return: Mapping from name of statistic to value. This may be NaN or infinite depending
             on the input arguments. If the number of eligible voxels is less than min_nvoxels,
             None is returned (not NaN or zero).
    """
    if list(img.shape) != list(roi.shape):
        raise ValueError("Image must have same dimensions as ROI")
    if var_img is not None and list(var_img.shape) != list(roi.shape):
        raise ValueError("Variance image must have same dimensions as ROI")
    if mask is not None and list(mask.shape) != list(roi.shape):
        raise ValueError("Mask must have same dimensions as ROI")

    if mask is None:
        mask = np.ones(roi.shape, dtype=int)
    if ignore_nan:
        mask = np.logical_and(mask, ~np.isnan(img))
    if ignore_inf:
        mask = np.logical_and(mask, np.isfinite(img))
    if var_img is not None and ignore_zerovar:
        mask = np.logical_and(mask, var_img > 0)

    effective_roi = np.logical_and(roi, mask)

    sample_data = img[effective_roi]
    if var_img is not None:
        sample_var = var_img[effective_roi]
        sample_var[sample_var == 0] = 1e-6
    else:
        sample_var = None
    nvoxels = len(sample_data)
    stats["Nvoxels" + suffix] = nvoxels
    for stat, fn in STATS_FNS.items():
        if nvoxels < min_nvoxels:
            stats[stat + suffix] = None
        else:
            stats[stat + suffix] = fn(sample_data, sample_var)

def standardise_fuzzy(wsp, roi_set, mode="expand"):
    """
    Ensure that fuzzy ROI set has PVs that sum to 1 in every voxel

    :param roi_set: 4D Numpy array with same dimensions as img and each volume
                    containing partial volumes for each ROI in the set
    :param mode: 'expand' adds an ROI to the set containing any non-included PV.
                 'normalise' scales PVs so they sum to 1
    """
    if mode not in ("normalise", "expand"):
        raise ValueError(f"Mode {mode} not supported - expected `normalise` or `expand`")

    # Warning if normalising and ROI only has one label
    if np.shape(roi_set)[1] == 1 and mode == "normalise":
        wsp.log.write("WARNING: Only one ROI in set with mode='normalise' - all nonzero voxels will have PV=1\n")

    # Get total PV in each voxel
    roi_sum = np.sum(roi_set, axis=1, keepdims=True)

    if mode == "normalise":
        roi_set = np.where(roi_sum!=0, roi_set/roi_sum, 0)
    elif mode == "expand":
        # TODO: insert a check to see if any diff_img < 0 and raise a warning
        roi_set = np.concatenate((roi_set, 1 - roi_sum), axis=-1)
    return roi_set

def get_stats_fuzzy(wsp, stats, img, var_img, roi_set, suffix="",
                    ignore_nan=True, ignore_inf=True, ignore_zerovar=True,
                    mask=None, pv_threshold=0.):
    """
    Get a set of statistics for a set of 'fuzzy' ROIs. Unless otherwise
    specified, parameters are as for ``get_stats_binary``

    :param roi_set: 4D Numpy array with same dimensions as img and each volume
                    containing partial volumes for each ROI in the set
    :param pv_threshold: Minimum total PV for voxel to be included

    :return: Mapping from name of statistic to sequence of values, one for each ROI in the set.
             This may be NaN or infinite depending on the input arguments. 
    """
    roi_shape = list(roi_set.shape)[:3]
    if list(img.shape) != roi_shape:
        raise ValueError("Image must have same dimensions as ROI")
    if var_img is not None and list(var_img.shape) != roi_shape:
        raise ValueError("Variance image must have same dimensions as ROI")
    if mask is not None and list(mask.shape) != roi_shape:
        raise ValueError("Mask must have same dimensions as ROI")

    if mask is None:
        mask = np.ones(roi_shape, dtype=np.int32)
    if ignore_nan:
        mask = np.logical_and(mask, ~np.isnan(img))
    if ignore_inf:
        mask = np.logical_and(mask, np.isfinite(img))
    if var_img is not None and ignore_zerovar:
        mask = np.logical_and(mask, var_img > 0)

    # Only take voxels where at least one of the ROIs has non-zero percentage
    mask = np.logical_and(mask, np.sum(roi_set, axis=3) > pv_threshold)

    # Flatten ROI PVs and data into masked 2D array
    roi_array = roi_set[mask]
    g = img[mask]

    # Standardize ROI set so total PV is 1
    roi_array = standardise_fuzzy(wsp, roi_array, mode='expand')

    # Ask Jack about this???
    #if var:
    #    roi_array = np.square(roi_array)

    HT = roi_array.T
    wsp.log.write(f" - Fuzzy ROI set: condition number for transfer matrix (unweighted) = {np.linalg.cond(HT):.2f}\n")

    # Calculate roi means by linear regression
    means_lstsq, _res, _rank, _s = np.linalg.lstsq(HT@roi_array, HT@g[..., np.newaxis],
                                                   rcond=None) # None uses future default
                                                               # and silences warning

    # Note that we do not report stats for the 'background' ROI added to ensure total PV of 1
    stats["Nvoxels" + suffix] = [np.count_nonzero(roi_array[:, idx] > pv_threshold) for idx in range(roi_set.shape[-1])]
    stats["Mean" + suffix] = np.atleast_1d(np.squeeze(means_lstsq[:-1]))

    # If variance has been supplied add a precision-weighted mean
    if var_img is not None:
        V_inv = scipy.sparse.diags(1/var_img[mask])
        HT = roi_array.T @ V_inv
        wsp.log.write(f" - Fuzzy ROI set: condition number for transfer matrix (prec-weighted) = {np.linalg.cond(HT):.2f}\n")

        # Calculate roi means by linear regression
        means_lstsq, _res, _rank, _s = np.linalg.lstsq(HT@roi_array, HT@g[..., np.newaxis],
                                                       rcond=None) # None uses future default
                                                                   # and silences warning
        stats["Precision-weighted mean" + suffix] = np.atleast_1d(np.squeeze(means_lstsq[:-1]))

def apply_psf(array, psf):
    """
    Apply PSF blurring to an array (typically a binary or PV mask)

    The PSF is assumed to act only along the Z axis
    """
    if psf is None:
        return array

    # Make sure array is 4D
    array = array.astype(np.float32)
    was_3d = False
    if array.ndim == 3:
        was_3d = True
        array = array[..., np.newaxis]

    # Detect Z dimension padding by comparing size of psf and data
    n_slices = psf.shape[0]
    padding_slices = n_slices - array.shape[2]
    if padding_slices < 0 or padding_slices % 2 != 0:
        raise ValueError("Invalid padding in psf: %i slices vs %i in data (difference must be even and > 0)" % (n_slices, array.shape[2]))
    padding_slices = int(padding_slices/2)
    array = np.pad(array, [(0, 0), (0, 0), (padding_slices, padding_slices), (0, 0)], 'edge')

    # Calculate mean along z direction for each (x, y, t) to demean volume
    zmean = np.expand_dims(np.mean(array, 2), 2)
    array = array - zmean

    # Apply blurring using multiplication in Fourier domain
    fftkern = fft(psf)[np.newaxis, np.newaxis, ..., np.newaxis]
    fftvol = fft(array, axis=2)

    # Get blurred volume in Image domain and add DC term back
    blurred_array = np.real(ifft(np.multiply(fftvol, fftkern), axis=2)) + zmean

    # Unpad and Unsqueeze extra dimension if original volume was 3D
    blurred_array = blurred_array[:, :, padding_slices:-padding_slices, :]
    if was_3d:
        blurred_array = np.squeeze(blurred_array, axis=3)
    return blurred_array

def get_stats(wsp, roi_stats, roi, data_item):
    f = data_item["f"].data
    var = data_item.get("var", None)
    if var is not None:
        var = var.data

    if "fuzzy_asl" in roi:
        get_stats_fuzzy(wsp, roi_stats, f, var, roi["fuzzy_asl"], mask=data_item["mask"])
    elif "mask_asl" in roi:
        get_stats_binary(wsp, roi_stats, f, var, roi["mask_asl"], mask=data_item["mask"], min_nvoxels=wsp.roi_min_nvoxels)
    else:
        # Should never happen
        raise RuntimeError("No ASL-space ROI to get stats: %s" % str(roi))

def add_roi(wsp, rois, name, roi, threshold=0.5):
    """
    Add an ROI

    :param rois: Current list of ROIs
    :param name: Name for ROI
    :param threshold: Threshold for generating binary mask
    """
    wsp.log.write(" - %s..." % name)
    roi_space = reg.get_img_space(wsp, roi)
    roi_asl = reg.change_space(wsp, roi, "asl")
    rois.append({
        "name" : name,
        "roi_asl" : roi_asl,
        "mask_asl" : (roi_asl.data > threshold).astype(np.int32),
        "roi_%s" % roi_space : roi,
    })
    wsp.log.write("DONE\n")

def add_roi_set(wsp, rois, names, roi_set, threshold=None):
    """
    Add an ROI set

    :param rois: Current list of ROIs
    :param names: Array of ROI names
    :param roi_set: 4D Image where volumes define disjoint masks
    :param threshold: Optional threshold for generating binary native space mask.
                      Note this is not required since roi sets can be 'fuzzy'
    """
    wsp.log.write(" - %s..." % ",".join(names))
    roi_set_space = reg.get_img_space(wsp, roi_set)
    roi_set_asl = reg.change_space(wsp, roi_set, "asl")
    if threshold:
        mask_set_asl = roi_set_asl.data > threshold
    else:
        mask_set_asl = roi_set_asl.data
    fuzzy_set_asl = apply_psf(mask_set_asl, wsp.psf)
    rois.append({
        "names" : names, 
        "roi_%s" % roi_set_space : roi_set,
        "roi_asl" : roi_set_asl, 
        "mask_asl" : mask_set_asl, 
        "fuzzy_asl" : fuzzy_set_asl
    })
    wsp.log.write("DONE\n")

def add_roi_set_from_fsl_atlas(wsp, rois, atlas_name, resolution=2, threshold=0.5):
    """
    Get ROIs from an FSL atlas
    
    :param rois: Current list of ROIs
    :param atlas_name: Name of the FSL atlas
    :param resolution: Resolution in mm
    :param threshold: Threshold for probabilistic atlases
    """
    wsp.log.write("\nAdding ROI set from standard atlas: %s (resolution=%imm, thresholding at %.2f)\n" % (atlas_name, resolution, threshold))
    registry = atlases.registry
    registry.rescanAtlases()
    desc = registry.getAtlasDescription(atlas_name)
    atlas = registry.loadAtlas(desc.atlasID, resolution=2)
    roi_set, names = [], []
    for label in desc.labels:
        roi_region = atlas.get(label=label)
        # Convert to probability
        roi_region = Image(roi_region.data / 100.0, header=roi_region.header)
        if not wsp.fuzzy_sets and wsp.psf is None:
            add_roi(wsp, rois, label.name, roi_region, threshold=threshold)
        else:
            roi_set.append(roi_region.data)
            names.append(label.name)

    if wsp.fuzzy_sets or wsp.psf is not None:
        # When treating as an ROI set, do not threshold, allow it to be fuzzy
        roi_set = Image(np.stack(roi_set, axis=3), header=roi_region.header)
        add_roi_set(rois, names, roi_set, names)

def add_rois_from_3d_label_atlas(wsp, rois, atlas_img, region_names):
    """
    Get ROIs from an atlas described by a 3D label image
 
    This is a 3D integer image where each unizue nonzero voxel value defines an
    ROI region

    :param rois: Current list of ROIs
    :param atlas_img: Atlas label image
    :param region_names: Atlas label image region names
    """
    wsp.log.write("\nAdding ROIs from 3D atlas label image: %s\n" % (atlas_img.name))
    if atlas_img.data.ndim != 3:
        raise RuntimeError("Label atlas is not 3D: %s" % atlas_img.name)
    labels = [idx for idx in np.unique(atlas_img.data) if idx != 0]
    if len(labels) != len(region_names):
        region_names = ["Region %i" % label for label in labels]

    for name, label in zip(region_names, labels):
        roi_data = atlas_img.data.copy()
        roi_bin = (roi_data == label).astype(np.int32)
        roi = Image(roi_bin, header=atlas_img.header)
        add_roi(wsp, rois, name, roi)

def add_roi_set_from_4d_atlas(wsp, rois, atlas_img, region_names, threshold=0.5):
    """
    Get ROIs from an atlas described by a 4D image
 
    Each volume defines an ROI, either as binary masks or a set of fuzzy ROIs

    :param rois: Current list of ROIs
    :param atlas_img: 4D Atlas ROI image
    :param region_names: Atlas label image region names
    :param threshold: Threshold for binarizing volumes
    """
    wsp.log.write("\nAdding ROI set from 4D atlas: %s\n" % (atlas_img.name))
    if atlas_img.data.ndim != 4:
        raise RuntimeError("Atlas is not 4D: %s" % atlas_img.name)
    labels = range(atlas_img.shape[-1])
    if len(labels) != len(region_names):
        region_names = ["Region %i" % label for label in labels]

    roi_set, names = [], []
    for name, label in zip(region_names, labels):
        roi_region = atlas_img.data[..., label]
        if not wsp.fuzzy_sets and wsp.psf is None:
            add_roi(rois, name, roi_region, threshold=threshold)
        else:
            roi_set.append(roi_region.data / 100)
            names.append(label.name)

    if wsp.fuzzy_sets or wsp.psf is not None:
        roi_set = Image(np.stack(roi_set, axis=3), header=roi_region.header)
        add_roi_set(rois, roi_set, names)

def get_perfusion_data(wsp):
    if wsp.perfusion_wm is not None:
        wsp.log.write(" - Found partial volume corrected results - will mask ROIs using 'base' GM/WM masks (PVE thresholds: %.2f / %.2f)\n" % (wsp.min_gm_thresh, wsp.min_wm_thresh))
        data = [
            {
                "suffix" : "_gm", 
                "f" : wsp.perfusion,
                "var" : wsp.perfusion_var,
                "mask" : np.logical_and(wsp.mask.data, wsp.structural.gm_pv_asl.data > wsp.min_gm_thresh),
            },
            {
                "suffix" : "_wm", 
                "f" : wsp.perfusion_wm,
                "var" : wsp.perfusion_wm_var,
                "mask" : np.logical_and(wsp.mask.data, wsp.structural.wm_pv_asl.data > wsp.min_wm_thresh),
            },
        ]
    else:
        wsp.log.write(" - No partial volume corrected results - will mask ROIs using 'pure' GM/WM masks (PVE thresholds: %.2f / %.2f)\n" % (wsp.pure_gm_thresh, wsp.pure_wm_thresh))
        data = [
            {
                "suffix" : "", 
                "f" : wsp.perfusion,
                "var" :  wsp.perfusion_var,
                "mask" : wsp.mask.data,
            },
            {
                "suffix" : "_gm",
                "f" : wsp.perfusion,
                "var" :  wsp.perfusion_var,
                "mask" : np.logical_and(wsp.mask.data, wsp.structural.gm_pv_asl.data > wsp.pure_gm_thresh),
            },
            {
                "suffix" : "_wm",
                "f" : wsp.perfusion,
                "var" :  wsp.perfusion_var,
                "mask" : np.logical_and(wsp.mask.data, wsp.structural.wm_pv_asl.data > wsp.pure_wm_thresh),
            },
        ]
    return data

def get_arrival_data(wsp):
    if wsp.arrival is None:
        return []

    # Note that for the arrival mask we also remove voxels which will be
    # eliminated from the perfusion stats because of nan/inf values or zero
    # variances. This is a bit of a hack but should help ensure that the
    # voxel set is consistent between the two measures. Perfusion can have
    # invalid values in voxels where the arrival time is valid because of
    # voxelwise calibration, however we do not expect the reverse to occur
    if wsp.perfusion_var is None:
        f_good = np.isfinite(wsp.perfusion.data)
    else:
        f_good = np.logical_and(np.isfinite(wsp.perfusion.data), wsp.perfusion_var.data > 0)
    effective_mask = np.logical_and(wsp.mask.data, f_good)

    if wsp.perfusion_wm is not None:
        wsp.log.write(" - Found partial volume corrected results - will mask ROIs using 'base' GM/WM masks (PVE thresholds: %.2f / %.2f)\n" % (wsp.min_gm_thresh, wsp.min_wm_thresh))
        arrival_data = [
            {
                "suffix" : "_arrival_gm",
                "f" : wsp.arrival,
                "var" : wsp.arrival_var,
                "mask" : np.logical_and(effective_mask, wsp.structural.gm_pv_asl.data > wsp.min_gm_thresh),
            },
            {
                "suffix" : "_arrival_wm",
                "f" : wsp.arrival_wm,
                "var" : wsp.arrival_wm_var,
                "mask" : np.logical_and(effective_mask, wsp.structural.wm_pv_asl.data > wsp.min_wm_thresh),
            },
        ]
    else:
        wsp.log.write(" - No partial volume corrected results - will mask ROIs using 'pure' GM/WM masks (PVE thresholds: %.2f / %.2f)\n" % (wsp.pure_gm_thresh, wsp.pure_wm_thresh))
        arrival_data = [
            {
                "suffix" : "_arrival", 
                "f" : wsp.arrival,
                "var" :  wsp.arrival_var,
                "mask" : effective_mask,
            },
            {
                "suffix" : "_arrival_gm",
                "f" : wsp.arrival,
                "var" :  wsp.arrival_var,
                "mask" : np.logical_and(effective_mask, wsp.structural.gm_pv_asl.data > wsp.pure_gm_thresh),
            },
            {
                "suffix" : "_arrival_wm",
                "f" : wsp.perfusion,
                "var" :  wsp.arrival_var,
                "mask" : np.logical_and(effective_mask, wsp.structural.wm_pv_asl.data > wsp.pure_wm_thresh),
            },
        ]
    return arrival_data

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
        wsp.structural.wm_pv_asl = reg.change_space(wsp, wsp.structural.wm_pv, "asl")

    if wsp.pvwm is not None:
        wsp.structural.gm_pv_asl = wsp.pvgm
    else:
        wsp.structural.gm_pv_asl = reg.change_space(wsp, wsp.structural.gm_pv, "asl")
    
    if wsp.pvcsf is not None:
        wsp.structural.csf_pv_asl = wsp.pvcsf
    else:
        wsp.structural.csf_pv_asl = reg.change_space(wsp, wsp.structural.csf_pv, "asl")

    wsp.pure_gm_thresh, wsp.pure_wm_thresh = wsp.rois.pure_gm_thresh, wsp.rois.pure_wm_thresh
    wsp.min_gm_thresh, wsp.min_wm_thresh = wsp.ifnone("min_gm_thresh", 0.1), wsp.ifnone("min_wm_thresh", 0.1)

    rois = []
    wsp.log.write("\nLoading generic ROIs\n")
    add_roi(wsp, rois, "%i%%+GM" % (wsp.min_gm_thresh*100), wsp.structural.gm_pv_asl, wsp.min_gm_thresh)
    add_roi(wsp, rois, "%i%%+WM" % (wsp.min_wm_thresh*100), wsp.structural.wm_pv_asl, wsp.min_wm_thresh)
    add_roi(wsp, rois, "%i%%+GM" % (wsp.pure_gm_thresh*100), wsp.structural.gm_pv_asl, wsp.pure_gm_thresh)
    add_roi(wsp, rois, "%i%%+WM" % (wsp.pure_wm_thresh*100), wsp.structural.wm_pv_asl, wsp.pure_wm_thresh)
    add_roi(wsp, rois, "Cortical %i%%+GM" % (wsp.pure_gm_thresh*100), wsp.rois.cortical_gm_asl)
    add_roi(wsp, rois, "Cerebral %i%%+WM" % (wsp.pure_wm_thresh*100), wsp.rois.cerebral_wm_asl)

    wsp.log.write("\nLoading tissue PV ROI set")
    roi_set = Image(np.stack([wsp.structural.gm_pv_asl.data, wsp.structural.wm_pv_asl.data, wsp.structural.csf_pv_asl.data], axis=-1), header=wsp.structural.csf_pv_asl.header)
    add_roi_set(wsp, rois, ["GM PV", "WM PV", "CSF PV"], roi_set)

    # Add ROIs from command line
    user_rois = [l.strip() for l in wsp.ifnone("add_roi", "").split(",") if l.strip() != ""]
    if user_rois:
        wsp.log.write("\nLoading user-specified ROIs\n")
        for fname in user_rois:
            add_roi(wsp, rois, os.path.basename(fname).split(".")[0], Image(fname), 0.5)

    add_atlas = [l.strip() for l in wsp.ifnone("add_atlas", "").split(",") if l.strip() != ""]
    atlas_labels = [l.strip() for l in wsp.ifnone("add_atlas_labels", "").split(",") if l.strip() != ""]
    for idx, fname in enumerate(add_atlas):
        if idx < len(atlas_labels):
            with open(atlas_labels[idx]) as f:
                names = [l.strip() for l in f.readlines()]
        else:
            names = [os.path.basename(fname).split(".")[0],]
        atlas = Image(fname)
        if atlas.data.ndim == 3:
            add_rois_from_3d_label_atlas(wsp, rois, atlas, names)
        else:
            add_roi_set_from_4d_atlas(wsp, rois, atlas, names)

    # Add ROIs from standard atlases
    add_roi_set_from_fsl_atlas(wsp, rois, "harvardoxford-cortical", threshold=0.5)
    add_roi_set_from_fsl_atlas(wsp, rois, "harvardoxford-subcortical", threshold=0.5)

    for calib_method in wsp.calibration.calib_method:
        wsp.log.write("\nCalibration method: %s\n" % calib_method)
        if calib_method != "prequantified":
            calib_wsp = getattr(wsp.native, "calib_%s" % calib_method)
        else:
            calib_wsp = wsp.native

        wsp.log.write("\nLoading perfusion images\n")
        stats_data = get_perfusion_data(calib_wsp)
        if wsp.arrival_stats:
            wsp.log.write(" - Also generating stats for arrival data\n")
            stats_data += get_arrival_data(calib_wsp)

        # Get stats in each ROI. Add name to stats dict to make TSV output easier
        wsp.log.write("\nGetting stats - minimum of %i voxels to report in region\n" % wsp.roi_min_nvoxels)
        for data_item in stats_data:
            stats = []
            for roi in rois:
                roi_stats = {"name" : roi.get("name", None), "names" : roi.get("names", None)}
                get_stats(wsp, roi_stats, roi, data_item)

                # The ROI might just be a single ROI or it might be a 'fuzzy set' - need
                # to handle both cases
                roi_stats = dict(roi_stats)
                if roi_stats["name"] is not None:
                    roi_stats.pop("names")
                    stats.append(roi_stats)
                else:
                    for idx in range(len(roi_stats["names"])):
                        roi_stats_copy = dict(roi_stats)
                        for k in list(roi_stats.keys()):
                            if roi_stats_copy[k] is not None:
                                roi_stats_copy[k] = roi_stats_copy[k][idx]
                        roi_stats_copy["name"] = roi_stats_copy.pop("names")
                        stats.append(roi_stats_copy)

            columns = list(stats[0].keys())
            setattr(calib_wsp, "roi_stats%s" % data_item["suffix"], pd.DataFrame(stats, columns=columns))
    
    # Save output masks/PVE maps
    if wsp.save_asl_rois or wsp.save_asl_masks or wsp.save_struct_rois or wsp.save_std_rois:
        wsp.sub("region_rois")
        for roi in rois:
            if "name" in roi:
                fname = roi["name"].replace(" ", "_").replace(",", "").lower()
            else:
                fname = "_".join(roi["names"]).replace(" ", "_").replace(",", "").lower()
            if wsp.save_asl_rois and "roi_asl" in roi:
                setattr(wsp.region_rois, fname + "_asl_roi", roi["roi_asl"])
            if wsp.save_asl_masks and "mask_asl" in roi:
                setattr(wsp.region_rois, fname + "_asl_mask", roi["mask_asl"])
            if wsp.save_struct_rois and "roi_struct" in roi:
                setattr(wsp.region_rois, fname + "_struct_roi", roi["roi_struct"])
            if wsp.save_std_rois and "roi_std" in roi:
                setattr(wsp.region_rois, fname + "_std_roi", roi["roi_std"])

    wsp.log.write("\nDONE\n")
