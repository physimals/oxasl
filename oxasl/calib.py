#!/bin/env python
"""
ASL_CALIB: Calibration for ASL data

Michael Chappell & Brad MacIntosh, FMRIB Image Analysis & Physics Groups

Copyright (c) 2008-2013 Univerisity of Oxford
"""

import sys
import os
import math
import traceback

import numpy as np
import scipy.ndimage

from fsl.data.image import Image
from fsl.data.atlases import AtlasRegistry

from oxasl import Workspace, struc, reg
from oxasl.image import summary
from oxasl.options import AslOptionParser, OptionCategory, IgnorableOptionGroup, GenericOptions
from oxasl.reporting import LightboxImage

def init(wsp):
    """ Initialize calibration sub-workspace """
    if wsp.calibration is None:
        wsp.sub("calibration")

def calculate_m0(wsp):
    """
    Calculate M0 value for use in calibration of perfusion images

    :param wsp: Workspace object

    Required workspace attributes
    -----------------------------

     - ``calib`` : Image containing voxelwise M0 map
     - ``calib_method`` : ``voxelwise`` or ``refregion``. Voxelwise calibration calibrates each voxel separately using
                          the corresponding voxel in the M0 map. Reference region calibration determines a single M0
                          value for the entire image by averaging over a region of the calibration image corresponding
                          to a known tissue type (e.g. CSF).

    Optional workspace attributes
    -----------------------------

     - ``calib_gain`` : Calibration gain (default 1.0)

    Additional optional and mandatory attributes may be required for different methods - see :ref get_m0_voxelwise:
    and :ref get_m0_refregion: functions for details.

    Workspace attributes updated
    -----------------------------

     - ``calibration.m0`` : M0 value, either a single value or an voxelwise Image
    """
    init(wsp)

    if wsp.calibration.m0 is None:
        wsp.log.write("\nCalibration - calculating M0\n")
        if wsp.calib_method in ("voxel", "voxelwise"):
            wsp.calibration.m0 = get_m0_voxelwise(wsp)
        elif wsp.calib_method in ("refregion", "single"):
            wsp.calibration.m0 = get_m0_refregion(wsp)
        elif wsp.calib_method == "wholebrain":
            wsp.calibration.m0 = get_m0_wholebrain(wsp)
        else:
            raise ValueError("Unknown calibration method: %s" % wsp.calib_method)

def calibrate(wsp, perf_img, multiplier=1.0, alpha=1.0, var=False):
    """
    Do calibration of a perfusion image from a calibration (M0) image

    :param wsp: Workspace object
    :param perf_img: Image containing perfusion data to calibrate
    :param multiplier: Scalar multiple to convert output to physical units
    :param alpha: Inversion efficiency
    :param var: If True, assume data represents variance rather than value

    :return: Image containing calibrated data

    Required workspace attributes
    -----------------------------

     - ``calibration.m0`` : M0 single value or voxelwise Image
    """
    if not perf_img:
        raise ValueError("Perfusion data cannot be None")
    if not wsp.calib:
        raise ValueError("No calibration data supplied")

    init(wsp)
    calculate_m0(wsp)
    wsp.log.write("\nCalibrating perfusion data: %s\n" % perf_img.name)
    m0 = wsp.calibration.m0
    if isinstance(m0, Image):
        m0 = m0.data

    if var:
        wsp.log.write(" - Treating data as variance - squaring M0 correction, multiplier and inversion efficiency\n")
        m0 = m0**2
        multiplier = multiplier**2
        alpha = alpha**2

    if isinstance(m0, np.ndarray):
        # If M0 is zero, make calibrated data zero
        calibrated = np.zeros(perf_img.shape)
        calibrated[m0 > 0] = perf_img.data[m0 > 0] / m0[m0 > 0]
    else:
        calibrated = perf_img.data / m0

    if alpha != 1.0:
        wsp.log.write(" - Using inversion efficiency correction: %f\n" % alpha)
        calibrated /= alpha

    if multiplier != 1.0:
        wsp.log.write(" - Using multiplier for physical units: %f\n" % multiplier)
        calibrated *= multiplier

    perf_calib = Image(calibrated, name=perf_img.name + "_calib", header=perf_img.header)
    return perf_calib

def get_m0_voxelwise(wsp):
    """
    Calculate M0 value using voxelwise calibration

    :param wsp: Workspace object
    :return: Image containing voxelwise M0 map

    Required workspace attributes
    -----------------------------

     - ``calib`` : Image containing calibration data

    Optional workspace attributes
    -----------------------------

     - ``calib_gain`` : Calibration gain (default 1.0)
     - ``tr`` : Sequence TR (s) (optional, if present used for short TR correction when tr<5)
     - ``t1t`` : Tissue T1 (optional, if present used for short TR correction)
     - ``pct`` : Partition coefficient used to convert M0 tissue into M0 arterial (default 0.9)
     - ``mask`` : Brain mask in calibration image space
     - ``calib_edgecorr`` : If True, and mask provided, apply edge correction
    """
    wsp.log.write(" - Doing voxelwise calibration\n")
    gain, pct = wsp.ifnone("calib_gain", 1), wsp.ifnone("pct", 0.9)
    wsp.log.write(" - Calibration gain: %f\n" % gain)

    # Calculate M0 value
    m0 = np.copy(wsp.calib.data).astype(np.float) * gain

    shorttr = 1
    if wsp.tr is not None and wsp.tr < 5:
        if wsp.t1t is not None:
            wsp.log.write(" - Correcting the calibration (M0) image for short TR (TR=%f, using T1 of tissue %f)\n" % (wsp.tr, wsp.t1t))
            shorttr = 1 / (1 - math.exp(-wsp.tr / wsp.t1t))
            m0 *= shorttr
        else:
            wsp.log.write("WARNING: tr < 5 (%f) but tissue T1 not provided so cannot apply short TR correction\n" % wsp.tr)

	# Include partiition co-effcient in M0 image to convert from M0 tissue to M0 arterial
    wsp.log.write(" - Using partition coefficient: %f\n" % pct)
    m0 /= pct

    if wsp.rois is not None and wsp.rois.mask is not None:
        if wsp.ifnone("calib_edgecorr", True):
            wsp.log.write(" - Doing edge correction\n")
            m0 = _edge_correct(m0, wsp.rois.mask)
        wsp.log.write(" - Masking M0 image")
        m0[wsp.rois.mask.data == 0] = 0

    m0img = Image(m0, header=wsp.calib.header)
    wsp.log.write(" - Mean M0: %f\n" % np.mean(m0))

    # Reporting
    page = wsp.report.page("m0")
    page.heading("Voxelwise M0 calculation")
    page.text("Voxelwise calibration calculates an M0 value for each voxel from the calibration image")
    page.heading("Correction factors", level=1)
    table = []
    table.append(["Calibration gain", "%.3g" % gain])
    table.append(["Blood/tissue partition coefficient", "%.3g" % pct])
    if wsp.tr is not None:
        table.append(["Sequence TR", "%.3g" % wsp.tr])
    if shorttr != 1:
        table.append(["Tissue T1", "%.3g" % wsp.t1t])
        table.append(["T1 Correction factor (for short TR)", "%.3g" % shorttr])
    table.append(["Overall correction factor", "%.3g" % (gain*shorttr/pct)])

    if wsp.rois is not None and wsp.rois.mask is not None:
        table.append(["Mean M0 (within mask)", "%.3g" % np.mean(m0[wsp.rois.mask.data > 0])])
        table.append(["Edge correction", "Enabled" if wsp.edgecorr else "Disabled"])
    else:
        table.append(["Mean M0", "%.3g" % np.mean(m0)])
    page.table(table)

    page.heading("M0 map", level=1)
    page.image("m0img", LightboxImage(m0img))

    return m0img

def _edge_correct(m0, brain_mask):
    """
    Correct for (partial volume) edge effects
    """
    brain_mask = brain_mask.data

    # Median smoothing
    m0 = scipy.ndimage.median_filter(m0, size=3)

    # Erode mask using 3x3x3 structuring element
    mask_ero = scipy.ndimage.morphology.binary_erosion(brain_mask, structure=np.ones([3, 3, 3]), border_value=1)
    m0[mask_ero == 0] = 0

    # Extrapolate remaining data to fit original mask
    # ASL_FILE works slicewise using a mean 5x5 filter on nonzero values, so we will do the same
    for z in range(m0.shape[2]):
        zslice = m0[..., z]
        zslice_extrap = scipy.ndimage.filters.generic_filter(zslice, _masked_mean, footprint=np.ones([5, 5]))
        m0[..., z] = zslice_extrap
    m0[brain_mask == 0] = 0

    return m0

def _masked_mean(vals):
    """
    Called by scipy.ndimage.filters.generic_filter

    :param vals: Values in the kernel (a 5x5 square patch centered on the voxel in this case)

    For nonzero voxels, returns the voxel value (i.e. the middle element of vals).
    For zero voxels, returns the mean of non zero voxels in the kernel
    """
    voxel_val = vals[int((len(vals)-1) / 2)]
    if voxel_val == 0:
        nonzero = vals[vals != 0]
        if np.any(nonzero):
            return np.mean(nonzero)
        else:
            return 0
    else:
        return voxel_val

def get_m0_wholebrain(wsp):
    """
    Get a whole-brain M0 value

    This calculates an M0 map for the whole brain using the T1, T2 and PC
    values for each tissue type with the contribution at each voxel weighted
    by the tissue partial volume

    :param wsp: Workspace object
    :return: Image containing voxelwise M0 map

    Required Workspace attributes
    -----------------------------

      - ``calib``     - Calibration Image in ASL native space
      - ``rois.mask`` - Brain mask Image in ASL native space
      - ``struc``     - Structural image
    """
    struc.segment(wsp)
    wsp.log.write("\n - Doing wholebrain region calibration\n")

    tr = wsp.ifnone("tr", 3.2)
    te = wsp.ifnone("te", 0)
    taq = wsp.ifnone("taq", 0)
    wsp.log.write(" - Using TE=%f, TR=%f, Readout time (TAQ)=%f\n" % (te, tr, taq))

    t2star = wsp.ifnone("t2star", False)
    if t2star:
        # From Petersen 2006 MRM 55(2):219-232 see discussion
        t2b = 50
    else:
        # From Lu et al. 2012 MRM 67:42-49 have 154ms at 3T during normoxia
        t2b = 150

    # Check the data and masks
    calib_data = np.copy(wsp.calib.data)
    if wsp.rois is not None and wsp.rois.mask is not None:
        brain_mask = wsp.rois.mask.data
    else:
        brain_mask = np.ones(wsp.calib.shape[:3])

    ### Sensitivity image calculation (if we have a sensitivity image)
    if wsp.sens:
        wsp.log.write(" - Using sensitivity image: %s\n" % wsp.sens.name)
        calib_data /= wsp.sens.data

    m0 = np.zeros(calib_data.shape, dtype=np.float)
    for tiss_type in ("wm", "gm", "csf"):
        pve_struc = getattr(wsp.structural, "%s_pv" % tiss_type)
        wsp.log.write(" - Transforming %s tissue PVE into ASL space\n" % tiss_type)
        pve = reg.struc2asl(wsp, pve_struc)
        t1r, t2r, t2rstar, pcr = tissue_defaults(tiss_type)
        if t2star:
            t2r = t2rstar
        t1_corr = 1 / (1 - math.exp(- (tr - taq) / t1r))
        t2_corr = 1 / math.exp(- te / t2r)
        wsp.log.write("Correction factors: T1: %f, T2 %f, PC: %f" % (t1_corr, t2_corr, pcr))
        tiss_m0 = calib_data * t1_corr * t2_corr / pcr
        wsp.log.write(" - Mean %s M0: %f (weighted by PV)\n" % (tiss_type, np.average(tiss_m0, weights=pve.data)))
        tiss_m0 *= pve.data
        setattr(wsp.calibration, "m0_img_%s" % tiss_type, Image(tiss_m0, header=wsp.calib.header))
        m0 += tiss_m0

    m0 = m0 * math.exp(- te / t2b)
    gain = wsp.ifnone("calib_gain", 1)
    wsp.log.write(" - Calibration gain: %f\n" % gain)
    if gain != 1:
        m0 = m0 * gain

    wsp.calibration.m0_img = Image(m0, header=wsp.calib.header)
    m0 = np.mean(m0[brain_mask != 0])
    wsp.log.write(" - M0 of brain: %f\n" % m0)

    # Reporting
    page = wsp.report.page("m0")
    page.heading("Whole-brain M0 calculation")
    page.text("Whole-brain calibration calculates an M0 value for each voxel from the calibration image based on partial volume estimates")
    page.text("- Calibration gain: %f" % gain)
    page.text(" - TR: %f" % tr)
    page.text(" - TE: %f" % te)

    page.heading("M0", level=1)
    page.text("Mean M0 value (within mask): %f" % m0)
    page.image("m0img", LightboxImage(wsp.calibration.m0_img))

    return float(m0)

def get_m0_refregion(wsp, mode="longtr"):
    """
    Do reference region calibration

    FIXME saturation recovery mode is not functional

    Required workspace attributes
    -----------------------------

     - ``calib`` : Image containing voxelwise M0 map
     - ``refmask`` : Reference region mask in calibration image space. Technically this is not required
                      as the brain mask will be used if not specified, however it is strongly recommended.

    Optional workspace attributes
    -----------------------------

     - ``calib_gain`` : Calibration gain (default 1.0)
     - ``tr`` : Sequence TR (s) (default 3.2)
     - ``te`` : Sequence TE (s) (default 0)
     - ``taq`` : Sequence TAQ (s) (default 0)
     - ``t2star`` : If True, correct for T2* rather than T2
     - ``t1r`` : Reference tissue T1 (default: see ``TISSUE_DEFAULTS``)
     - ``t2r`` : Reference tissue T2 (default: see ``TISSUE_DEFAULTS``)
     - ``t2rstar`` : Reference tissue T2* (default: see ``TISSUE_DEFAULTS``)
     - ``pcr`` : Reference tissue partition coefficient (default: see ``TISSUE_DEFAULTS``)
     - ``mask`` : Brain mask in calibration image space
    """
    wsp.log.write(" - Doing reference region calibration\n")
    gain = wsp.ifnone("calib_gain", 1)

    tr = wsp.ifnone("tr", 3.2)
    te = wsp.ifnone("te", 0)
    taq = wsp.ifnone("taq", 0)
    wsp.log.write(" - Acquisition: TE=%f, TR=%f, Readout time (TAQ)=%f\n" % (te, tr, taq))

    t2star = wsp.ifnone("t2star", False)
    if t2star:
        # From Petersen 2006 MRM 55(2):219-232 see discussion
        t2b = 50
    else:
        # From Lu et al. 2012 MRM 67:42-49 have 154ms at 3T during normoxia
        t2b = 150

    # Parameters for reference tissue type: T1, T2, T2*, partition coeffs, FAST seg ID
    # Partition coeffs based on Herscovitch and Raichle 1985 with a blood water density of 0.87
    # GM/WM T2* from Foucher 2011 JMRI 34:785-790
    wsp.tissref = wsp.ifnone("tissref", "csf")
    wsp.log.write(" - Using tissue reference type: %s\n" % wsp.tissref)

    t1r, t2r, t2rstar, pcr = tissue_defaults(wsp.tissref)
    if t2star:
        t2r = t2rstar

    # Command line override of default T1, T2, PC
    t1r_img, t2r_img = False, False

    if wsp.t1r is not None:
        t1r = wsp.t1r
        if isinstance(t1r, Image):
            wsp.log.write(" - Using T1 image for reference region: %s\n" % t1r.name)
            t1r_img = True
        else:
            wsp.log.write(" - Using user-specified T1r value: %f\n" % t1r)

    if wsp.t2r is not None:
        t2r = wsp.t2r
        if isinstance(t2r, Image):
            wsp.log.write(" - Using T2 image for reference region: %s\n" % t2r.name)
            t2r_img = True
        else:
            wsp.log.write(" - Using user-specified T2r value: %f\n" % t2r)

    if wsp.t2b is not None:
        t2b = wsp.t2b
        wsp.log.write(" - Using user-specified T2b value: %f\n" % t2b)

    if wsp.pcr is not None:
        pcr = wsp.pcr
        wsp.log.write(" - Using user-specified partition coefficient: %f\n" % pcr)

    if t1r is None:
        raise ValueError("T1 for reference tissue has not been set")
    if t2r is None:
        if te != 0:
            raise ValueError("T2 for reference tissue has not been set")
        else:
            t2r = 1.0
    if pcr is None:
        raise ValueError("Partition coefficient for reference tissue has not been set")

    wsp.log.write(" - T1r: %f; T2r: %f; T2b: %f; Part co-eff: %f\n" % (t1r, t2r, t2b, pcr))

    # Check the data and masks
    calib_data = np.copy(wsp.calib.data).astype(np.float)
    wsp.calibration.calib_img = wsp.calib
    if calib_data.ndim == 4:
        wsp.log.write(" - Taking mean across calibration images\n")
        calib_data = np.mean(calib_data, -1)

    if wsp.rois is not None and wsp.rois.mask is not None:
        brain_mask = wsp.rois.mask.data
    else:
        brain_mask = np.ones(wsp.calib.shape[:3])

    if wsp.refmask is not None:
        wsp.log.write(" - Using supplied reference tissue mask: %s\n" % wsp.refmask.name)
        wsp.calibration.refmask = Image(wsp.refmask.data.astype(np.int), header=wsp.refmask.header)
        wsp.calibration.refmask_trans = reg.calib2asl(wsp, wsp.calibration.refmask, mask=True)
        refmask = wsp.calibration.refmask_trans.data
    elif wsp.tissref.lower() in ("csf", "wm", "gm"):
        get_tissrefmask(wsp)
        refmask = wsp.calibration.refmask.data

    nonzero = np.count_nonzero(refmask)
    if nonzero < 1:
        raise ValueError("Reference mask does not contain any unmasked voxels")
    else:
        wsp.log.write(" - Number of voxels in tissue reference mask: %i\n" % nonzero)

    ### Sensitivity image calculation (if we have a sensitivity image)
    sens_corr = False
    if wsp.sens:
        wsp.log.write(" - Using sensitivity image: %s\n" % wsp.sens.name)
        sens_corr = True
        sens_data = wsp.sens.data

    wsp.log.write(" - MODE: %s\n" % mode)
    wsp.log.write(" - Calibration gain: %f\n" % gain)

    if mode == "longtr":
        if sens_corr:
            wsp.log.write(" - Applying sensitivity image\n")
            calib_data /= sens_data

        # Mask M0 map with tissue reference
        calib_data[refmask == 0] = 0

        # calcualte T1 of reference region (if a T1 image has been supplied)
        if t1r_img:
            t1r = np.mean(t1r.data[refmask != 0])
            wsp.log.write(" - Calculated T1 of reference tissue: %f\n" % t1r)

        # calcualte T2 of reference region (if a T2 image has been supplied)
        if t2r_img:
            t2r = np.mean(t2r.data[refmask != 0])
            wsp.log.write(" - Calculated T2 of reference tissue: %f\n" % t2r)

        # calculate M0_ref value
        mean_sig = np.mean(calib_data[refmask != 0])
        wsp.log.write(" - mean signal in reference tissue: %f\n" % mean_sig)
        t1_corr = 1 / (1 - math.exp(- (tr - taq) / t1r))
        wsp.log.write(" - T1 correction factor: %f\n" % t1_corr)

    elif mode == "satrecov":
        # Calibration image is control images and we want to do a saturation recovery fit
        # NB only do the fit in the CSF mask
        # FIXME this is not functional at the moment
        options = {
            "data" : wsp.calib,
            "mask" : refmask,
            "method" : "vb",
            "noise" : "white",
            "model" : "satrecov",
            "t1" : t1r,
        }

        #deal with TIs
        tis = wsp.ifnone("satrecov_tis", [])
        wsp.log.write(" - TIs: %s\n" % str(tis))
        for idx, ti in enumerate(tis):
            options["ti%i" % (idx+1)] = ti

        # Extra options for Look Locker
        if wsp.FA:
            options["FA"] = wsp.FA
        if wsp.calib_nphases:
            options["phases"] = wsp.calib_nphases
        if wsp.lfa:
            options["LFA"] = wsp.lfa

        # Extra sat recovery options
        if wsp.fixa:
            options["fixa"] = True

        # Do fabber within the tissue reference mask with a sensible T1 prior mean
        wsp.log.write(" - Running FABBER within reference tissue mask\n")
        from .wrappers import fabber
        fabber_result = fabber(options)
        mean_m0 = fabber_result["mean_M0t"]

        # Calculate M0 value - this is mean M0 of CSF at the TE of the sequence
        m0_value = np.mean(mean_m0.data[refmask != 0])

        wsp.log.write(" - M0 of reference tissue: %f\n" % m0_value)

        # Save useful results
        wsp.calibration.t1_ref = Image("%s/mean_T1t" % wsp.workdir)
        wsp.calibration.m0_ref = Image("%s/mean_M0t" % wsp.workdir)

        # Do fabber again within whole brain to get estimated T1 of tissue and FA correction (if LL)
        # (note that we do not apply sensitivity correction to the data here - thius is 'built-into' the M0t map)
        wsp.log.write(" - FABBER (again) within whole brain mask\n")
        options["mask"] = brain_mask
        fabber(options)
        # $fabber --data=$calib --mask=$bmask --output=$temp_calib/satrecovT --data-order=singlefile --model=satrecov --noise=white --method=vb $tislist $llopts $sropts

        # save useful results to specified output directory
        #imcp $temp_calib/satrecovT/mean_T1t $outdir/T1t
        #imcp $temp_calib/satrecovT/mean_M0t $outdir/M0t
        #if [ ! -z $lfa ]; then
        #imcp $temp_calib/satrecovT/mean_g $outdir/facorr
    else:
        raise ValueError("Unknown reference region mode: %s (Should be satrecov or longtr)" % mode)

    # Use equation to get the M0 value that is needed

    #  T2 correction
    t2_corr = math.exp(- te / t2b) / math.exp(- te / t2r)
    wsp.log.write(" - T2 correction factor: %f\n" % t2_corr)
    m0 = mean_sig * gain * t1_corr * t2_corr / pcr

    # Apply calibration to input image
    if sens_corr:
        # Apply sensitivity image
        # fslmaths $infile -div $temp_calib/sens $temp_calib/infile
        pass
    else:
        # imcp $infile $temp_calib/infile
        pass

    wsp.log.write(" - M0: %f\n" % m0)

    # Reporting
    page = wsp.report.page("m0")
    page.heading("Reference region M0 calculation")
    page.text("Reference region calibration calculates a single M0 value for a region of known tissue type")
    page.heading("Calculation details", level=1)
    page.table([
        ["Sequence TR (s)", tr],
        ["Sequence TE (ms)", te],
        ["T1 reference tissue (s)", t1r],
        ["T2 reference tissue (ms)", t2r],
        ["Blood T2 (ms)", t2b],
        ["Number of voxels in reference region", np.count_nonzero(refmask)],
        ["Mean signal in reference region", mean_sig],
        ["T1 correction factor", t1_corr],
        ["T2 correction factor", t2_corr],
        ["Partition coefficient", pcr],
        ["Calibration gain", gain],
        ["M0", m0],
    ])
    page.heading("Reference tissue mask", level=1)
    page.image("refmask", LightboxImage(wsp.calibration.refmask, bgimage=wsp.calib))

    return float(m0)

def get_tissrefmask(wsp):
    """
    Calculate a calibration reference mask for a particular known tissue type
    """
    struc.segment(wsp)

    page = wsp.report.page("auto_calib_mask")
    page.heading("Calibration reference region")
    page.text("Reference region was automatically generated for tissue type: %s" % wsp.tissref.upper())
    page.heading("Partial volume map for %s tissue (from structural segmentation)" % wsp.tissref.upper(), level=1)
    wsp.calibration.refpve = getattr(wsp.structural, "%s_pv" % wsp.tissref.lower())
    page.image("refpve", LightboxImage(wsp.calibration.refpve, bgimage=wsp.structural.brain))

    if wsp.tissref == "csf" and not wsp.csfmaskingoff:
        wsp.log.write(" - Doing automatic ventricle selection using standard atlas\n")
        # By deafult now we do FNRIT transformation of ventricle mask
        # FIXME disabled as not being used in ASL_CALIB at present
        #wsp.calibration.struc2stdfnirt = wsp.ifnone("stdmaskfnirt", True)

        # Select ventricles based on standard space atlas
        page.heading("Automatic ventricle selection", level=1)
        page.text("Standard space ventricles mask (from Harvard-Oxford atlas) eroded by 1 pixel")
        atlases = AtlasRegistry()
        atlases.rescanAtlases()
        atlas = atlases.loadAtlas("harvardoxford-subcortical", loadSummary=False, resolution=2)
        ventricles = ((atlas.data[..., 2] + atlas.data[..., 13]) > 0.1).astype(np.int)
        wsp.calibration.ventricles = Image(scipy.ndimage.binary_erosion(ventricles, structure=np.ones([3, 3, 3]), border_value=1).astype(np.int), header=atlas.header)
        std_img = Image(os.path.join(os.environ["FSLDIR"], "data", "standard", 'MNI152_T1_2mm_brain'))
        page.image("ventricles_std", LightboxImage(wsp.calibration.ventricles, bgimage=std_img))

        page.heading("Structural space ventricles mask", level=1)
        reg.reg_struc2std(wsp)
        # FIXME nearest neighbour interpolation?
        wsp.calibration.ventricles_struc = reg.std2struc(wsp, wsp.calibration.ventricles)
        page.text("This is the above image transformed into structural space. The transformation was obtained by registering the structural image to the standard brain image")
        page.image("ventricles_struc", LightboxImage(wsp.calibration.ventricles_struc, bgimage=wsp.structural.brain))

        wsp.log.write(" - Masking FAST output with standard space derived ventricle mask\n")
        wsp.calibration.refpve_pre_mask = wsp.calibration.refpve
        refpve_data = np.copy(wsp.calibration.refpve.data)
        refpve_data[wsp.calibration.ventricles_struc.data == 0] = 0
        wsp.calibration.refpve = Image(refpve_data, header=wsp.calibration.refpve.header)
        wsp.calibration.refpve_post = wsp.calibration.refpve

        page.heading("Structural space ventricles PVE", level=1)
        page.text("This is the CSF partial volume masked by the ventricles mask. It should select only the ventricles from the original partial volume image.")
        page.image("refpve_post", LightboxImage(wsp.calibration.refpve, bgimage=wsp.structural.brain))

    wsp.log.write(" - Transforming tissue reference mask into ASL space\n")
    # FIXME calibration image may not be in ASL space! Oxford_asl does not handle this currently
    wsp.calibration.refpve_calib = reg.struc2asl(wsp, wsp.calibration.refpve)
    #wsp.calibration.refpve_calib.data[wsp.calibration.refpve_calib.data < 0.001] = 0 # Better for display
    page.heading("Reference region in ASL space", level=1)
    page.text("Partial volume map")
    page.image("refpve_calib", LightboxImage(wsp.calibration.refpve_calib, bgimage=wsp.calib))

    # Threshold reference mask conservatively to select only reference tissue
    wsp.log.write(" - Thresholding reference mask\n")
    wsp.calibration.refmask = Image((wsp.calibration.refpve_calib.data > 0.9).astype(np.int), header=wsp.calibration.refpve_calib.header)

    page.text("Reference Mask (thresholded at 0.9")
    page.image("refmask", LightboxImage(wsp.calibration.refmask, bgimage=wsp.calib))

TISS_DEFAULTS = {
    "csf" : [4.3, 750, 400, 1.15],
    "gm" : [1.3, 100, 50, 0.98],
    "wm" : [1.0, 50, 60, 0.82],
}

def tissue_defaults(tiss_type=None):
    """
    Get default T1, T2, T2* and PC for different tissue types

    Parameters for reference tissue type: T1, T2, T2*, partition coeffs, FAST seg ID
    Partition coeffs based on Herscovitch and Raichle 1985 with a blood water density of 0.87
    GM/WM T2* from Foucher 2011 JMRI 34:785-790

    :return: If tiss_type given, tuple of T1, T2, T2* and PC for tissue type.
             Otherwise return dictionary of tissue type name (csf, wm, gm) to tuple
             of T1, T2, T2* and PC for all known types
    """
    if tiss_type is None:
        return TISS_DEFAULTS
    elif tiss_type.lower() in TISS_DEFAULTS:
        return TISS_DEFAULTS[tiss_type.lower()]
    else:
        raise ValueError("Invalid tissue type: %s" % tiss_type)

class CalibOptions(OptionCategory):
    """
    OptionCategory which contains options for preprocessing ASL data
    """

    def __init__(self, **kwargs):
        OptionCategory.__init__(self, "calib", **kwargs)

    def groups(self, parser):
        groups = []
        group = IgnorableOptionGroup(parser, "Calibration", ignore=self.ignore)
        group.add_option("--calib", "-c", help="Calibration image", type="image")
        group.add_option("--perf", "-i", help="Perfusion image for calibration, in same image space as calibration image", type="image")
        group.add_option("--calib-method", "--cmethod", help="Calibration method: voxelwise or refregion")
        group.add_option("--calib-alpha", "--alpha", help="Inversion efficiency", type=float, default=None)
        group.add_option("--calib-gain", "--cgain", help="Relative gain between calibration and ASL data", type=float, default=1.0)
        group.add_option("--tr", help="TR used in calibration sequence (s)", type=float, default=3.2)
        groups.append(group)

        group = IgnorableOptionGroup(parser, "Voxelwise calibration", ignore=self.ignore)
        group.add_option("--pct", help="Tissue/arterial partition coefficiant", type=float, default=0.9)
        group.add_option("--t1t", help="T1 of tissue (s)", type=float, default=1.3)
        groups.append(group)

        group = IgnorableOptionGroup(parser, "Reference region calibration", ignore=self.ignore)
        group.add_option("--mode", help="Calibration mode (longtr or satrevoc)", default="longtr")
        group.add_option("--tissref", help="Tissue reference type (csf, wm, gm or none)", default="csf")
        group.add_option("--te", help="Sequence TE (ms)", type=float, default=0.0)
        group.add_option("--t1r", help="T1 of reference tissue (s) - defaults: csf 4.3, gm 1.3, wm 1.0", type=float, default=None)
        group.add_option("--t2r", help="T2/T2* of reference tissue (ms) - defaults T2/T2*: csf 750/400, gm 100/60,  wm 50/50", type=float, default=None)
        group.add_option("--t2b", help="T2/T2* of blood (ms) - default T2/T2*: 150/50)", type=float, default=None)
        group.add_option("--refmask", "--csf", help="Reference tissue mask in calibration image space", type="image")
        group.add_option("--t2star", action="store_true", default=False, help="Correct with T2* rather than T2 (alters the default T2 values)")
        group.add_option("--pcr", help="Reference tissue partition coefficiant (defaults csf 1.15, gm 0.98,  wm 0.82)", type=float, default=None)
        groups.append(group)

        group = IgnorableOptionGroup(parser, "longtr mode (calibration image is a control image with a long TR)", ignore=self.ignore)
        groups.append(group)

        group = IgnorableOptionGroup(parser, "satrecov mode (calibration image is a sequnce of control images at various TIs)", ignore=self.ignore)
        group.add_option("--tis", help="Comma separated list of inversion times, e.group. --tis 0.2,0.4,0.6")
        group.add_option("--fa", help="Flip angle (in degrees) for Look-Locker readouts", type=float)
        group.add_option("--lfa", help="Lower flip angle (in degrees) for dual FA calibration", type=float)
        group.add_option("--calib-nphases", help="Number of phases (repetitions) of higher FA", type=int)
        group.add_option("--fixa", action="store_true", default=False, help="Fix the saturation efficiency to 100% (useful if you have a low number of samples)")
        groups.append(group)

        #group = IgnorableOptionGroup(parser, "Coil sensitivity correction, either using existing sensitivity image or reference images collected using same parameters", ignore=self.ignore)
        #groups.append(group)

        #group = IgnorableOptionGroup(parser, "CSF masking options (only for --tissref csf)", ignore=self.ignore)
        #group.add_option("--csfmaskingoff", action="store_true", default=False, help="Turn off the ventricle masking, reference is based on segmentation only.")
        #group.add_option("--str2std", help="Structural to MNI152 linear registration (.mat)")
        #group.add_option("--warp", help="Structural to MNI152 non-linear registration (warp)")

        return groups

def main():
    """
    Entry point for oxasl_calib command line program
    """

    debug = False
    try:
        parser = AslOptionParser(usage="oxasl_calib -i <perfusion image> -c <calibration image> --calib-method <voxelwise|refregion> -o <output filename> [options]")
        parser.add_category(CalibOptions())
        parser.add_category(GenericOptions(output_type="file"))
        options, _ = parser.parse_args(sys.argv)

        if not options.perf:
            sys.stderr.write("Perfusion input file not specified\n")
            parser.print_help()
            sys.exit(1)

        if not options.calib:
            sys.stderr.write("Calibration input file not specified\n")
            parser.print_help()
            sys.exit(1)

        wsp = Workspace(**vars(options))

        summary(wsp.perf)
        summary(wsp.calib)

        if wsp.output is None:
            wsp.output = "%s_calib" % wsp.perf.name

        calibrated_img = calibrate(wsp, wsp.perf)
        summary(calibrated_img)
        calibrated_img.save(wsp.output)

    except ValueError as exc:
        sys.stderr.write("ERROR: " + str(exc) + "\n")
        if debug:
            traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()
