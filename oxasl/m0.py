#!/bin/env python
"""
OXASL - MOdule to calculate M0 value or image for calibration

Copyright (c) 2008 University of Nottingham
"""

import os
import math
from optparse import OptionGroup

import numpy as np
import scipy.ndimage

from fsl.data.image import Image
from fsl.data.atlases import AtlasRegistry

from oxasl import reg
from oxasl.reporting import LightboxImage

def add_options(parser):
    group = OptionGroup(parser, "Calibration")
    group.add_option("--calib", "-c", help="Calibration image", type="image")
    group.add_option("--calib-method", "--cmethod", help="Calibration method: voxelwise / refregion. Multiple values may be specified (comma separated)")
    group.add_option("--calib-m0", help="Specify a precalculated M0 value for reference region calibration", type=float)
    group.add_option("--calib-m0img", help="Specify a precalculated voxelwise M0 image for calibration", type="image")
    group.add_option("--calib-alpha", "--alpha", help="Inversion efficiency", type=float, default=None)
    group.add_option("--calib-gain", "--cgain", help="Relative gain between calibration and ASL data", type=float, default=1.0)
    group.add_option("--calib-aslreg", help="Calibration image is already aligned with ASL image", action="store_true", default=False)
    group.add_option("--tr", help="TR used in calibration sequence (s)", type=float, default=3.2)
    parser.add_option_group(group)

    group = OptionGroup(parser, "Voxelwise calibration")
    group.add_option("--pct", help="Tissue/arterial partition coefficiant", type=float, default=0.9)
    group.add_option("--calib-no-edgecorr", dest="calib_edgecorr", help="Disable edge correction on M0 image", action="store_false", default=True)
    parser.add_option_group(group)

    # FIXME multiple values for refmask, t1r, t2r, pcr
    group = OptionGroup(parser, "Reference region calibration")
    group.add_option("--mode", help="Calibration mode (longtr or satrevoc)", default="longtr")
    group.add_option("--tissref", help="Tissue reference type (csf, wm, gm). Multiple values may be specified (comma separated)", default="csf")
    group.add_option("--te", help="Sequence TE (ms)", type=float, default=0.0)
    group.add_option("--refmask", help="Reference tissue mask in calibration image space for when single tissref is used", type="image")
    group.add_option("--csf", help="Reference CSF tissue mask in calibration image space", type="image")
    group.add_option("--wm", help="Reference WM tissue mask in calibration image space", type="image")
    group.add_option("--gm", help="Reference GM tissue mask in calibration image space", type="image")
    group.add_option("--t1r", help="T1 of reference tissue (s) - defaults: csf 4.3, gm 1.3, wm 1.0.", type=float, default=None)
    group.add_option("--t2r", help="T2/T2* of reference tissue (ms) - defaults T2/T2*: csf 750/400, gm 100/60,  wm 50/50", type=float, default=None)
    group.add_option("--t2star", action="store_true", default=False, help="Correct with T2* rather than T2 (alters the default T2 values)")
    group.add_option("--pcr", help="Reference tissue partition coefficiant (defaults csf 1.15, gm 0.98,  wm 0.82)", type=float, default=None)
    parser.add_option_group(group)

    group = OptionGroup(parser, "longtr mode (calibration image is a control image with a long TR)")
    parser.add_option_group(group)

    group = OptionGroup(parser, "satrecov mode (calibration image is a sequnce of control images at various TIs)")
    #group.add_option("--tis", help="Comma separated list of inversion times, e.group. --tis 0.2,0.4,0.6")
    group.add_option("--fa", help="Flip angle (in degrees) for Look-Locker readouts", type=float)
    group.add_option("--lfa", help="Lower flip angle (in degrees) for dual FA calibration", type=float)
    group.add_option("--calib-nphases", help="Number of phases (repetitions) of higher FA", type=int)
    group.add_option("--fixa", action="store_true", default=False, help="Fix the saturation efficiency to 100% (useful if you have a low number of samples)")
    parser.add_option_group(group)

    #group = OptionGroup(parser, "CSF masking options (only for --tissref csf)")
    #group.add_option("--csfmaskingoff", action="store_true", default=False, help="Turn off the ventricle masking, reference is based on segmentation only.")
    #group.add_option("--str2std", help="Structural to MNI152 linear registration (.mat)")
    #group.add_option("--warp", help="Structural to MNI152 non-linear registration (warp)")

def run(wsp):
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

     - ``calibration.<method>.m0`` : M0 value, either a single value or an voxelwise Image
    """
    wsp.sub("calibration")

    if wsp.asldata.iaf == "quant":
        wsp.log.write("\nCalibration\n")
        wsp.log.write(" - No calibration, data is already quantified\n")
        wsp.calibration.calib_method = ["prequantified"]
        wsp.calibration.sub("prequantified")
        wsp.calibration.prequantified.m0 = 1.0
        return

    wsp.calibration.calib_method = [m.strip().lower() for m in wsp.ifnone("calib_method", "voxelwise").split(",")]
    # Synonyms
    wsp.calibration.calib_method = ["voxelwise" if w == "voxel" else w for w in wsp.calibration.calib_method]
    wsp.calibration.calib_method = ["refregion" if w == "single" else w for w in wsp.calibration.calib_method]
    if any([c not in ("voxelwise", "refregion", "wholebrain") for c in wsp.calibration.calib_method]):
        raise ValueError("Unknown calibration method found in: %s" % wsp.calib_method)

    wsp.calibration.tissref = [t.strip().lower() for t in wsp.ifnone("tissref", "").split(",") if t.strip().lower() != ""]

    wsp.log.write("\nCalibration - calculating M0\n")
    if wsp.calib_m0 is not None:
        wsp.log.write(" - Using precalculated M0 value: %f\n" % wsp.calib_m0)
        wsp.calibration.sub("user")
        wsp.calibration.user.m0 = wsp.calib_m0
        wsp.calibration.calib_method = ["user"]
    elif wsp.calib_m0img is not None:
        wsp.log.write(" - Using precalculated voxelwise M0 image: %s\n" % wsp.calib_m0img.name)
        wsp.calibration.sub("user")
        wsp.calibration.user.m0 = wsp.calib_m0img
        wsp.calibration.calib_method = ["user"]
    elif wsp.calib is None:
        wsp.log.write(" - No calibration image and no M0 specified - calibration will not be performed\n")
        wsp.calibration.calib_method = []

    if "voxelwise" in wsp.calibration.calib_method:
        wsp.calibration.sub("voxelwise")
        get_m0_voxelwise(wsp.calibration.voxelwise)

    if "refregion" in wsp.calibration.calib_method:
        # FIXME if multiple tissue references are given there might also be multiple
        # user-specified override parameters (t1, t2 etc)
        wsp.calibration.calib_method.remove("refregion")
        for tissref in wsp.calibration.tissref:
            sub_wsp = wsp.calibration.sub("refregion_%s" % tissref)
            sub_wsp.tissref = tissref
            get_m0_refregion(sub_wsp)
            wsp.calibration.calib_method.append("refregion_%s" % tissref)

    if "wholebrain" in wsp.calibration.calib_method:
        wsp.calibration.sub("wholebrain")
        get_m0_wholebrain(wsp.calibration.wholebrain)

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
     - ``t1`` : Tissue T1 (optional, if present used for short TR correction)
     - ``pct`` : Partition coefficient used to convert M0 tissue into M0 arterial (default 0.9)
     - ``mask`` : Brain mask in calibration image space
     - ``calib_edgecorr`` : If True, and mask provided, apply edge correction
    """
    wsp.log.write(" - Doing voxelwise calibration\n")
    gain, pct = wsp.ifnone("calib_gain", 1), wsp.ifnone("pct", 0.9)
    wsp.log.write(" - Calibration gain: %f\n" % gain)

    # Calculate M0 value
    wsp.calib_img = wsp.corrected.calib
    calib_data = wsp.calib_img.data
    m0 = calib_data.astype(np.float32) * gain

    shorttr = 1
    if wsp.tr is not None and wsp.tr < 5:
        if wsp.t1 is not None:
            wsp.log.write(" - Correcting the calibration (M0) image for short TR (TR=%f, using T1 of tissue %f)\n" % (wsp.tr, wsp.t1))
            shorttr = 1 / (1 - math.exp(-wsp.tr / wsp.t1))
            m0 *= shorttr
        else:
            wsp.log.write("WARNING: tr < 5 (%f) but tissue T1 not provided so cannot apply short TR correction\n" % wsp.tr)

	# Include partiition co-effcient in M0 image to convert from M0 tissue to M0 arterial
    wsp.log.write(" - Using partition coefficient: %f\n" % pct)
    m0 /= pct

    edgecorr = wsp.ifnone("calib_edgecorr", True)
    if wsp.rois is not None and wsp.rois.mask is not None:
        if edgecorr:
            wsp.log.write(" - Doing edge correction\n")
            m0 = _edge_correct(m0, wsp.rois.mask)
        wsp.log.write(" - mean in mask: %f\n" % np.mean(m0[wsp.rois.mask.data > 0]))
    else:
        wsp.log.write(" - Mean M0 (no mask): %f\n" % np.mean(m0))

    wsp.m0 = Image(m0, header=wsp.calib_img.header)

    # Reporting
    page = wsp.report.page("m0_voxelwise")
    page.heading("Voxelwise M0 calculation")
    page.text("Voxelwise calibration calculates an M0 value for each voxel from the calibration image")
    page.heading("Correction factors", level=1)
    table = []
    table.append(["Calibration gain", "%.3g" % gain])
    table.append(["Blood/tissue partition coefficient", "%.3g" % pct])
    if wsp.tr is not None:
        table.append(["Sequence TR", "%.3g" % wsp.tr])
    if shorttr != 1:
        table.append(["Tissue T1", "%.3g" % wsp.t1])
        table.append(["T1 Correction factor (for short TR)", "%.3g" % shorttr])
    table.append(["Overall correction factor", "%.3g" % (gain*shorttr/pct)])

    if wsp.rois is not None and wsp.rois.mask is not None:
        table.append(["Mean M0 (within mask)", "%.3g" % np.mean(m0[wsp.rois.mask.data > 0])])
        table.append(["Edge correction", "Enabled" if edgecorr else "Disabled"])
    else:
        table.append(["Mean M0", "%.3g" % np.mean(m0)])
    page.table(table)

    page.heading("M0 map", level=1)
    page.image("m0img", LightboxImage(wsp.m0))

def _edge_correct(m0, brain_mask):
    """
    Correct for (partial volume) edge effects
    """
    brain_mask = brain_mask.data

    # Median smoothing. Note that we do this with a generic filter to ignore values
    # beyond the boundary (by setting them to nan)
    m0 = scipy.ndimage.generic_filter(m0, np.nanmedian, size=3, mode='constant', cval=np.nan)

    # Erode mask using 3x3x3 structuring element and zero M0 values outside the mask (i.e. edge voxels)
    mask_ero = scipy.ndimage.morphology.binary_erosion(brain_mask, structure=np.ones([3, 3, 3]), border_value=1)
    m0[mask_ero == 0] = 0

    # Extrapolate remaining data
    # ASL_FILE works slicewise using a mean 5x5 filter on nonzero values, so we will do the same as far as possible
    #
    # Note that we run the filter continuously until the whole volume is full. Subsequent runs do not
    # affect existing nonzero data, but we want to fill everything so if the mask changes (e.g. prior to PVC)
    # we don't end up with zero calibration data in an unmasked voxel
    #
    # Finally, there is an edge case where a z-slice is entirely zero after erosion. In this case we have nothing to
    # extrapolate from, so we avoid the 2D extrapolation step. A final 3D extrapolation picks up these slices.
    for z in range(m0.shape[2]):
        zslice_extrap = np.copy(m0[..., z])
        if np.all(np.isclose(zslice_extrap, 0)):
            continue
        while np.any(np.isclose(zslice_extrap, 0)):
            zslice_extrap = scipy.ndimage.filters.generic_filter(zslice_extrap, _masked_mean, footprint=np.ones([5, 5]))
        m0[..., z] = zslice_extrap

    # If any zeros remain because of completely empty slices, fill them in using 3D extrapolation
    while np.any(np.isclose(m0, 0)):
        m0 = scipy.ndimage.filters.generic_filter(m0, _masked_mean, footprint=np.ones([5, 5, 5]))

    return m0

    # Extrapolate remaining data
    # ASL_FILE works slicewise using a mean 5x5 filter on nonzero values, so we will do the same
    # Note that we run the filter continuously until the whole volume is full. Subsequent runs do not
    # affect existing nonzero data, but we want to fill everything so if the mask changes (e.g. prior to PVC)
    # we don't end up with zero calibration data in an unmasked voxel
    for z in range(m0.shape[2]):
        zslice_extrap = np.copy(m0[..., z])
        while np.any(np.isclose(zslice_extrap, 0)):
            zslice_extrap = scipy.ndimage.filters.generic_filter(zslice_extrap, _masked_mean, footprint=np.ones([5, 5]))
        m0[..., z] = zslice_extrap

    return m0

def _masked_mean(vals):
    """
    Called by scipy.ndimage.filters.generic_filter

    :param vals: Values in the kernel (a 5x5 square patch centered on the voxel in this case)

    For nonzero voxels, returns the voxel value (i.e. the middle element of vals).
    For zero voxels, returns the mean of non zero voxels in the kernel
    """
    voxel_val = vals[int((len(vals)-1) / 2)]
    if np.isclose(voxel_val, 0):
        nonzero = vals[~np.isclose(vals, 0)]
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

      - ``calib``     - Calibration Image in ASL space
      - ``rois.mask`` - Brain mask Image in ASL space
      - ``struc``     - Structural image
    """
    wsp.log.write("\n - Doing wholebrain region calibration\n")

    tr = wsp.ifnone("tr", 3.2)
    te = wsp.ifnone("te", 0)
    taq = wsp.ifnone("taq", 0)
    wsp.log.write(" - Using TE=%f, TR=%f, Readout time (TAQ)=%f\n" % (te, tr, taq))

    t2star = wsp.ifnone("t2star", False)
    if t2star:
        t2b = wsp.ifnone("t2sb", 50)
    else:
        t2b = wsp.ifnone("t2b", 150)

    # Check the data and masks
    wsp.calib_img = wsp.corrected.calib
    calib_data = wsp.calib_img.data

    ### Sensitivity image calculation (if we have a sensitivity image)
    if wsp.sens:
        wsp.log.write(" - Using sensitivity image: %s\n" % wsp.sens.name)
        calib_data /= wsp.sens.data

    m0 = np.zeros(calib_data.shape, dtype=np.float32)
    for tiss_type in ("wm", "gm", "csf"):
        pve_struc = getattr(wsp.structural, "%s_pv" % tiss_type)
        wsp.log.write(" - Transforming %s tissue PVE into ASL space\n" % tiss_type)
        pve = reg.change_space(wsp, pve_struc, "asl")
        t1r, t2r, t2sr, pcr = tissue_defaults(tiss_type)
        if t2star:
            t2r = t2sr
        t1_corr = 1 / (1 - math.exp(- (tr - taq) / t1r))
        t2_corr = math.exp(- te / t2b) / math.exp(- te / t2r)
        wsp.log.write("Correction factors: T1: %f, T2 %f, PC: %f" % (t1_corr, t2_corr, pcr))
        tiss_m0 = calib_data * t1_corr * t2_corr / pcr
        wsp.log.write(" - Mean %s M0: %f (weighted by PV)\n" % (tiss_type, np.average(tiss_m0, weights=pve.data)))
        tiss_m0 *= pve.data
        setattr(wsp, "m0_img_%s" % tiss_type, Image(tiss_m0, header=wsp.calib_img.header))
        m0 += tiss_m0

    gain = wsp.ifnone("calib_gain", 1)
    wsp.log.write(" - Calibration gain: %f\n" % gain)
    if gain != 1:
        m0 = m0 * gain

    wsp.m0_img = Image(m0, header=wsp.calib_img.header)
    wsp.m0 = np.mean(m0[wsp.rois.mask.data > 0])
    wsp.log.write(" - M0 of brain: %f\n" % wsp.m0)

    # Reporting
    page = wsp.report.page("m0_wholebrain")
    page.heading("Whole-brain M0 calculation")
    page.text("Whole-brain calibration calculates an M0 value for each voxel from the calibration image based on partial volume estimates")
    page.text("- Calibration gain: %f" % gain)
    page.text(" - TR: %f" % tr)
    page.text(" - TE: %f" % te)

    page.heading("M0", level=1)
    page.text("Mean M0 value (within mask): %f" % wsp.m0)
    page.image("m0img", LightboxImage(wsp.m0_img))

def get_m0_refregion(wsp, mode="longtr"):
    """
    Do reference region calibration

    FIXME saturation recovery mode is not functional

    Required workspace attributes
    -----------------------------

     - ``calib`` : Image containing voxelwise M0 map

    Optional workspace attributes
    -----------------------------

     - ``calib_gain`` : Calibration gain (default 1.0)
     - ``tr`` : Sequence TR (s) (default 3.2)
     - ``te`` : Sequence TE (s) (default 0)
     - ``taq`` : Sequence TAQ (s) (default 0)
     - ``t2star`` : If True, correct for T2* rather than T2 (i.e. use T2* defaults not T2 defaults)
     - ``t1r`` : Reference tissue T1 (default: see ``TISSUE_DEFAULTS``)
     - ``t2r`` : Reference tissue T2/T2* (default: see ``TISSUE_DEFAULTS``)
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
        t2b = wsp.ifnone("t2sb", 50)
    else:
        t2b = wsp.ifnone("t2b", 150)

    # Parameters for reference tissue type: T1, T2, T2*, partition coeffs, FAST seg ID
    # Partition coeffs based on Herscovitch and Raichle 1985 with a blood water density of 0.87
    # GM/WM T2* from Foucher 2011 JMRI 34:785-790
    wsp.tissref = wsp.ifnone("tissref", "csf")
    wsp.log.write(" - Using tissue reference type: %s\n" % wsp.tissref)
    if wsp.refmask is None:
        wsp.refmask = getattr(wsp, wsp.tissref)

    t1r, t2r, t2sr, pcr = tissue_defaults(wsp.tissref)
    if t2star:
        t2r = t2sr

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
    wsp.calib_img = wsp.corrected.calib
    calib_data = wsp.calib_img.data

    if wsp.rois is not None and wsp.rois.mask is not None:
        brain_mask = wsp.rois.mask.data
    else:
        brain_mask = np.ones(wsp.calib_img.shape[:3])

    if wsp.refmask is not None:
        wsp.log.write(" - Using supplied reference tissue mask: %s" % wsp.refmask.name)
        wsp.refmask = Image(wsp.refmask.data.astype(np.int32), header=wsp.refmask.header)
        if wsp.calib_aslreg:
            wsp.log.write(" (Aligned to ASL image already)\n")
            wsp.refmask_trans = wsp.refmask
        else:
            wsp.log.write(" (Transforming to ASL image space)\n")
            wsp.refmask_trans = reg.change_space(wsp, wsp.refmask, "asl", source_space="calib", mask=True)
        refmask = wsp.refmask_trans.data
    elif wsp.tissref.lower() in ("csf", "wm", "gm"):
        get_tissrefmask(wsp)
        refmask = wsp.refmask.data

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
            "data" : calib_data,
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
        wsp.t1_ref = Image("%s/mean_T1t" % wsp.workdir)
        wsp.m0_ref = Image("%s/mean_M0t" % wsp.workdir)

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
    wsp.m0 = mean_sig * gain * t1_corr * t2_corr / pcr

    # Apply calibration to input image
    if sens_corr:
        # Apply sensitivity image
        # fslmaths $infile -div $temp_calib/sens $temp_calib/infile
        pass
    else:
        # imcp $infile $temp_calib/infile
        pass

    wsp.log.write(" - M0: %f\n" % wsp.m0)

    # Reporting
    page = wsp.report.page("m0_refregion_%s" % wsp.tissref)
    page.heading("Reference region M0 calculation for %s" % wsp.tissref)
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
        ["M0", wsp.m0],
    ])
    page.heading("Reference tissue mask", level=1)
    page.image("refmask", LightboxImage(wsp.refmask, bgimage=wsp.calib_img))

def get_tissrefmask(wsp):
    """
    Calculate a calibration reference mask for a particular known tissue type
    """
    page = wsp.report.page("auto_calib_mask_%s" % wsp.tissref)
    page.heading("Calibration reference region: %s" % wsp.tissref)
    page.text("Reference region was automatically generated for tissue type: %s" % wsp.tissref.upper())
    page.heading("Partial volume map for %s tissue (from structural segmentation)" % wsp.tissref.upper(), level=1)
    refpve = getattr(wsp.structural, "%s_pv" % wsp.tissref.lower())
    wsp.refpve = reg.change_space(wsp, refpve, "struc")
    page.image("refpve", LightboxImage(wsp.refpve, bgimage=wsp.structural.brain))

    if wsp.tissref == "csf" and not wsp.csfmaskingoff:
        wsp.log.write(" - Doing automatic ventricle selection using standard atlas\n")
        # By deafult now we do FNRIT transformation of ventricle mask
        # FIXME disabled as not being used in ASL_CALIB at present
        #wsp.struc2stdfnirt = wsp.ifnone("stdmaskfnirt", True)

        # Select ventricles based on standard space atlas
        page.heading("Automatic ventricle selection", level=1)
        page.text("Standard space ventricles mask (from Harvard-Oxford atlas) eroded by 1 pixel")
        atlases = AtlasRegistry()
        atlases.rescanAtlases()
        atlas = atlases.loadAtlas("harvardoxford-subcortical", loadSummary=False, resolution=2)
        ventricles = ((atlas.data[..., 2] + atlas.data[..., 13]) > 0.1).astype(np.int32)
        wsp.ventricles = Image(scipy.ndimage.binary_erosion(ventricles, structure=np.ones([3, 3, 3]), border_value=1).astype(np.int32), header=atlas.header)
        std_img = Image(os.path.join(os.environ["FSLDIR"], "data", "standard", 'MNI152_T1_2mm_brain'))
        page.image("ventricles_std", LightboxImage(wsp.ventricles, bgimage=std_img))

        page.heading("Structural space ventricles mask", level=1)
        # FIXME nearest neighbour interpolation?
        wsp.ventricles_struc = reg.change_space(wsp, wsp.ventricles, "struc")
        page.text("This is the above image transformed into structural space. The transformation was obtained by registering the structural image to the standard brain image")
        page.image("ventricles_struc", LightboxImage(wsp.ventricles_struc, bgimage=wsp.structural.brain))

        wsp.log.write(" - Masking FAST output with standard space derived ventricle mask\n")
        wsp.refpve_pre_mask = wsp.refpve
        refpve_data = np.copy(wsp.refpve.data)
        refpve_data[wsp.ventricles_struc.data == 0] = 0
        wsp.refpve = Image(refpve_data, header=wsp.refpve.header)
        wsp.refpve_post = wsp.refpve

        page.heading("Structural space ventricles PVE", level=1)
        page.text("This is the CSF partial volume masked by the ventricles mask. It should select only the ventricles from the original partial volume image.")
        page.image("refpve_post", LightboxImage(wsp.refpve, bgimage=wsp.structural.brain))

    wsp.log.write(" - Transforming tissue reference mask into ASL space\n")
    # FIXME calibration image may not be in ASL space! Oxford_asl does not handle this currently
    wsp.refpve_calib = reg.change_space(wsp, wsp.refpve, "asl")
    #wsp.refpve_calib.data[wsp.refpve_calib.data < 0.001] = 0 # Better for display
    page.heading("Reference region in ASL space", level=1)
    page.text("Partial volume map")
    page.image("refpve_calib", LightboxImage(wsp.refpve_calib, bgimage=wsp.calib_img))

    # Threshold reference mask conservatively to select only reference tissue
    wsp.log.write(" - Thresholding reference mask\n")
    wsp.refmask = Image((wsp.refpve_calib.data > 0.9).astype(np.int32), header=wsp.refpve_calib.header)

    page.text("Reference Mask (thresholded at 0.9")
    page.image("refmask", LightboxImage(wsp.refmask, bgimage=wsp.calib_img))

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
