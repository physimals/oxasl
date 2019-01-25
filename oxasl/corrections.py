#!/bin/env python
"""
Motion and distortion corrections for ASL

The functions in this module calculate linear or non-linear correction transformations
to apply to the ASL and calibration images. Once calculated the ``apply_corrections``
function generates corrected images with the minimum of interpolation.

Currently four sources of transformation exist:

 - Motion correction of the ASL data. This generates a series of linear (rigid body)
   transformations in ASL space, one for each ASL volume. If calibration data is also 
   present a calibration->ASL transform is also generated as part of this process

 - Fieldmap-based distortion correction. This generates a nonlinear warp image
   in structural space which is then transformed to ASL space

 - Phase encoding reversed (CBLIP) distortion correction using TOPUP. This generates
   a nonlinear warp image in ASL space FIXME calibration space?

 - User-supplied nonlinear warp image for gradient distortion corection

Except for the TOPUP correction, all of the above can be combined in a single
transformation to minimise interpolation of the ASL data

Copyright (c) 2008-2013 Univerisity of Oxford
"""
from __future__ import unicode_literals

import tempfile
import shutil

import numpy as np

import fsl.wrappers as fsl
from fsl.data.image import Image

from oxasl import reg, struc
from oxasl.options import OptionCategory, IgnorableOptionGroup
from oxasl.reporting import LightboxImage, LineGraph
from oxasl.wrappers import epi_reg, fnirtfileutils

class DistcorrOptions(OptionCategory):
    """
    OptionCategory which contains options for distortion correction
    """

    def __init__(self, **kwargs):
        OptionCategory.__init__(self, "distcorr", **kwargs)

    def groups(self, parser):
        ret = []
        g = IgnorableOptionGroup(parser, "Distortion correction using fieldmap")
        g.add_option("--fmap", help="fieldmap image (in rad/s)", type="image")
        g.add_option("--fmapmag", help="fieldmap magnitude image - wholehead extracted", type="image")
        g.add_option("--fmapmagbrain", help="fieldmap magnitude image - brain extracted", type="image")
        g.add_option("--nofmapreg", help="Do not perform registration of fmap to T1 (use if fmap already in T1-space)", action="store_true", default=False)
        ret.append(g)

        g = IgnorableOptionGroup(parser, "Distortion correction using phase-encode-reversed calibration image (TOPUP)")
        g.add_option("--cblip", help="phase-encode-reversed (blipped) calibration image", type="image")
        ret.append(g)

        g = IgnorableOptionGroup(parser, "General distortion correction options")
        g.add_option("--echospacing", help="Effective EPI echo spacing (sometimes called dwell time) - in seconds", type=float)
        g.add_option("--pedir", help="Phase encoding direction, dir = x/y/z/-x/-y/-z")
        g.add_option("--gdcwarp", help="Additional warp image for gradient distortion correction - will be combined with fieldmap or TOPUP distortion correction", type="image")
        ret.append(g)

        g = IgnorableOptionGroup(parser, "Sensitivity correction")
        g.add_option("--cref", help="Reference image for sensitivity correction", type="image")
        g.add_option("--cact", help="Image from coil used for actual ASL acquisition (default: calibration image - only in longtr mode)", type="image")
        g.add_option("--isen", help="User-supplied sensitivity correction in ASL space")
        g.add_option("--senscorr-auto", "--senscorr", help="Apply automatic sensitivity correction using bias field from FAST", action="store_true", default=False)
        g.add_option("--senscorr-off", help="Do not apply any sensitivity correction", action="store_true", default=False)
        ret.append(g)

        g = IgnorableOptionGroup(parser, "Partial volume correction")
        g.add_option("--pvcorr", help="Apply partial volume correction", action="store_true", default=False)
        ret.append(g)

        return ret

def single_volume(wsp, img, moco=True, discard_first=True):
    """
    Convert a potentially 4D image into a single 3D volume

    :param moco: If True, perform basic motion correction
    :param discard_first: If True, discard first volume if nvols > 1

    """
    if img is not None:
        wsp.log.write(" - Pre-processing image: %s\n" % img.name)
        if img.ndim == 4:
            if discard_first and img.shape[3] > 1:
                wsp.log.write("   - Removing first volume to ensure data is in steady state\n")
                img = Image(img.data[..., :-1], header=img.header)
            
            if moco and img.shape[3] > 1:
                if moco:
                    wsp.log.write("   - Motion correcting\n")
                    img = fsl.mcflirt(img, out=fsl.LOAD, log=wsp.fsllog)["out"]
    
            wsp.log.write("   - Taking mean across time axis\n")
            img = Image(np.mean(img.data, axis=-1), header=img.header)

        return img
    else:
        return None

def get_cblip_correction(wsp):
    """
    Get the cblip based distortion correction warp

    Required workspace attributes
    -----------------------------

     - ``calib`` : Calibration image
     - ``cblip`` : Phase-encode-reversed calibration image
     - ``echospacing`` :
     - ``pedir`` : 

    Optional workspace attributes
    -----------------------------

    Updated workspace attributes
    ----------------------------

     - ``cblip_warp``    : CBLIP Distortion correction warp image
     
    """
    if wsp.topup is not None:
        return
    elif wsp.cblip is None:
        wsp.log.write("\nNo CBLIP images provided for distortion correction\n")
        return

    wsp.sub("topup")
    wsp.log.write("\nCalculating distortion Correction using TOPUP\n")

    topup_params = {
        "x"  : [[1, 0, 0, -99], [-1, 0, 0, -99]],
        "-x" : [[-1, 0, 0, -99], [1, 0, 0, -99]],
        "y"  : [[0, 1, 0, -99], [0, -1, 0, -99]],
        "-y" : [[0, -1, 0, -99], [0, 1, 0, -99]],
        "z"  : [[0, 0, 1, -99], [0, 0, -1, -99]],
        "-z" : [[0, 0, -1, -99], [0, 0, 1, -99]],
    }
    dim_idx = {
        "x"  : 0, "-x" : 0,
        "y"  : 1, "-y" : 1,
        "z"  : 2, "-z" : 2,
    }
    my_topup_params = np.array(topup_params[wsp.pedir], dtype=np.float)
    dimsize = wsp.asldata.shape[dim_idx[wsp.pedir]]
    my_topup_params[:, 3] = wsp.echospacing * (dimsize - 1)
    wsp.topup.params = my_topup_params
    
    # Run TOPUP to calculate correction
    wsp.topup.calib_blipped = Image(np.stack((wsp.calib.data, wsp.cblip.data), axis=-1), header=wsp.calib.header)
    topup_result = fsl.topup(imain=wsp.topup.calib_blipped, datain=wsp.topup.params, config="b02b0.cnf", out=fsl.LOAD, iout=fsl.LOAD, fout=fsl.LOAD, log=wsp.fsllog)
    wsp.topup.fieldcoef, wsp.topup.movpar = topup_result["out_fieldcoef"], topup_result["out_movpar"]
    wsp.topup.iout = topup_result["iout"]
    wsp.topup.fout = topup_result["fout"]

    page = wsp.report.page("topup")
    page.heading("TOPUP distortion correction", level=0)
    page.text("PE direction: %s" % wsp.pedir)
    page.text("Echo spacing: %f s" % wsp.echospacing)
    page.heading("Correction image", level=1)
    for dim in range(3):
        img = Image(wsp.topup.fieldcoef.data[..., dim], header=wsp.topup.fieldcoef.header)
        page.text("Dimension %i" % dim)
        page.image("fmap_warp%i" % dim, LightboxImage(img))

def get_fieldmap_correction(wsp):
    """
    Get the fieldmap based distortion correction warp

    Required workspace attributes
    -----------------------------

     - ``pwi``          : Perfusion weighted image (generated by preproc_asl)
     - ``fmap``         : Fieldmap image
     - ``fmapmag``      : Fieldmap magnitude image
     - ``fmapmagbrain`` : Fieldmap magnitude brain image
     - ``echospacing``  :
     - ``pedir``        :

    Optional workspace attributes
    -----------------------------

     - ``nofmapreg``       : If True assume fieldmap in structural space

    Updated workspace attributes
    ----------------------------

     - ``fmap_warp``    : Fieldmap distortion correction warp image in ASL space
    """
    if wsp.fieldmap is not None:
        return
    elif wsp.fmap is None or wsp.fmapmag is None or wsp.fmapmagbrain is None:
        wsp.log.write("\nNo fieldmap images for distortion correction\n")
        return
    elif wsp.pedir is None or wsp.echospacing is None:
        wsp.log.write("\nWARNING: Fieldmap images supplied but pedir and echospacing required for distortion correction\n")
        return

    struc.segment(wsp)
    wsp.sub("fieldmap")
    wsp.log.write("\nCalculating distortion correction from fieldmap images using EPI_REG\n")

    epi_reg_opts = {
        "inweight" : wsp.inweight,
        "init" : wsp.reg.asl2struc,
        "fmap" : wsp.fmap,
        "fmapmag" : wsp.fmapmag,
        "fmapmagbrain" : wsp.fmapmagbrain,
        "pedir" : wsp.pedir,
        "echospacing" : wsp.echospacing,
        "nofmapreg" : wsp.ifnone("nofmapreg", False),
    }
    
    result = epi_reg(epi=wsp.asldata.perf_weighted(), t1=wsp.structural.struc, t1brain=wsp.structural.brain, out=fsl.LOAD, wmseg=wsp.structural.wm_seg, log=wsp.fsllog, **epi_reg_opts)
    wsp.fieldmap.warp_struc = result["out_warp"]
    wsp.fieldmap.asl2struc = result["out"]
    wsp.fieldmap.struc2asl = np.linalg.inv(wsp.fieldmap.asl2struc)

    result = fsl.convertwarp(out=fsl.LOAD, ref=wsp.nativeref, warp1=wsp.fieldmap.warp_struc, postmat=wsp.fieldmap.struc2asl, rel=True, log=wsp.fsllog)
    wsp.fieldmap.warp = result["out"]
        
    page = wsp.report.page("fmap")
    page.heading("Fieldmap distortion correction", level=0)
    page.text("PE direction: %s" % wsp.pedir)
    page.text("Echo spacing: %f s" % wsp.echospacing)
    page.heading("Correction warps", level=1)
    for dim in range(3):
        img = Image(wsp.fieldmap.warp.data[..., dim], header=wsp.fieldmap.warp.header)
        page.text("Dimension %i" % dim)
        page.image("fmap_warp%i" % dim, LightboxImage(img))

def get_motion_correction(wsp):
    """
    Calculate motion correction transforms for ASL data
    
    Note simple motion correction of multi-volume calibration data is done in preprocessing.

    The reference volume for motion correction is the calibration image, if supplied, or
    otherwise the middle volume of the ASL data is used. 
    
    If the calibration image is used, the inverse of the middle ASL volume -> calibration
    transform is applied to each transform matrix. This ensures that the middle volume of 
    the ASL data is unchanged and interpolation on the other volumes is also minimised.
    In this case, all calibration images are also themselves transformed to bring them in
    to ASL middle volume space.

    Required workspace attributes
    -----------------------------

     - ``asldata`` : ASL data image

    Optional workspace attributes
    -----------------------------

     - ``calib``    : Calibration image

    Updated workspace attributes
    ----------------------------

     - ``asldata_mc_mats`` : Sequence of matrices giving motion correction transform for each ASL volume
     - ``asl2calib``       : ASL->calibration image transformation
     - ``calib2asl``       : Calibration->ASL image transformation
    """
    if wsp.moco is not None:
        return
    elif not wsp.mc:
        wsp.log.write("\nNo motion correction\n")
        return

    reg.init(wsp)

    wsp.sub("moco")
    wsp.log.write("\nCalculating Motion Correction\n")
    # If available, use the calibration image as reference since this will be most consistent if the data has a range 
    # of different TIs and background suppression etc. This also removes motion effects between asldata and calibration image
    if wsp.input.regfrom is not None:
        wsp.log.write(" - Using user-specified regfrom as reference\n")
        ref_source = "User specified: %s" % (wsp.input.regfrom.name)
        mcflirt_result = fsl.mcflirt(wsp.input.asldata, reffile=wsp.input.regfrom, out=fsl.LOAD, mats=fsl.LOAD, log=wsp.fsllog)
        mats = [mcflirt_result["out.mat/MAT_%04i" % vol] for vol in range(wsp.asldata.shape[3])]
    elif wsp.input.calib is not None:
        wsp.log.write(" - Using calibration image as reference\n")
        ref_source = "Calibration image"
        wsp.moco.ref = single_volume(wsp, wsp.input.calib)
        wsp.moco.input = wsp.input.asldata
        mcflirt_result = fsl.mcflirt(wsp.moco.input, reffile=wsp.moco.ref, out=fsl.LOAD, mats=fsl.LOAD, log=wsp.fsllog)
        mats = [mcflirt_result["out.mat/MAT_%04i" % vol] for vol in range(wsp.asldata.shape[3])]

        # To reduce interpolation of the ASL data change the transformations so that we end up in the space of the central volume of asldata
        wsp.reg.asl2calib = mats[int(float(len(mats))/2)]
        wsp.reg.calib2asl = np.linalg.inv(wsp.reg.asl2calib)
        mats = [np.dot(wsp.reg.calib2asl, mat) for mat in mats]
        
        wsp.log.write("   ASL middle volume->Calib:\n%s\n" % str(wsp.reg.asl2calib))
        wsp.log.write("   Calib->ASL middle volume:\n%s\n" % str(wsp.reg.calib2asl))
    else:
        wsp.log.write(" - Using ASL data middle volume as reference\n")
        ref_source = "ASL data middle volume: %i" % int(float(wsp.asldata.shape[3])/2)
        mcflirt_result = fsl.mcflirt(wsp.input.asldata, out=fsl.LOAD, mats=fsl.LOAD, log=wsp.fsllog)
        mats = [mcflirt_result["out.mat/MAT_%04i" % vol] for vol in range(wsp.asldata.shape[3])]
        
    # Convert motion correction matrices into single (4*nvols, 4) matrix - convenient for writing
    # to file, and same form that applywarp expects
    wsp.moco.mc_mats = np.concatenate(mats, axis=0)

    page = wsp.report.page("moco")
    page.heading("Motion correction", level=0)
    page.heading("Reference volume", level=1)
    page.text(ref_source)
    page.heading("Motion parameters", level=1)
    moco_params = [reg.get_motion_params(mat) for mat in mats]
    trans = [p[0] for p in moco_params]
    abstrans = np.fabs(trans)
    rot = [p[1] for p in moco_params]
    absrot = np.fabs(rot)
    page.table([
        ["Mean translation", "%.3g mm" % np.mean(trans)],
        ["Translation std.dev.", "%.3g mm" % np.std(trans)],
        ["Absolute maximum translation", "%.3g mm (volume %i)" % (np.max(abstrans), np.argmax(abstrans))],
        ["Mean rotation", "%.3g \N{DEGREE SIGN}" % np.mean(rot)],
        ["Rotation std.dev.", "%.3g \N{DEGREE SIGN}" % np.std(rot)],
        ["Absolute maximum rotation", "%.3g \N{DEGREE SIGN} (volume %i)" % (np.max(absrot), np.argmax(absrot))],       
    ])
    page.image("moco_trans", LineGraph(trans, "Volume number", "Translation (mm)"))
    page.image("moco_rot", LineGraph(rot, "Volume number", "Rotation relative to reference (\N{DEGREE SIGN})"))

def get_sensitivity_correction(wsp):
    """
    Get sensitivity correction image
    
    Required workspace attributes
    -----------------------------

     - ``asldata`` : ASL data

    Optional workspace attributes
    -----------------------------

     - ``isen`` : User supplied sensitivity image
     - ``cact`` : Calibration image. Used in conjunction with ``cref`` to calculate sensitivity map
     - ``calib`` : Calibration image. Used as alternative to cact provided ``mode`` is ``longtr``
     - ``cref`` : Calibration reference image 
     - ``senscorr_auto`` : If True, automatically calculate sensitivity correction using FAST
     - ``senscorr_off`` If True, do not apply sensitivity correction

    Updated workspace attributes
    ----------------------------

     - ``sensitivity``    : Sensitivity correction image in ASL space
    """
    if wsp.senscorr is not None:
        return

    wsp.log.write("\nCalculating Sensitivity correction\n")
    sensitivity = None
    bias = None
    if wsp.senscorr_off:
        wsp.log.write(" - Sensitivity correction disabled\n")
    elif wsp.isen is not None:
        wsp.log.write(" - Sensitivity image supplied by user\n")
        sensitivity = wsp.isen
    elif wsp.cact is not None and wsp.cref is not None:
        wsp.log.write(" - Sensitivity image calculated from calibration actual and reference images\n")
        cref_data = np.copy(wsp.cref.data)
        cref_data[cref_data == 0] = 1
        sensitivity = Image(wsp.cact.data.astype(np.float) / cref_data, header=wsp.calib.header)
    elif wsp.calib is not None and wsp.cref is not None:
        if wsp.ifnone("mode", "longtr") != "longtr":
            raise ValueError("Calibration reference image specified but calibration image was not in longtr mode - need to provided additional calibration image using the ASL coil")
        wsp.log.write(" - Sensitivity image calculated from calibration and reference images\n")
        cref_data = np.copy(wsp.cref.data)
        cref_data[cref_data == 0] = 1
        sensitivity = Image(wsp.calib.data.astype(np.float) / cref_data, header=wsp.calib.header)
    elif wsp.senscorr_auto and wsp.structural.bias is not None:
        struc.segment(wsp)
        wsp.log.write(" - Sensitivity image calculated from bias field\n")
        bias = reg.struc2asl(wsp, wsp.structural.bias)
        sensitivity = Image(np.reciprocal(bias.data), header=bias.header)
    else:
        wsp.log.write(" - No source of sensitivity correction was found\n")

    if sensitivity is not None:
        sdata = sensitivity.data
        sdata[sdata < 1e-12] = 1
        sdata[np.isnan(sdata)] = 1
        sdata[np.isinf(sdata)] = 1
        wsp.sub("senscorr")
        wsp.senscorr.sensitivity = Image(sdata, header=sensitivity.header)

        page = wsp.report.page("sensitivity")
        page.heading("Sensitivity correction", level=0)
        page.heading("Sensitivity map", level=1)
        page.image("sensitivity", LightboxImage(wsp.senscorr.sensitivity))

    if bias is not None:
        wsp.senscorr.bias = bias

def apply_corrections(wsp):
    """
    Apply distortion and motion corrections to ASL and calibration data

    Required workspace attributes
    -----------------------------

     - ``asldata_orig`` : Uncorrected ASL data image

    Optional workspace attributes
    -----------------------------

     - ``calib_orig``      : Calibration image
     - ``cref_orig``       : Calibration reference image
     - ``cblip_orig``      : Calibration BLIP image
     - ``asldata_mc_mats`` : ASL motion correction matrices
     - ``calib2asl``       : Calibration -> ASL transformation matrix
     - ``distcorr_warp``   : Distortion correction warp image
     - ``gdc_warp``        : Gradient distortion correction warp image

    Updated workspace attributes
    ----------------------------

     - ``asldata``    : Corrected ASL data
     - ``calib``      : Corrected calibration image
     - ``cref``       : Corrected calibration reference image
     - ``cblip``      : Corrected calibration BLIP image
    """
    wsp.log.write("\nApplying preprocessing corrections\n")
    if wsp.corrected is None:
        wsp.sub("corrected")

    wsp.corrected.asldata = wsp.input.asldata
    wsp.corrected.nativeref = wsp.input.asldata.mean()
    wsp.corrected.calib = single_volume(wsp, wsp.input.calib)
    wsp.corrected.cref = single_volume(wsp, wsp.input.cref)
    wsp.corrected.cact = single_volume(wsp, wsp.input.cact)
    wsp.corrected.cblip = single_volume(wsp, wsp.input.cblip)
    
    wsp.log.write(" - Data transformations\n")
    if wsp.moco is not None:
        wsp.log.write("   - Using motion correction\n")

    warps, moco_mats = [], None
    if wsp.fieldmap is not None:
        wsp.log.write("   - Using fieldmap distortion correction\n")
        warps.append(wsp.fieldmap.warp)
    
    if wsp.gdc_warp:
        wsp.log.write("   - Using user-supplied GDC warp\n")
        warps.append(wsp.gdc_warp)
        
    if wsp.moco is not None: 
        moco_mats = wsp.moco.mc_mats

    if warps:
        kwargs = {}
        for idx, warp in enumerate(warps):
            kwargs["warp%i" % (idx+1)] = warp
                
        wsp.log.write("   - Converting all warps to single transform and extracting Jacobian\n")
        result = fsl.convertwarp(ref=wsp.nativeref, out=fsl.LOAD, rel=True, jacobian=fsl.LOAD, log=wsp.fsllog, **kwargs)
        wsp.corrected.total_warp = result["out"]

        # Calculation of the jacobian for the warp - method suggested in:
        # https://www.jiscmail.ac.uk/cgi-bin/webadmin?A2=FSL;d3fee1e5.0908
        wsp.corrected.warp_coef = fnirtfileutils(wsp.corrected.total_warp, outformat="spline", out=fsl.LOAD, log=wsp.fsllog)["out"]
        jacobian = fnirtfileutils(wsp.corrected.warp_coef, jac=fsl.LOAD, log=wsp.fsllog)["jac"]
        wsp.corrected.jacobian = Image(jacobian.data, header=wsp.corrected.total_warp.header)
        
    if not warps and moco_mats is None:
        wsp.log.write("   - No corrections to apply\n")
    else:
        # Apply all corrections to ASL data - note that we make sure the output keeps all the ASL metadata
        wsp.log.write("   - Applying to ASL data\n")
        asldata_corr = correct_img(wsp, wsp.input.asldata, moco_mats)
        wsp.corrected.asldata = wsp.input.asldata.derived(asldata_corr.data)

        # Apply corrections to calibration images
        if wsp.input.calib is not None:
            wsp.log.write("   - Applying to calibration data\n")
            wsp.corrected.calib = correct_img(wsp, wsp.corrected.calib, wsp.reg.calib2asl)
        
            if wsp.cref is not None:
                wsp.corrected.cref = correct_img(wsp, wsp.corrected.cref, wsp.reg.calib2asl)
            if wsp.cblip is not None:
                wsp.corrected.cblip = correct_img(wsp, wsp.corrected.cblip, wsp.reg.calib2asl)

    if wsp.topup is not None:
        wsp.log.write(" - Adding TOPUP distortion correction\n")
        # This can't currently be done using the FSL wrappers - we need the TOPUP output as two prefixed files
        # Only workaround currently is to create a temp directory to store appropriately named input files
        topup_input = tempfile.mkdtemp(prefix="topup_input")
        try:
            wsp.topup.fieldcoef.save("%s/topup_fieldcoef" % topup_input)
            movpar_file = open("%s/topup_movpar.txt" % topup_input, "w")
            for row in wsp.topup.movpar:
                movpar_file.write("\t".join([str(val) for val in row]) + "\n")
            movpar_file.close()
            # TOPUP does not do the jacboian magntiude correction - so only okay if using voxelwise calibration
            wsp.corrected.calib = fsl.applytopup(wsp.corrected.calib, datain=wsp.topup.params, index=1, topup="%s/topup" % topup_input, out=fsl.LOAD, method="jac", log=wsp.fsllog)["out"]
            wsp.corrected.cblip = fsl.applytopup(wsp.corrected.cblip, datain=wsp.topup.params, index=2, topup="%s/topup" % topup_input, out=fsl.LOAD, method="jac", log=wsp.fsllog)["out"]
            if wsp.cref:
                wsp.corrected.cref = fsl.applytopup(wsp.corrected.cref, datain=wsp.topup.params, index=1, topup="%s/topup" % topup_input, out=fsl.LOAD, method="jac", log=wsp.fsllog)["out"]
            post_topup = fsl.applytopup(wsp.corrected.asldata, datain=wsp.topup.params, index=1, topup="%s/topup" % topup_input, out=fsl.LOAD, method="jac", log=wsp.fsllog)["out"]
            wsp.corrected.asldata = wsp.corrected.asldata.derived(post_topup.data)
            #if wsp.calib_method != "voxel":
            #    wsp.log.write("WARNING: Using TOPUP does not correct for magntiude using the jocbian in distortion correction")
            #    wsp.log.write("         This is not optimal when not using voxelwise calibration\n")
            #    wsp.log.write("         To avoid this supply structural image(s)\n")
        finally:
            shutil.rmtree(topup_input)

    if wsp.senscorr and wsp.corrected.calib:
        # Apply sensitivity correction to calibration image only. In principle we could
        # apply it to the ASL image, but in keeping with OXFORD_ASL we apply it to the 
        # perfusion maps instead at output time. Note that this means the sensitivity
        # correction cancels out of the calibrated outputs when using voxelwise calibration
        wsp.corrected.calib, = apply_sensitivity_correction(wsp, wsp.corrected.calib)

def correct_img(wsp, img, linear_mat):
    """
    Apply combined warp/linear transformations to an image
    
    :param img: fsl.data.image.Image to correct
    :param linear_mat: img->ASL space linear transformation matrix.
    :return: Corrected Image

    If a jacobian is present, also corrects for quantitative signal magnitude as volume has been locally scaled

    FIXME there are slight differences to oxford_asl here due to use of spline interpolation rather than
    applyxfm4D which uses sinc interpolation.

    Required workspace attributesd
    -----------------------------

     - ``asldata_mean`` : Mean ASL image used as reference space

    Optional workspace attributes
    -----------------------------

     - ``total_warp``      : Combined warp image
     - ``jacobian``        : Jacobian associated with warp image
    """
    if wsp.corrected.total_warp is not None:
        img = reg.transform(wsp, img, trans=wsp.corrected.total_warp, ref=wsp.nativeref, premat=linear_mat)
    else:
        img = reg.transform(wsp, img, trans=linear_mat, ref=wsp.nativeref)

    if wsp.corrected.total_jacobian is not None:
        wsp.log.write("   - Correcting for local volume scaling using Jacobian\n")
        jdata = wsp.corrected.total_jacobian.data
        if img.data.ndim == 4:
            # Required to make broadcasting work
            jdata = jdata[..., np.newaxis]
        img = Image(img.data * jdata, header=img.header)
    return img

def apply_sensitivity_correction(wsp, *imgs):
    """
    Apply sensitivity correction

    :param imgs: Sequence of Image objects

    :return: Tuple of corrected Image objects corresponding to input.
             If no sensitivity correction is defined, returns the same
             images as input.

    Optional workspace attributes
    -----------------------------

     - ``sensitivity``  : Sensitivity correction image
     - ``senscorr_off`` : If True, no correction will be applied even if ``sensitivity`` image exists
    """
    if wsp.senscorr is not None:
        wsp.log.write(" - Applying sensitivity correction\n")
        ret = []
        for img in imgs:
            ret.append(Image(img.data / wsp.senscorr.sensitivity.data, header=img.header))
        return tuple(ret)
    else:
        return tuple(imgs)
