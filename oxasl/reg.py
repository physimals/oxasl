#!/bin/env python
"""
OXASL - Registration for ASL data

Copyright (c) 2008 University of Nottingham
"""
from __future__ import unicode_literals

import os
import sys
import math
import warnings
import traceback

import numpy as np

from fsl.data.image import Image, defaultExt
import fsl.wrappers as fsl

try:
    import regtricks
except ImportError:
    regtricks = None

from oxasl import __version__, Workspace, struc, brain
from oxasl.options import AslOptionParser, GenericOptions, OptionCategory, OptionGroup, load_matrix
from oxasl.wrappers import epi_reg
from oxasl.reporting import LightboxImage

class Options(OptionCategory):
    """
    OptionCategory which contains options for registration of ASL data to structural image
    """

    def __init__(self, **kwargs):
        OptionCategory.__init__(self, "reg", **kwargs)

    def groups(self, parser):
        groups = []

        group = OptionGroup(parser, "Registration")
        group.add_option("--aslref", "--nativeref", "--regfrom", help="ASL Registration image (e.g. perfusion weighted image)", type="image")
        group.add_option("--aslref-method", "--nativeref-method", "--regfrom-method", help="How to choose the ASL registration reference image - calib=use calibration image, mean=use mean ASL data, pwi=use PWI (mean differenced ASL data)")
        group.add_option("--reg-init-bbr", help="Use BBR for initial as well as final registration - generally used in conjunction with --aslref-method=pwi", action="store_true", default=False)
        group.add_option("--struc-std-fnirt", help="Use FNIRT to do nonlinear structural/std space registration when FSL_ANAT output is not available", action="store_true", default=False)
        #group.add_option("--bbr", dest="do_bbr", help="Include BBR registration step using EPI_REG", action="store_true", default=False)
        #group.add_option("--flirt", dest="do_flirt", help="Include rigid-body registration step using FLIRT", action="store_true", default=True)
        #group.add_option("--flirtsch", help="user-specified FLIRT schedule for registration")
        groups.append(group)

        #group = OptionGroup(parser, "Extra BBR registration refinement")
        #group.add_option("-c", dest="cfile", help="ASL control/calibration image for initial registration - brain extracted")
        #group.add_option("--wm_seg", dest="wm_seg", help="tissue segmenation image for bbr (in structural image space)")
        #groups.append(group)

        #group = OptionGroup(parser, "Distortion correction using fieldmap (see epi_reg)")
        #g.add_option("--nofmapreg", dest="nofmapreg", help="do not perform registration of fmap to T1 (use if fmap already registered)", action="store_true", default=False)
        #groups.append(group)

        #group = OptionGroup(parser, "Deprecated")
        #g.add_option("-r", dest="lowstruc", help="extra low resolution structural image - brain extracted")
        #g.add_option("--inweight", dest="inweight", help="specify weights for input image - same functionality as the flirt -inweight option", type="float")
        #groups.append(group)

        return groups

def run(wsp, redo=False, struc_flirt=True, struc_bbr=False, use_quantification_wsp=None):
    if wsp.reg is None:
        wsp.sub("reg")
        redo = True

    if redo:
        wsp.log.write("\nPerforming registration\n")

        if not struc_bbr and wsp.reg_init_bbr:
            wsp.log.write(" - Using BBR for initial structural registration\n")
            struc_bbr = True

        # Save any previous registration data
        if wsp.reg.aslref is not None:
            idx = 1
            while getattr(wsp.reg, "aslref_old_%i" % idx) is not None:
                idx += 1
            for attr in ("aslref", "asl2struc", "struc2asl"):
                setattr(wsp.reg, "%s_old_%i" % (attr, idx), getattr(wsp.reg, attr))

        get_ref_imgs(wsp, use_quantification_wsp)
        reg_asl2calib(wsp)
        reg_asl2struc(wsp, flirt=struc_flirt, bbr=struc_bbr)
        reg_asl2custom(wsp)

def get_ref_imgs(wsp, use_quantification_wsp=None):
    """
    Get the images that define the various processing 'spaces' and are used for registration
    to/from these spaces. The built in spaces are 'asl' (aka 'native') 'struc' and 'std' (MNI).

    Note that the 'custom' space requires a user-specified reference image and transformation
    from structural space

    aslref defines the 'asl' space

    Optional workspace attributes
    -----------------------------

     - ``aslref`` : User-supplied registration reference image
     - ``aslref_method`` : Method for choosing registration reference image
     - ``asldata`` : Raw ASL data
     - ``calib``   : Calibration image

    Updated workspace attributes
    ----------------------------

     - ``aslref``    : Registration reference image in ASL space
     - ``strucref``     : Registration reference image in structural space
     - ``stdref``       : Registration reference image in standard space
    """
    if wsp.input is not None and wsp.input.aslref_method is None:
        if wsp.asldata is not None and wsp.asldata.iaf in ("tc", "ct"):
            wsp.reg.aslref_method = "mean"
        elif wsp.calib is not None and wsp.asldata is not None and wsp.calib.sameSpace(wsp.asldata):
            wsp.reg.aslref_method = "calib"
        else:
            wsp.reg.aslref_method = "mean"

    if use_quantification_wsp:
        wsp.log.write(" - ASL Registration reference image is PWI image generated by quantification\n")
        pwi = np.copy(use_quantification_wsp.finalstep.mean_ftiss.data)
        pwi[pwi < 0] = 0
        wsp.reg.aslref = Image(pwi, header=use_quantification_wsp.finalstep.mean_ftiss.header)
    elif wsp.input.aslref is not None:
        wsp.log.write(" - ASL Registration reference image supplied by user\n")
        wsp.reg.aslref = wsp.input.aslref
    elif wsp.reg.aslref_method == "mean":
        wsp.log.write(" - ASL Registration reference is mean ASL signal (brain extracted)\n")
        wsp.reg.aslref = brain.brain(wsp, wsp.asldata.mean(), thresh=0.2)
    elif wsp.reg.aslref_method == "calib":
        wsp.log.write(" - ASL Registration reference is calibration image (brain extracted)\n")
        if not wsp.calib.sameSpace(wsp.asldata):
            raise ValueError("Calibration image is not in same space as ASL data - cannot use as registration reference")
        wsp.reg.aslref = brain.brain(wsp, wsp.calib, thresh=0.2)
    elif wsp.reg.aslref_method == "pwi":
        wsp.log.write(" - ASL Registration reference is PWI (brain extracted)\n")
        wsp.reg.aslref = brain.brain(wsp, wsp.asldata.perf_weighted(), thresh=0.2)
    else:
        raise ValueError("Unrecognized aslref_method: %s" % wsp.aslref_method)

    if wsp.input is not None: wsp.reg.calibref = wsp.input.calib
    if wsp.structural is not None: wsp.reg.strucref = wsp.structural.struc
    wsp.reg.stdref = Image(os.path.join(os.environ["FSLDIR"], "data/standard/MNI152_T1_2mm_brain"))

def reg_asl2calib(wsp):
    """
    Register calibration image to ASL space

    Note that this might already have been done as part of motion correction
    """
    if wsp.calib_aslreg:
        wsp.log.write(" - Calibration image already registered to ASL image\n")
    elif wsp.moco is not None and wsp.moco.asl2calib is not None:
        wsp.log.write(" - Calibration image registered to ASL image as part of motion correction\n")
        wsp.reg.asl2calib = wsp.moco.asl2calib
        wsp.reg.calib2asl = wsp.moco.calib2asl
    elif wsp.calib is not None:
        wsp.log.write(" - Registering calibration image to ASL image\n")
        _, wsp.reg.asl2calib = reg_flirt(wsp, wsp.reg.aslref, wsp.calib)
        wsp.reg.calib2asl = np.linalg.inv(wsp.reg.asl2calib)

    if wsp.reg.asl2calib is not None:
        wsp.log.write(" - ASL->Calibration transform\n")
        wsp.log.write(str(wsp.reg.asl2calib) + "\n")
        wsp.log.write(" - Calibration->ASL transform\n")
        wsp.log.write(str(wsp.reg.calib2asl) + "\n")

def reg_asl2custom(wsp):
    """
    Register custom output image to ASL space, via structural. 

    If no output_custom_mat (struc -> custom) has been provided, then
    FLIRT will be used to generate this. The transformation from ASL space
    is the concatenation of asl2struc and struc2custom. 
    """
    if wsp.input.output_custom:
        wsp.reg.customref = wsp.input.output_custom

        if wsp.input.output_custom_mat is not None:
            wsp.reg.struc2custom = np.loadtxt(wsp.input.output_custom_mat)
            wsp.reg.custom2struc = np.linalg.inv(wsp.reg.struc2custom)

        else:  
            wsp.log.write(" - Registering calibration image to ASL image\n")
            _, wsp.reg.custom2struc = reg_flirt(wsp, wsp.reg.customref, wsp.structural.struc)
            wsp.reg.struc2custom = np.linalg.inv(wsp.reg.custom2struc)

        wsp.reg.asl2custom = wsp.reg.struc2custom @ wsp.reg.asl2struc
        wsp.reg.custom2asl = np.linalg.inv(wsp.reg.asl2custom)

def reg_asl2struc(wsp, flirt=True, bbr=False, name="initial"):
    """
    Registration of ASL images to structural image

    :param flirt: If provided, sets whether to use FLIRT registration
    :param bbr: If provided, sets whether to use BBR registration

    Required workspace attributes
    -----------------------------

     - ``aslref``            : Registration reference image in ASL space
     - ``struc``              : Structural image

    Updated workspace attributes
    ----------------------------

     - ``asl2struc``    : ASL->structural transformation matrix
     - ``struc2asl``    : Structural->ASL transformation matrix
     - ``regto``        : ``aslref`` image transformed to structural space
    """
    if wsp.structural is not None and wsp.structural.struc is not None:
        if wsp.struc2asl is not None or wsp.asl2struc is not None:
            wsp.log.write("\nASL->Structural registration provided by user\n")

            wsp.reg.struc2asl = wsp.struc2asl
            wsp.reg.asl2struc = wsp.asl2struc
            if wsp.reg.asl2struc is None:
                wsp.reg.asl2struc = np.linalg.inv(wsp.reg.struc2asl)
            if wsp.reg.struc2asl is None:
                wsp.reg.struc2asl = np.linalg.inv(wsp.reg.asl2struc)
            wsp.reg.regto = change_space(wsp, wsp.reg.aslref, "struc")
        else:
            wsp.log.write("\nRegistering ASL data to structural data\n")
            if flirt:
                wsp.reg.regto, wsp.reg.asl2struc = reg_flirt(wsp, wsp.reg.aslref, wsp.structural.brain, wsp.reg.asl2struc)
            if bbr:
                wsp.reg.regto, wsp.reg.asl2struc = reg_bbr(wsp)

            wsp.reg.struc2asl = np.linalg.inv(wsp.reg.asl2struc)

        wsp.log.write(" - ASL->Structural transform\n")
        wsp.log.write(str(wsp.reg.asl2struc) + "\n")
        wsp.log.write(" - Structural->ASL transform\n")
        wsp.log.write(str(wsp.reg.struc2asl) + "\n")

        name = "final" if bbr else "initial"
        page = wsp.report.page("asl2struc_%s" % name)
        page.heading("%s ASL -> Structural registration" % name.title(), level=0)
        page.heading("Transformation parameters", level=1)
        motion_params = get_transform_params(wsp.reg.asl2struc)
        page.table([
            ["Translation magnitude", "%.3g mm" % motion_params[0]],
            ["Rotation magnitude", "%.3g \N{DEGREE SIGN}" % motion_params[1]],
        ])
        page.heading("ASL->Structural transformation matrix", level=1)
        page.matrix(wsp.reg.asl2struc)
        page.heading("Structural->ASL transformation matrix", level=1)
        page.matrix(wsp.reg.struc2asl)

        if wsp.structural.gm_seg is not None:
            gm_asl = change_space(wsp, wsp.structural.gm_seg, "asl", interp="nn")
            page.heading("GM mask aligned with ASL data", level=1)
            page.image("gm_reg_%s" % name, LightboxImage(gm_asl, bgimage=wsp.reg.aslref))
            wm_asl = change_space(wsp, wsp.structural.wm_seg, "asl", interp="nn")
            page.heading("WM mask aligned with ASL data", level=1)
            page.image("wm_reg_%s" % name, LightboxImage(wm_asl, bgimage=wsp.reg.aslref))

def reg_struc2std(wsp, **kwargs):
    """
    Determine structural -> standard space registration

    Optional workspace attributes
    -----------------------------

     - ``structural.struc``   : Structural image
     - ``fslanat``            : Path to existing FSLANAT data

    Updated workspace attributes
    ----------------------------

     - ``reg.struc2std``    : Structural->MNI transformation matrix - either warp image or FLIRT matrix
     - ``reg.std2struc``    : MNI->structural transformation - either warp image or FLIRT matrix
    """
    if wsp.reg.std2struc is not None:
        return

    if wsp.struc2std_warp is not None:
        wsp.log.write(" - Using user-specified structural->std nonlinear transformation warp\n")
        wsp.reg.struc2std = wsp.struc2std_warp

    elif wsp.struc2std is not None:
        wsp.log.write(" - Using user-specified structural->std linear transformation matrix\n")
        wsp.reg.struc2std = wsp.struc2std
        wsp.log.write(str(wsp.reg.struc2std) + "\n")

    elif wsp.fslanat:
        warp = os.path.join(wsp.fslanat, "T1_to_MNI_nonlin_coeff.nii.gz")
        mat = os.path.join(wsp.fslanat, "T1_to_MNI_lin.mat")
        if os.path.isfile(warp):
            wsp.log.write(" - Using structural->std nonlinear transformation from FSL_ANAT\n")
            wsp.reg.struc2std = Image(warp, loadData=False)
        elif os.path.isfile(mat):
            wsp.log.write(" - Using structural->std linear transformation from FSL_ANAT\n")
            wsp.reg.struc2std = load_matrix(mat)
            wsp.log.write(str(wsp.reg.struc2std) + "\n")

    if wsp.reg.struc2std is None:
        wsp.log.write(" - Registering structural image to standard space using FLIRT\n")
        flirt_result = fsl.flirt(wsp.structural.brain, os.path.join(os.environ["FSLDIR"], "data/standard/MNI152_T1_2mm_brain"), omat=fsl.LOAD)
        wsp.reg.struc2std = flirt_result["omat"]

        if wsp.struc_std_fnirt:
            wsp.log.write(" - Registering structural image to standard space using FNIRT\n")
            fnirt_result = fsl.fnirt(wsp.structural.brain, aff=wsp.reg.struc2std, config="T1_2_MNI152_2mm.cnf", cout=fsl.LOAD)
            wsp.reg.struc2std = fnirt_result["cout"]

    if isinstance(wsp.reg.struc2std, Image):
        # Calculate the inverse warp using INVWARP
        invwarp_result = fsl.invwarp(wsp.reg.struc2std, wsp.structural.struc, out=fsl.LOAD)
        wsp.reg.std2struc = invwarp_result["out"]
    else:
        wsp.reg.std2struc = np.linalg.inv(wsp.reg.struc2std)

def get_img_space(wsp, img):
    """
    Find out what image space an image is in
    
    Note that this only compares the voxel->world transformation matrix to the
    reference image for each space. It is quite possible for two images to be in
    the same space but not be registered to one another. In this case, 
    the returned space may not be accurate when determining whether a registration
    is required.

    :param wsp: Workspace object
    :param img: Image
    :return: Name of image space for ``img``, e.g. ``asl``, ``struc``
    """ 
    img_space = None
    for space in ('asl', 'calib', 'struc', 'std', 'custom'):
        ref = getattr(wsp.reg, "%sref" % space)
        if ref is not None and img.sameSpace(ref):
            img_space = space
            break

    if img_space is None:
        raise RuntimeError("Could not determine space for image: %s" % str(img))
    return img_space

def change_space(wsp, img, target_space, source_space=None, **kwargs):
    """
    Convert an image to a different space
    
    Note that while the source space can be determined from the image, this may
    not be correct if images (e.g. ASL and calibration) share the same voxel->world
    transformation but still need registration to one another

    :param wsp: Workspace object
    :param img: Image
    :param target_space: Either an Image in the target space, or the name of the target space
    :param src_space: If specified, explicit indication of source image space
    """
    # Source and target space can be specified directly or determined from image
    if source_space is None:
        source_space = get_img_space(wsp, img)
    if isinstance(target_space, Image):
        target_space = get_img_space(wsp, target_space)

    # Backwards compatiblity for naming convention
    if source_space == "native":
        source_space = "asl"
    if target_space == "native":
        target_space = "asl"

    target_ref = getattr(wsp.reg, "%sref" % target_space)
    if target_ref is None:
        raise RuntimeError("Couldn't find reference image for target space: %s" % target_space)

    if source_space == target_space:
        # Nothing to be done - image is already in target space
        return img

    if source_space == "std" or target_space == "std":
        # Calculating the nonlinear standard space registration is slow so we only do
        # it if necessary
        reg_struc2std(wsp, **kwargs)

    # For ASL to/from std space, go via structural image using a pre/post matrix
    # Not necessary for custom space because we already have asl2custom and custom2asl
    # from matrix multiplication
    if source_space == "asl" and target_space == "std":
        tform = wsp.reg.struc2std
        kwargs["premat"] = wsp.reg.asl2struc
    elif source_space == "std" and target_space == "asl":
        tform = wsp.reg.std2struc
        kwargs["postmat"] = wsp.reg.struc2asl
    else:
        tform = getattr(wsp.reg, "%s2%s" % (source_space, target_space))

    if tform is None:
        raise RuntimeError("No registration available for transform %s->%s" % (source_space, target_space))

    return transform(wsp, img, tform, target_ref, **kwargs)

def transform(wsp, img, trans, ref, use_flirt=False, interp="trilinear", paddingsize=1, premat=None, postmat=None, mask=False, mask_thresh=0.5):
    """
    Transform an image

    :param wsp: Workspace, used for logging only
    :param img: Image to transform
    :param trans: Transformation matrix or warp image
    :param ref: Target space reference image
    :param use_flirt: Use flirt to apply the transformation which must be a matrix
    :param interp: Interpolation method
    :param paddingsize: Padding size in pixels
    :param premat: If trans is a warp, this can be set to a pre-warp affine transformation matrix

    :return: Transformed Image object
    """
    if trans is None:
        raise ValueError("Transformation matrix not available - has registration been performed?")

    have_warp = isinstance(trans, Image)
    is_moco = premat is not None and max(premat.shape) > 4
    if use_flirt and have_warp:
        raise ValueError("Cannot transform using Flirt when we have a warp")

    if regtricks is not None and have_warp and not is_moco:
        # For now we only use regtricks for nonlinear transformations. It could be used for any but doing so leads to significant
        # numerical differences compared to previous releases. Part of the reason may be some artifacts in the regtricks output
        return transform_regtricks(wsp, img, trans, ref, use_flirt=use_flirt, interp=interp, premat=premat, postmat=postmat, mask=mask, mask_thresh=mask_thresh)

    if use_flirt:
        if interp == "nn":
            interp = "nearestneighbour"
        ret = fsl.applyxfm(img, ref, trans, out=fsl.LOAD, interp=interp, paddingsize=paddingsize, log=wsp.fsllog)["out"]
    else:
        if have_warp:
            kwargs = {"warp" : trans, "premat" : premat, "rel" : True, "postmat" : postmat}
        else:
            kwargs = {"premat" : trans}
            if premat is not None:
                kwargs["premat"] = np.dot(kwargs["premat"], premat)
            elif postmat is not None:
                kwargs["premat"] = np.dot(postmat, kwargs["premat"])
        ret = fsl.applywarp(img, ref, out=fsl.LOAD, interp=interp, paddingsize=paddingsize, super=True, superlevel="a", log=wsp.fsllog, **kwargs)["out"]

    if mask:
        # Binarise mask images
        ret = Image((ret.data > mask_thresh).astype(np.int32), header=ret.header)
    return ret

def transform_regtricks(wsp, img, trans, ref, use_flirt=False, interp="trilinear", premat=None, postmat=None, mask=False, mask_thresh=0.5):
    """
    Transform an image using Regtricks for better performance relative to applywarp

    :param wsp: Workspace, used for logging only
    :param img: Image to transform
    :param trans: Transformation matrix or warp image
    :param ref: Reference image
    :param use_flirt: Use flirt to apply the transformation which must be a matrix
    :param interp: Interpolation method
    :param premat: If trans is a warp, this can be set to a pre-warp affine transformation matrix

    :return: Transformed Image object
    """
    if interp == "nn":
        order = 0
    else:
        order = 1

    have_warp = isinstance(trans, Image)
    ref_img_space = get_img_space(wsp, ref)
    if use_flirt and have_warp:
        raise ValueError("Cannot transform using Flirt when we have a warp")
    elif use_flirt:
        rt_trans = regtricks.Registration.from_flirt(trans, src=img.nibImage, ref=ref.nibImage)
    else:
        if have_warp:
            if premat is not None:
                if ref_img_space == "std":
                    # asl2std
                    rt_premat = regtricks.Registration.from_flirt(premat, src=img.nibImage, ref=wsp.reg.strucref.nibImage)
                    rt_warp = regtricks.NonLinearRegistration.from_fnirt(trans.nibImage, src=wsp.reg.strucref.nibImage, ref=ref.nibImage)
                    rt_trans = regtricks.chain(rt_premat, rt_warp)
                else:
                    # asl2asl
                    rt_premat = regtricks.Registration.from_flirt(premat, src=img.nibImage, ref=img.nibImage)
                    rt_warp = regtricks.NonLinearRegistration.from_fnirt(trans.nibImage, src=img.nibImage, ref=ref.nibImage)
                    rt_trans = regtricks.chain(rt_premat, rt_warp)
            elif postmat is not None:
                # std2asl
                rt_warp = regtricks.NonLinearRegistration.from_fnirt(trans.nibImage, src=img.nibImage, ref=wsp.reg.strucref.nibImage)
                rt_postmat = regtricks.Registration.from_flirt(postmat, src=wsp.reg.strucref.nibImage, ref=ref.nibImage)
                rt_trans = regtricks.chain(rt_warp, rt_postmat)
            else:
                rt_trans = regtricks.NonLinearRegistration.from_fnirt(trans.nibImage, src=img.nibImage, ref=ref.nibImage)
        else:
            if premat is not None:
                # asl2std
                rt_premat = regtricks.Registration.from_flirt(premat, src=img.nibImage, ref=wsp.reg.strucref.nibImage)
                rt_mat = regtricks.Registration.from_flirt(trans, src=wsp.reg.strucref.nibImage, ref=ref.nibImage)
                rt_trans = regtricks.chain(rt_premat, rt_mat)
            elif postmat is not None:
                # std2asl
                rt_mat = regtricks.Registration.from_flirt(trans, src=img.nibImage, ref=wsp.reg.strucref.nibImage)
                rt_postmat = regtricks.Registration.from_flirt(postmat, src=wsp.reg.strucref.nibImage, ref=ref.nibImage)
                rt_trans = regtricks.chain(rt_mat, rt_postmat)
            else:
                if trans.shape[0] != 4:
                    rt_trans = regtricks.MotionCorrection.from_mcflirt(trans, src=img.nibImage, ref=ref.nibImage)
                else:
                    rt_trans = regtricks.Registration.from_flirt(trans, src=img.nibImage, ref=ref.nibImage)

    wsp.log.write(" - Transforming image using regtricks\n")
    ret = Image(rt_trans.apply_to_image(img.nibImage, ref=ref.nibImage), order=order)
    if mask:
        # Binarise mask images
        ret = Image((ret.data > mask_thresh).astype(np.int32), header=ret.header)
    return ret

def reg_flirt(wsp, img, ref, initial_transform=None):
    """
    Register low resolution ASL or calibration data to a high resolution
    structural image using Flirt rigid-body registration

    The brain extracted structural image is used as the reference image. If
    this is not supplied, BET will be run on the whole head structural image.

    :param reg_img: Data to register, e.g. PWI or calibration image. Normally would be brain extracted
    :param struc_brain_img: Brain-extracted structural image

    Optional keyword arguments:

    :param inweight:
    :param init: Initial transform matrix
    :param schedule: FLIRT transform schedule file (default: xyztrans.sch")
    :param dof: FLIRT degrees of freedom

    :return Tuple of registered image, transform matrix
    """
    wsp.log.write(" - Registering image: %s using FLIRT\n" % img.name)

    # Step 1: 3D translation only
    flirt_opts = {
        "schedule" : os.path.join(os.environ["FSLDIR"], "etc", "flirtsch", "xyztrans.sch"),
        "init" : initial_transform,
        "inweight" : wsp.inweight,
        "log" : wsp.fsllog,
    }
    step1_trans = fsl.flirt(img, ref, omat=fsl.LOAD, **flirt_opts)["omat"]

    # Step 2: 6 DOF transformation with small search region
    flirt_opts.update({
        "schedule" : os.path.join(os.environ["FSLDIR"], "etc", "flirtsch", wsp.ifnone("flirtsch", "simple3D.sch")),
        "init" : step1_trans,
        "dof" : wsp.ifnone("dof", 6),
    })
    flirt_result = fsl.flirt(img, ref, out=fsl.LOAD, omat=fsl.LOAD, **flirt_opts)

    return flirt_result["out"], flirt_result["omat"]

def reg_bbr(wsp):
    """
    Perform BBR registration

    :param reg_img: Data to register, e.g. PWI or calibration image. Normally would be brain extracted
    :param struc_img: Structural image
    :param struc_brain_img: Brain-extracted structural image

    Optional keyword arguments:

    :param inweight:
    :param init: Initial transform matrix

    Optional keyword arguments for fieldmap distortion correction:

    :param fmap: Fieldmap image
    :param fmapmag: Fieldmap magnitude image
    :param fmapmagbrain: Fieldmap magnitude image - brain extracted
    :param pedir: Phase encoding direction (x, -x, y, -y, z, -z)
    :param echospacing: Echo spacing

    :return Tuple of registered image, transform matrix
    """
    wsp.log.write(" - BBR registration using epi_reg...")
    # Windows can't run epi_reg as it's a batch script. Use our experimental python
    # implementation but use the standard epi_reg on other platforms until the python
    # version is better tested
    if sys.platform.startswith("win"):
        import oxasl.epi_reg as pyepi
        result = pyepi.epi_reg(wsp, wsp.reg.aslref)
    else:
        result = epi_reg(epi=wsp.reg.aslref, t1=wsp.structural.struc, t1brain=wsp.structural.brain, out=fsl.LOAD, wmseg=wsp.structural.wm_seg, init=wsp.reg.asl2struc, inweight=wsp.inweight, log=wsp.fsllog)
    return result["out%s" % defaultExt()], result["out"]

    #OUTPUT
    #echo "Saving FINAL output"
    #if [ -z $finalonly ]; then
	#cp $outdir/asl2struct.mat $outdir/asl2struct_init.mat # save the initial transformation matrix to allow chekcing if this part failed
    #fi
    #cp $tempdir/low2high_final.mat $outdir/asl2struct.mat #the transformation matrix from epi_reg - this overwrites the version from MAIN registration
    #convert_xfm -omat $outdir/struct2asl.mat -inverse $outdir/asl2struct.mat #often useful to have the inverse transform, so calcuate it
    #if [ ! -z $fmap ]; then
	#imcp $tempdir/low2high_final_warp $outdir/asl2struct_warp #the warp from epi_reg
    #fi
    #imcp $tempdir/low2high_final $outdir/asl2struct # save the transformed image to check on the registration
    #
    # # copy the edge image from epi_reg output as that is good for visualisation
    # imcp $wm_seg $outdir/wm_seg
    #imcp $tempdir/low2high_final_fast_wmedge $outdir/tissedge

def get_transform_params(mat):
    """
    Get motion parameters from a Flirt motion correction matrix

    This is done under the assumption that the matrix may contain
    rotation, translation and possibly minor scaling but no reflection,
    shear etc. So the output could be incorrect for some extreme
    correction matrices, but this probably indicates an error in the
    registration process. We wrap the whole thing in a try block so
    if anything goes horribly wrong it does not at least stop the
    pipeline running

    See http://en.wikipedia.org/wiki/Rotation_matrix for details
    of the rotation calculation.

    :return: Tuple of magnitude of translation, angle and rotation axis
    """
    if tuple(mat.shape) != (4, 4):
        raise ValueError("Not a 4x4 Flirt matrix")

    try:
        # Extract scales - last one is the magnitude of the translation
        scales = np.linalg.norm(mat[:3, :], axis=0)

        # Normalise unit vectors by scaling before extracting rotation
        mat[:, 0] /= scales[0]
        mat[:, 1] /= scales[1]
        mat[:, 2] /= scales[2]

        # Rotation axis
        rot_axis = np.array([
            mat[2, 1] - mat[1, 2],
            mat[0, 2] - mat[2, 0],
            mat[1, 0] - mat[0, 1],
        ], dtype=np.float32)

        # Rotation angle - note that we need to check the sign
        trace = np.trace(mat[:3, :3])
        costheta = (float(trace)-1) / 2
        sintheta = math.sqrt(1-costheta*costheta)
        theta = math.acos(costheta)
        test_element = rot_axis[1]*rot_axis[0]*(1-costheta) + rot_axis[2]*sintheta
        if np.abs(test_element - mat[1, 0]) > np.abs(test_element - mat[0, 1]):
            theta = -theta
        return scales[-1], math.degrees(theta), rot_axis
    except:
        warnings.warn("Error extracting motion parameters from transformation matrix - check registration/moco looks OK!")
        traceback.print_exc()
        return 1, 0, [0, 0, 1]

def main():
    """
    Entry point for command line tool
    """
    try:
        parser = AslOptionParser(usage="asl_reg [options]", version=__version__)
        parser.add_category(Options())
        parser.add_category(struc.StructuralImageOptions())
        parser.add_category(GenericOptions())

        options, _ = parser.parse_args(sys.argv)
        wsp = Workspace(**vars(options))

        if not options.aslref:
            sys.stderr.write("Input file not specified\n")
            parser.print_help()
            sys.exit(1)

        reg_asl2struc(wsp, wsp.do_flirt, wsp.do_bbr)
        if wsp.output:
            wsp.reg.regto.save(wsp.output)
        if wsp.reg.asl2struc:
            with open(wsp.omat, "w") as transform_file:
                for row in wsp.reg.asl2struc:
                    transform_file.write(" ".join(["%f" % val for val in row]) + "\n")

    except ValueError as exc:
        sys.stderr.write("ERROR: " + str(exc) + "\n")
        sys.exit(1)

if __name__ == "__main__":
    main()
