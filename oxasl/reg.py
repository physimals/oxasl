#!/bin/env python
"""
ASL_REG: Registration for ASL data

Michael Chappell, IBME QuBIc & FMRIB Image Analysis Groups

Copyright (c) 2008-2018 University of Oxford
"""
from __future__ import unicode_literals

import os
import sys
import math
import warnings
import traceback

import numpy as np

from fsl.data.image import Image
import fsl.wrappers as fsl

from oxasl import __version__, Workspace, struc, brain
from oxasl.options import AslOptionParser, GenericOptions, OptionCategory, IgnorableOptionGroup, load_matrix
from oxasl.wrappers import epi_reg
from oxasl.reporting import LightboxImage

def init(wsp):
    """
    Create registration sub-workspace if not already there
    """
    if wsp.reg is None:
        wsp.sub("reg")

def get_regfrom(wsp):
    """
    Set the 3D image to be used as the ASL registration target for structural->ASL registration

    Regfrom defines the 'native' space

    Optional workspace attributes
    -----------------------------

     - ``regfrom`` : User-supplied registration reference image
     - ``asldata`` : Raw ASL data
     - ``calib``   : Calibration image

    Updated workspace attributes
    ----------------------------

     - ``regfrom``    : Registration reference image in ASL space
    """
    init(wsp)
    if wsp.reg.regfrom is None:
        wsp.log.write("\nGetting the ASL image to use for registration)\n")
        if wsp.regfrom is not None:
            wsp.log.write(" - Registration reference image supplied by user\n")
            wsp.reg.regfrom = wsp.regfrom
        elif wsp.asldata.iaf in ("tc", "ct"):
            wsp.log.write(" - Registration reference is mean ASL signal (brain extracted)\n")
            wsp.reg.regfrom = brain.brain(wsp, wsp.asldata.mean(), thresh=0.2)
        elif wsp.calib is not None and wsp.calib.sameSpace(wsp.asldata):
            wsp.log.write(" - Registration reference is calibration image (brain extracted)\n")
            wsp.reg.regfrom = brain.brain(wsp, wsp.calib, thresh=0.2)
        else:
            wsp.log.write(" - Registration reference is mean ASL image (brain extracted)\n")
            wsp.reg.regfrom = brain.brain(wsp, wsp.asldata.mean(), thresh=0.2)

def get_motion_params(mat):
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

    :return: magnitude of translation, angle and rotation axis
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
        ], dtype=np.float)

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

def reg_asl2calib(wsp):
    """
    Register calibration image to ASL space

    Note that this might already have been done as part of motion correction
    """
    init(wsp)
    if wsp.calib is not None and wsp.reg.asl2calib is None:
        get_regfrom(wsp)
        wsp.log.write("Registering calibration image to ASL image\n")
        _, wsp.reg.asl2calib = reg_flirt(wsp, wsp.reg.regfrom, wsp.calib)
        wsp.reg.calib2asl = np.linalg.inv(wsp.reg.asl2calib)

def reg_asl2struc(wsp, flirt=True, bbr=False, name="initial"):
    """
    Registration of ASL images to structural image

    :param flirt: If provided, sets whether to use FLIRT registration
    :param bbr: If provided, sets whether to use BBR registration

    Required workspace attributes
    -----------------------------

     - ``regfrom``            : Registration reference image in ASL space
     - ``struc``              : Structural image

    Updated workspace attributes
    ----------------------------

     - ``asl2struc``    : ASL->structural transformation matrix
     - ``struc2asl``    : Structural->ASL transformation matrix
     - ``regto``        : ``regfrom`` image transformed to structural space
    """
    init(wsp)
    struc.init(wsp)
    if wsp.structural.struc is not None:
        get_regfrom(wsp)
        wsp.log.write("\nRegistering ASL data to structural data\n")
        if flirt:
            wsp.reg.regto, wsp.reg.asl2struc = reg_flirt(wsp, wsp.reg.regfrom, wsp.structural.brain, wsp.reg.asl2struc)
        if bbr:
            wsp.reg.regto, wsp.reg.asl2struc = reg_bbr(wsp)

        wsp.reg.struc2asl = np.linalg.inv(wsp.reg.asl2struc)

        wsp.log.write(" - ASL->Structural transform\n")
        wsp.log.write(str(wsp.reg.asl2struc) + "\n")
        wsp.log.write(" - Structural->ASL transform\n")
        wsp.log.write(str(wsp.reg.struc2asl) + "\n")

        page = wsp.report.page("asl2struc_%s" % name)
        page.heading("%s ASL -> Structural registration" % name.title(), level=0)
        page.heading("Transformation parameters", level=1)
        motion_params = get_motion_params(wsp.reg.asl2struc)
        page.table([
            ["Translation magnitude", "%.3g mm" % motion_params[0]],
            ["Rotation magnitude", "%.3g \N{DEGREE SIGN}" % motion_params[1]],
        ])
        page.heading("ASL->Structural transformation matrix", level=1)
        page.matrix(wsp.reg.asl2struc)
        page.heading("Structural->ASL transformation matrix", level=1)
        page.matrix(wsp.reg.struc2asl)

        if wsp.structural.gm_seg is not None:
            gm_asl = struc2asl(wsp, wsp.structural.gm_seg, interp="nn")
            page.heading("GM mask aligned with ASL data", level=1)
            page.image("gm_reg_%s" % name, LightboxImage(gm_asl, bgimage=wsp.reg.regfrom))
            wm_asl = struc2asl(wsp, wsp.structural.wm_seg, interp="nn")
            page.heading("WM mask aligned with ASL data", level=1)
            page.image("wm_reg_%s" % name, LightboxImage(wm_asl, bgimage=wsp.reg.regfrom))

def reg_struc2std(wsp, fnirt=False):
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
    init(wsp)

    if wsp.reg.std2struc is not None:
        return

    if wsp.fslanat:
        warp = os.path.join(wsp.fslanat, "T1_to_MNI_nonlin_coeff.nii.gz")
        mat = os.path.join(wsp.fslanat, "T1_to_MNI_lin.mat")
        if os.path.isfile(warp):
            wsp.log.write(" - Using structural->std nonlinear transformation from FSL_ANAT\n")
            wsp.reg.struc2std = Image(warp, loadData=False)
        elif os.path.isfile(mat):
            wsp.log.write(" - Using structural->std linear transformation from FSL_ANAT\n")
            wsp.reg.struc2std = load_matrix(mat)

    if wsp.reg.struc2std is None:
        struc.init(wsp)
        wsp.log.write(" - Registering structural image to standard space using FLIRT\n")
        flirt_result = fsl.flirt(wsp.structural.brain, os.path.join(os.environ["FSLDIR"], "data/standard/MNI152_T1_2mm_brain"), omat=fsl.LOAD)
        wsp.reg.struc2std = flirt_result["omat"]

        if fnirt:
            wsp.log.write(" - Registering structural image to standard space using FNIRT\n")
            fnirt_result = fsl.fnirt(wsp.structural.brain, aff=wsp.reg.struc2std, config="T1_2_MNI152_2mm.cnf", cout=fsl.LOAD)
            wsp.reg.struc2std = fnirt_result["cout"]

    if isinstance(wsp.reg.struc2std, Image):
        # Calculate the inverse warp using INVWARP
        invwarp_result = fsl.invwarp(wsp.reg.struc2std, wsp.structural.struc, out=fsl.LOAD)
        wsp.reg.std2struc = invwarp_result["out"]
    else:
        wsp.reg.std2struc = np.linalg.inv(wsp.reg.struc2std)

def std2struc(wsp, img, **kwargs):
    """
    Transform an image from standard space to structural space
    """
    return transform(wsp, img, wsp.reg.std2struc, wsp.structural.struc, **kwargs)

def struc2std(wsp, img, **kwargs):
    """
    Transform an image from structural space to standard space
    """
    ref = Image(os.path.join(os.environ["FSLDIR"], "data/standard/MNI152_T1_2mm_brain"))
    return transform(wsp, img, wsp.reg.struc2std, ref, **kwargs)

def struc2asl(wsp, img, **kwargs):
    """
    Convert an image from structural to ASL space

    :param img: Image object in structural space
    :return: Transformed Image object in ASL (native) space
    """
    init(wsp)
    return transform(wsp, img, wsp.reg.struc2asl, wsp.nativeref, **kwargs)

def asl2struc(wsp, img, **kwargs):
    """
    Convert an image from ASL to structural space

    Keyword arguments are passed to ``transform``

    :param img: Image object in native (ASL) space
    :return: Transformed Image object in structural space
    """
    init(wsp)
    return transform(wsp, img, wsp.reg.asl2struc, wsp.structural.struc, **kwargs)

def calib2asl(wsp, img, **kwargs):
    """
    Convert an image from calibration space to ASL space

    :param img: Image object in calibration space
    :return: Transformed Image object in ASL (native) space
    """
    init(wsp)
    return transform(wsp, img, wsp.reg.calib2asl, wsp.nativeref, **kwargs)

def asl2calib(wsp, img, **kwargs):
    """
    Convert an image from ASL to calibration space

    Keyword arguments are passed to ``transform``

    :param img: Image object in native (ASL) space
    :return: Transformed Image object in calibration space
    """
    init(wsp)
    return transform(wsp, img, wsp.reg.asl2calib, wsp.structural.struc, **kwargs)

def transform(wsp, img, trans, ref, use_flirt=False, interp="trilinear", paddingsize=1, premat=None, mask=False, mask_thresh=0.5):
    """
    Transform an image

    :param wsp: Workspace, used for logging only
    :param img: Image to transform
    :param trans: Transformation matrix or warp image
    :param ref: Reference image
    :param use_flirt: Use flirt to apply the transformation which must be a matrix
    :param interp: Interpolation method
    :param paddingsize: Padding size in pixels
    :param premat: If trans is a warp, this can be set to a pre-warp affine transformation matrix

    :return: Transformed Image object
    """
    if trans is None:
        raise ValueError("Transformation matrix not available - has registration been performed?")

    have_warp = isinstance(trans, Image)
    if use_flirt and have_warp:
        raise ValueError("Cannot transform using Flirt when we have a warp")
    elif use_flirt:
        if interp == "nn":
            interp = "nearestneighbour"
        ret = fsl.applyxfm(img, ref, trans, out=fsl.LOAD, interp=interp, paddingsize=paddingsize, log=wsp.fsllog)["out"]
    else:
        if have_warp:
            kwargs = {"warp" : trans, "premat" : premat, "rel" : True}
        elif premat:
            raise ValueError("Can't set a pre-transformation matrix unless using a warp")
        else:
            kwargs = {"premat" : trans}
        ret = fsl.applywarp(img, ref, out=fsl.LOAD, interp=interp, paddingsize=paddingsize, super=True, superlevel="a", log=wsp.fsllog, **kwargs)["out"]
    if mask:
        # Binarise mask images
        ret = Image((ret.data > mask_thresh).astype(np.int), header=ret.header)
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
    struc.segment(wsp)

    wsp.log.write("  - BBR registration using epi_reg\n")
    # Windows can't run epi_reg as it's a batch script. Use our experimental python
    # implementation but use the standard epi_reg on other platforms until the python
    # version is better tested
    if sys.platform.startswith("win"):
        import oxasl.epi_reg as pyepi
        result = pyepi.epi_reg(wsp, wsp.reg.regfrom)
    else:
        result = epi_reg(epi=wsp.reg.regfrom, t1=wsp.structural.struc, t1brain=wsp.structural.brain, out=fsl.LOAD, wmseg=wsp.structural.wm_seg, init=wsp.reg.asl2struc, inweight=wsp.inweight, log=wsp.fsllog)
    return result["out.nii.gz"], result["out"]

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

class RegOptions(OptionCategory):
    """
    OptionCategory which contains options for registration of ASL data to structural image
    """

    def __init__(self, **kwargs):
        OptionCategory.__init__(self, "reg", **kwargs)

    def groups(self, parser):
        groups = []

        group = IgnorableOptionGroup(parser, "Registration", ignore=self.ignore)
        group.add_option("--regfrom", help="Registration image (e.g. perfusion weighted image)", type="image")
        #group.add_option("--omat", help="Output file for transform matrix", default=None)
        #group.add_option("--bbr", dest="do_bbr", help="Include BBR registration step using EPI_REG", action="store_true", default=False)
        #group.add_option("--flirt", dest="do_flirt", help="Include rigid-body registration step using FLIRT", action="store_true", default=True)
        #group.add_option("--flirtsch", help="user-specified FLIRT schedule for registration")
        groups.append(group)

        #group = IgnorableOptionGroup(parser, "Extra BBR registration refinement", ignore=self.ignore)
        #group.add_option("-c", dest="cfile", help="ASL control/calibration image for initial registration - brain extracted")
        #group.add_option("--wm_seg", dest="wm_seg", help="tissue segmenation image for bbr (in structural image space)")
        #groups.append(group)

        #group = IgnorableOptionGroup(parser, "Distortion correction using fieldmap (see epi_reg)", ignore=self.ignore)
        #g.add_option("--nofmapreg", dest="nofmapreg", help="do not perform registration of fmap to T1 (use if fmap already registered)", action="store_true", default=False)
        #groups.append(group)

        #group = IgnorableOptionGroup(parser, "Deprecated", ignore=self.ignore)
        #g.add_option("-r", dest="lowstruc", help="extra low resolution structural image - brain extracted")
        #g.add_option("--inweight", dest="inweight", help="specify weights for input image - same functionality as the flirt -inweight option", type="float")
        #groups.append(group)

        return groups

def main():
    """
    Entry point for command line tool
    """
    try:
        parser = AslOptionParser(usage="asl_reg [options]", version=__version__)
        parser.add_category(RegOptions())
        parser.add_category(struc.StructuralImageOptions())
        parser.add_category(GenericOptions())

        options, _ = parser.parse_args(sys.argv)
        wsp = Workspace(**vars(options))

        if not options.regfrom:
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
