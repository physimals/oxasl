#!/bin/env python
"""
ASL_REG: Registration for ASL data

Michael Chappell, IBME QuBIc & FMRIB Image Analysis Groups

Copyright (c) 2008-2018 University of Oxford
"""

import os
import sys

import numpy as np

import fsl.wrappers as fsl

from oxasl.workspace import Workspace
from . import calib, preproc, struc
from .options import AslOptionParser, GenericOptions, OptionCategory, IgnorableOptionGroup
from ._version import __version__
from .wrappers import epi_reg

def get_regfrom(wsp):
    """
    Set the 3D image to be used as the ASL registration target for structural->ASL registration
    
    Optional workspace attributes
    -----------------------------

     - ``regfrom`` : User-supplied registration reference image
     - ``asldata`` : Raw ASL data
     - ``calib``   : Calibration image

    Updated workspace attributes
    ----------------------------

     - ``regfrom``    : Registration reference image in ASL space
    """
    if wsp.isdone("get_regfrom"):
        return

    wsp.log.write("\nGetting image to use for ASL->structural registration)\n")
    calib.preproc_calib(wsp)
    preproc.preproc_asl(wsp)
    if wsp.regfrom is not None:
        wsp.log.write(" - Registration reference image supplied by user\n")
    elif wsp.asldata_mean_brain:
        wsp.log.write(" - Registration reference is mean raw ASL image (brain extracted)\n")
        wsp.regfrom = wsp.asldata_mean_brain
    elif wsp.calib_brain is not None:
        wsp.log.write(" - Registration reference is calibration image (brain extracted)\n")
        wsp.regfrom = wsp.calib_brain
    elif wsp.diffasl_mean_brain:
        wsp.log.write(" - Registration reference is mean subtracted ASL image (brain extracted)\n")
        wsp.regfrom = wsp.diffasl_mean_brain
    elif wsp.pwi_brain:
        wsp.log.write(" - Registration reference is perfusion weighted image (brain extracted)\n")
        wsp.regfrom = wsp.pwi_brain

    wsp.done("get_regfrom")

def reg_asl2struc(wsp):
    """
    Registration of ASL images to structural image
    
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
    if wsp.isdone("reg_asl2struc"):
        return

    get_regfrom(wsp)
    struc.preproc_struc(wsp)
    wsp.log.write("\nRegistering ASL data to structural data\n")
    if wsp.do_flirt:
        wsp.regto, wsp.asl2struc = reg_flirt(wsp)
    if wsp.do_bbr:
        wsp.regto, wsp.asl2struc = reg_bbr(wsp)
    
    wsp.struc2asl = np.linalg.inv(wsp.asl2struc)

    wsp.log.write(" - ASL->Structural transform\n")
    wsp.log.write(str(wsp.asl2struc) + "\n")
    wsp.log.write(" - Structural->ASL transform\n")
    wsp.log.write(str(wsp.struc2asl) + "\n")

    wsp.done("reg_asl2struc")

def reg_struc2std(wsp, flirt=True):
    """
    Determine structural -> standard space registration
    
    Required workspace attributes
    -----------------------------

     - ``struc``              : Structural image

    Updated workspace attributes
    ----------------------------

     - ``std2struc``    : MNI->structural transformation - either warp image or FLIRT matrix
     - ``struc2std``    : Structural->MNI transformation matrix - either warp image or FLIRT matrix
    """
    if wsp.isdone("reg_struc2std"):
        return

    if flirt:
        wsp.log.write(" - Registering structural image to standard space using FLIRT\n")
        flirt_result = fsl.flirt(wsp.struc, os.path.join(os.environ["FSLDIR"], "data/standard/MNI152_T1_2mm_brain"), omat=fsl.LOAD)
        wsp.struc2std = flirt_result["omat"]
        wsp.std2struc = np.linalg.inv(wsp.struc2std)
    else:
        # FIXME expecting flirt already run?
        wsp.log.write(" - Registering structural image to standard space using FNIRT\n")
        fnirt_result = fsl.fnirt(wsp.struc, aff=wsp.struc2std, config="T1_2_MNI152_2mm.cnf", cout=fsl.LOAD)
        wsp.struc2std = fnirt_result["cout"]
    
        # Calculate the inverse warp using INVWARP
        invwarp_result = fsl.invwarp(wsp.struc, wsp.struc2std_warp, out=fsl.LOAD)
        wsp.std2struc = invwarp_result["out"]

    wsp.done("reg_struc2std")

def struc2asl(wsp, img):
    """
    Convert an image from structural to ASL space

    :param img: Image object in structural space
    :return: Transformed Image object
    """
    get_regfrom(wsp)
    reg_asl2struc(wsp)
    return fsl.applyxfm(img, wsp.regfrom, wsp.struc2asl, out=fsl.LOAD, interp="trilinear", log=wsp.fsllog)["out"]

def reg_flirt(wsp):
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
    # This is done to avoid the contrast enhanced rim resulting from low intensity ref image
    # FIXME from ENABLE!
    #ref_data = ref_img.data
    #threshold = np.percentile(ref_data[ref_data != 0], 10.0)
    #ref_data[ref_data < threshold] = 0
    #ref_img = Image(ref_data, name=ref_img.name + "_thr", header=ref_img.header)
    
    wsp.log.write(" - Registering image: %s using FLIRT\n" % wsp.regfrom.name)
    
    # Step 1: 3D translation only
    flirt_opts = {
        "schedule" : os.path.join(os.environ["FSLDIR"], "etc", "flirtsch", "xyztrans.sch"),
        "init" : wsp.asl2struc,
        "inweight" : wsp.inweight,
        "log" : wsp.fsllog,
    }
    step1_trans = fsl.flirt(wsp.regfrom, wsp.struc_brain, omat=fsl.LOAD, **flirt_opts)["omat"]

    # Step 2: 6 DOF transformation with small search region
    flirt_opts.update({
        "schedule" : os.path.join(os.environ["FSLDIR"], "etc", "flirtsch", wsp.ifnone("flirtsch", "simple3D.sch")),
        "init" : step1_trans,
        "dof" : wsp.ifnone("dof", 6),
    })
    flirt_result = fsl.flirt(wsp.regfrom, wsp.struc_brain, out=fsl.LOAD, omat=fsl.LOAD, **flirt_opts)

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
    result = epi_reg(epi=wsp.regfrom, t1=wsp.struc, t1brain=wsp.struc_brain, out=fsl.LOAD, wmseg=wsp.wm_seg_struc, init=wsp.asl2struc, inweight=wsp.inweight)
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
        group.add_option("--omat", help="Output file for transform matrix", default=None)
        group.add_option("--bbr", dest="do_bbr", help="Include BBR registration step using EPI_REG", action="store_true", default=False)
        group.add_option("--flirt", dest="do_flirt", help="Include rigid-body registration step using FLIRT", action="store_true", default=True)
        group.add_option("--flirtsch", help="user-specified FLIRT schedule for registration")
        groups.append(group)
        
        #group = IgnorableOptionGroup(parser, "Extra BBR registration refinement", ignore=self.ignore)
        #group.add_option("-c", dest="cfile", help="ASL control/calibration image for initial registration - brain extracted")
        #group.add_option("--wm_seg", dest="wm_seg", help="tissue segmenation image for bbr (in structural image space)")
        #groups.append(group)

        #group = IgnorableOptionGroup(parser, "Distortion correction using fieldmap (see epi_reg)", ignore=self.ignore)
        #g.add_option("--fmap", dest="fmap", help="fieldmap image (in rad/s)")
        #g.add_option("--fmapmag", dest="fmapmag", help="fieldmap magnitude image - wholehead extracted")
        #g.add_option("--fmapmagbrain", dest="fmapmagbrain", help="fieldmap magnitude image - brain extracted")
        #g.add_option("--wmseg", dest="wmseg", help="white matter segmentation of T1 image")
        #g.add_option("--echospacing", dest="echospacing", help="Effective EPI echo spacing (sometimes called dwell time) - in seconds", type="float")
        #g.add_option("--pedir", dest="pedir", help="phase encoding direction, dir = x/y/z/-x/-y/-z")
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

        reg_asl2struc(wsp)
        if wsp.output:
            wsp.regto.save(wsp.output)
        if wsp.asl2struc:
            with open(wsp.omat, "w") as transform_file:
                for row in wsp.asl2struc:
                    transform_file.write(" ".join(["%f" % val for val in row]) + "\n")

    except ValueError as exc:
        sys.stderr.write("ERROR: " + str(exc) + "\n")
        sys.exit(1)

if __name__ == "__main__":
    main()
