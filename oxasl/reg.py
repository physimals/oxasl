#!/bin/env python
"""
ASL_REG: Registration for ASL data

Michael Chappell, IBME QuBIc & FMRIB Image Analysis Groups

Copyright (c) 2008-2018 University of Oxford
"""

import os
import sys

from fsl.wrappers import LOAD, flirt

import oxasl.struc as struc

from oxasl.workspace import AslWorkspace
from .options import AslOptionParser
from ._version import __version__

def reg_flirt(reg_img, struc_brain_img, log=sys.stdout, **kwargs):
    """ 
    Register low resolution ASL or calibration data to a high resolution
    structural image using Flirt rigid-body registration

    The brain extracted structural image is used as the reference image. If
    this is not supplied, BET will be run on the whole head structural image.

    :param reg_img: Data to register, e.g. PWI or calibration image. Normally would be brain extracted
    :param struc_brain_img: Brain-extracted structural image
    :return Tuple of registered image, transform matrix
    """
    log.write("Registration to structural data...\n")
        
    # This is done to avoid the contrast enhanced rim resulting from low intensity ref image
    # FIXME from ENABLE!
    #ref_data = ref_img.data
    #threshold = np.percentile(ref_data[ref_data != 0], 10.0)
    #ref_data[ref_data < threshold] = 0
    #ref_img = Image(ref_data, name=ref_img.name + "_thr", header=ref_img.header)
    
    log.write(" - Registering image: %s\n" % reg_img.name)

    # Step 1: 3D translation only
    flirt_opts = {
        "schedule" : os.path.join(os.environ["FSLDIR"], "etc", "flirtsch", "xyztrans.sch"),
        "init" : kwargs.get("initmat", None),
        "inweight" : kwargs.get("inweight", None),
        "log" : kwargs.get("fsllog", {"cmd" : log}),
    }
    step1_trans = flirt(reg_img, struc_brain_img, omat=LOAD, **flirt_opts)["omat"]

    # Step 2: 6 DOF transformation with small search region
    flirt_opts.update({
        "schedule" : os.path.join(os.environ["FSLDIR"], "etc", "flirtsch", kwargs.get("flirtsch", "simple3D.sch")),
        "init" : step1_trans,
        "dof" : kwargs.get("dof", 6),
    })
    flirt_result = flirt(reg_img, struc_brain_img, out=LOAD, omat=LOAD, **flirt_opts)

    log.write("DONE registration\n\n")
    return flirt_result["out"], flirt_result["omat"]

def reg_bbr(reg_img, struc_img, struc_brain_img, wm_seg, log=sys.stdout, **kwargs):
    """
    Perform BBR registration
    """
    log.write("BBR registration using epi_reg\n")

    # brain extract the perfsion image - using supplied mask or mask derived from the strctural BET
    #mask = kwargs.get("mask", None)
    #if mask is None:
    #    log.write(" - Generating brain mask\n")
    #    #convert_xfm -omat $tempdir/high2low.mat -inverse $tempdir/low2high.mat
    #    raise NotImplementedError("Mask generation disabled")
        #struc_brain_mask = fslmaths(struc_brain_img).thr(0).bin().run()
	    #flirt(struc_brain_mask, $infile -applyxfm -init $tempdir/high2low.mat -out $tempdir/mask -interp trilinear
	    #fslmaths $tempdir/mask -thr 0.25 -bin -fillh $tempdir/mask
	    #fslcpgeom $infile $tempdir/mask
	    #mask=$tempdir/mask
        
    # Refinement of the registration using perfusion and the white matter segmentation
    # using epi_reg to get BBR (and allow for fielmap correction in future)
    epi_reg_opts = {
        "inweight" : kwargs.get("inweight", None),
        "init" : kwargs.get("initmat", None),
    }
    #if fmap:
    if False:
    	# With fieldmap
        # FIXME
	    #fmapregstr=""
        #if [ ! -z $nofmapreg ]; then
	    #    fmapregstr="--nofmapreg"
        epi_reg_opts.update({
            "fmap" : kwargs.get("fmap", None),
            "fmapmag" : kwargs.get("fmapmag", None),
            "fmapmagbrain" : kwargs.get("fmapmagbrain", None),
            "pedir" : kwargs.get("pedir", None),
            "echospacing" : kwargs.get("echospacing", None),
        })
    
    trans = epi_reg(reg_img, t1=struc_img, t1brain=struc_brain_img, wmseg=wm_seg, **epi_reg_opts)
    log.write("BBR end\n")
    return trans

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

def main():
    """
    Entry point for command line tool
    """
    try:
        parser = AslOptionParser(usage="asl_reg -i <aslimg> [options]", version=__version__)
        parser.add_option("-i", dest="regfrom", help="Registration image (e.g. perfusion weighted image)", type="image")
        parser.add_option("-o", "--output", dest="output", help="Output file for registered image", default=None)
        parser.add_option("-m", "--mask", dest="mask", help="brain mask for brain extraction of the input image")
        parser.add_option("--omat", dest="omat", help="Output file for transform matrix", default=None)
        parser.add_option("--debug", dest="debug", help="Debug mode", action="store_true", default=False)
        parser.add_option("--bbr", help="Include BBR registration step using EPI_REG", action="store_true", default=False)
        parser.add_option("--flirt", help="Include rigid-body registration step using FLIRT", action="store_true", default=True)
        parser.add_option("--flirtsch", dest="flirtsch", help="user-specified FLIRT schedule for registration")
        
        parser.add_category(struc.StructuralImageOptions())

        options, _ = parser.parse_args(sys.argv)
        wsp = AslWorkspace(**vars(options))
        
        if not options.regfrom:
            sys.stderr.write("Input file not specified\n")
            parser.print_help()
            sys.exit(1)

        wsp.reg_asl2struc(do_bbr=options.bbr, do_flirt=options.flirt)
        if wsp.output:
            wsp.regto.save(wsp.output)
        if wsp.omat:
            with open(wsp.omat, "w") as transform_file:
                for row in wsp.asl2struc:
                    transform_file.write(" ".join(["%f" % val for val in row]) + "\n")

    except ValueError as exc:
        sys.stderr.write("ERROR: " + str(exc) + "\n")
        sys.exit(1)

    #g = OptionGroup(p, "Extra 'final' registration refinement")
    #g.add_option("-c", dest="cfile", help="ASL control/calibration image for initial registration - brain extracted")
    #g.add_option("--wm_seg", dest="wm_seg", help="tissue segmenation image for bbr (in structural image space)")
    #p.add_option_group(g)

    #g = OptionGroup(p, "Distortion correction using fieldmap (see epi_reg)")
    #g.add_option("--fmap", dest="fmap", help="fieldmap image (in rad/s)")
    #g.add_option("--fmapmag", dest="fmapmag", help="fieldmap magnitude image - wholehead extracted")
    #g.add_option("--fmapmagbrain", dest="fmapmagbrain", help="fieldmap magnitude image - brain extracted")
    #g.add_option("--wmseg", dest="wmseg", help="white matter segmentation of T1 image")
    #g.add_option("--echospacing", dest="echospacing", help="Effective EPI echo spacing (sometimes called dwell time) - in seconds", type="float")
    #g.add_option("--pedir", dest="pedir", help="phase encoding direction, dir = x/y/z/-x/-y/-z")
    #g.add_option("--nofmapreg", dest="nofmapreg", help="do not perform registration of fmap to T1 (use if fmap already registered)", action="store_true", default=False)
    #p.add_option_group(g)

    #g = OptionGroup(p, "Deprecated")
    #g.add_option("-r", dest="lowstruc", help="extra low resolution structural image - brain extracted")
    #g.add_option("--inweight", dest="inweight", help="specify weights for input image - same functionality as the flirt -inweight option", type="float")
    #p.add_option_group(g)

if __name__ == "__main__":
    main()
