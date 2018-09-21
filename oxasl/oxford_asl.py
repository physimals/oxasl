#!/bin/env python
"""
OXFORD_ASL: Converts ASL images into perfusion maps
===================================================

Michael Chappell, FMRIB Image Analysis Group & IBME

Copyright (c) 2008-2013 Univerisity of Oxford

Overview
--------

This module (and the associated command line tool) provide a 'one stop shop'
for analysing ASL data. They call the various other modules to preprocess,
correct, register, model and calibrate the data, outputting the perfusion
data, and also other estimated parameters.

The processing sequence is approximately:

1. Preprocess ASL data. This involves generation of mean and subtracted images
   and brain extracted versions. Not all of these images will be used later, but
   some might be.

2. Preprocess calibration data. This involves internal motion correction (for
   multi-volume calibration images), and the generation of mean and brain extracted
   versions of the calibration data.

3. Preprocess structural data. This involves generating brain extracted version of
   the structural data, and extracting any existing structural->standard space transformations
   (e.g. if an FSL_ANAT directory is the source of structural information)

4. If requested, the motion correction transformations for the ASL data are determined. 
   This is always done so that the middle volume of the
   ASL data is unchanged, although this volume is not always used as the reference.
   If applicable, the calibration->ASL registration is also determined.

5. If requested, the fieldmap distortion correction is determined from the fieldmap
   images

6. All available corrections are applied simultaneously to minimise interpolation

7. If required, a mask is generated from the corrected data

8. Model fitting is carried out to determine native space perfusion and other parameter maps

9. The ASL->structural registration is re-done using the perfusion image as the ASL reference.

10. If applicable, calibration is applied to the raw perfusion images

11. Raw and calibrated output images are generated in all output spaces (native, structural standard)
    
"""
from __future__ import print_function

import sys
import os
import traceback
import collections

import numpy as np

from fsl.data.image import Image

from oxasl import Workspace, AslImage, __version__, image, calib, struc, basil, preproc, mask, corrections, reg
from oxasl.options import AslOptionParser, GenericOptions, OptionCategory, IgnorableOptionGroup

class OxfordAslOptions(OptionCategory):
    """
    OptionCategory which contains options for preprocessing ASL data
    """

    def __init__(self, **kwargs):
        OptionCategory.__init__(self, "oxford_asl", **kwargs)

    def groups(self, parser):
        ret = []
        g = IgnorableOptionGroup(parser, "Main Options")
        g.add_option("--wp", dest="wp", help="Analysis which conforms to the 'white papers' (Alsop et al 2014)", action="store_true", default=False)
        g.add_option("--mc", dest="mc", help="Motion correct data", action="store_true", default=False)
        g.add_option("--fixbat", dest="inferbat", help="Fix bolus arrival time", action="store_false", default=True)
        ret.append(g)
        g = IgnorableOptionGroup(parser, "Acquisition/Data specific")
        g.add_option("--bat", dest="bat", help="Bolus arrival time (default=0.7 (pASL), 1.3 (cASL)", type=float)
        g.add_option("--t1", dest="t1", help="Tissue T1 value", type=float, default=1.3)
        g.add_option("--t1b", dest="t1b", help="Blood T1 value", type=float, default=1.65)
        ret.append(g)
        return ret

def main():
    """
    Entry point for oxasl command line tool
    """
    debug = True
    try:
        parser = AslOptionParser(usage="oxford_asl -i <asl_image> [options]", version=__version__)
        parser.add_category(image.AslImageOptions())
        parser.add_category(struc.StructuralImageOptions())
        parser.add_category(OxfordAslOptions())
        parser.add_category(calib.CalibOptions(ignore=["perf", "tis"]))
        parser.add_category(corrections.DistcorrOptions())
        parser.add_category(GenericOptions())

        options, _ = parser.parse_args(sys.argv)
        if not options.output:
            options.output = "oxasl"

        if not options.asldata:
            sys.stderr.write("Input file not specified\n")
            parser.print_help()
            sys.exit(1)

        # Some oxford_asl command-line specific defaults
        options.spatial = True
        if options.calib is not None and options.calib_method is None:
            if options.struc is not None:
                options.calib_method = "refregion"
            else:
                options.calib_method = "voxelwise"

        wsp = Workspace(savedir=options.output, **vars(options))
        wsp.asldata = AslImage(options.asldata, **parser.filter(vars(wsp), "image"))

        oxasl(wsp)
        report_build_dir = None
        if wsp.debug:
            report_build_dir = os.path.join(options.output, "report_build")
        wsp.report.generate_html(os.path.join(options.output, "report"), report_build_dir)

    except Exception as e:
        sys.stderr.write("ERROR: " + str(e) + "\n")
        if debug:
            traceback.print_exc()
        sys.exit(1)

def oxasl(wsp):
    """
    Run standard processing on ASL data

    This method requires wsp to be a Workspace containing certain standard information.
    As a minimum, the attribute ``asldata`` must contain an AslImage object.
    """
 
    wsp.log.write("OXFORD_ASL - running\n")
    wsp.log.write("Version: %s\n" % __version__)
    wsp.log.write("Input ASL data: %s\n" % wsp.asldata.name)
    wsp.asldata.summary(wsp.log)

    # Always output in native space
    wsp.output_spaces = ["native", ]
    wsp.do_flirt, wsp.do_bbr = True, False # FIXME

    preproc.preproc_asl(wsp)
    calib.preproc_calib(wsp)
    struc.preproc_struc(wsp)

    if wsp.mc: 
        corrections.get_motion_correction(wsp)

    if wsp.fmap:
        corrections.get_fieldmap_correction(wsp)

    corrections.apply_corrections(wsp)
    mask.generate_mask(wsp)

    basil.basil(wsp, output_wsp=wsp.sub("basil"))

    # Re-do ASL->structural registration using perfusion image
    if wsp.struc is not None:
        wsp.do_flirt, wsp.do_bbr = False, True # FIXME
        new_regfrom = wsp.basil.main.finalstep.mean_ftiss.data
        new_regfrom[new_regfrom < 0] = 0
        wsp.regfrom = Image(new_regfrom, header=wsp.regfrom.header)
        wsp.asl2struc_initial = wsp.asl2struc
        wsp.struc2asl_initial = wsp.struc2asl
        wsp.done("reg_asl2struc", status=False)
        reg.reg_asl2struc(wsp)

    if wsp.pvcorr:
        if wsp.mask_src == "struc":
            # Re-calculate mask as registration has changed
            wsp.mask = None
            mask.generate_mask(wsp)

        # Generate PVM and PWM maps for Basil
        struc.segment(wsp)
        wsp.wm_pv_asl = reg.struc2asl(wsp, wsp.wm_pv_struc, use_applywarp=True)
        wsp.gm_pv_asl = reg.struc2asl(wsp, wsp.gm_pv_struc, use_applywarp=True)
        wsp.basil_options = {"pwm" : wsp.wm_pv_asl, "pgm" : wsp.gm_pv_asl}
        wsp.initmvn = None
        basil.basil(wsp, output_wsp=wsp.sub("basil_pvcorr"), prefit=False)

    do_output(wsp)
    if not wsp.debug:
        cleanup(wsp)

    wsp.log.write("\nOutput is %s\n" % wsp.savedir)
    wsp.log.write("OXFORD_ASL - done\n")

def do_output(wsp):
    for space in wsp.output_spaces:
        output_wsp = wsp.sub(space)
        # Make negative values = 0
        img = wsp.basil.main.finalstep.mean_ftiss
        img.data[img.data < 0] = 0
        output_wsp.perfusion = img
        if wsp.calib is not None:
            output_wsp.perfusion_calib = calib.calibrate(wsp, output_wsp.perfusion)

def _deletable(value):
    return isinstance(value, (Image, np.ndarray)) or isinstance(value, collections.Sequence) or isinstance(value, collections.Mapping) 

def cleanup(wsp):
    output_spaces = wsp.output_spaces
    if output_spaces is None: output_spaces = []
    for attr in dir(wsp):
        value = getattr(wsp, attr)
        if not attr.startswith("__") and attr not in output_spaces:
            if isinstance(value, Workspace):
                print("Deleting workspace %s" % attr)
                cleanup(value)
            elif _deletable(value):
                print("Deleting %s" % attr)
                delattr(wsp, attr)
            else:
                print("Keeping %s (%s)" % (attr, value))
