#!/bin/env python
"""
OXASL: ASL processing pipeline
==============================

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

# Quick-and-dirty plugin imports - to be replaced with entry points and a defined plugin api
# at some point
try:
    import oxasl_ve
except ImportError:
    oxasl_ve = None

try:
    import oxasl_mp
except ImportError:
    oxasl_mp = None

try:
    import oxasl_deblur
except ImportError:
    oxasl_deblur = None

try:
    import oxasl_enable
except ImportError:
    oxasl_enable = None

try:
    import oxasl_surfpvc
except ImportError:
    oxasl_surfpvc = None

try:
    import oxasl_multite
except ImportError:
    oxasl_multite = None

from oxasl import *
from oxasl.options import AslOptionParser, GenericOptions, OptionCategory, OptionGroup
from oxasl.reporting import LightboxImage

def add_options(parser):
    g = OptionGroup(parser, "General Pipeline Options")
    g.add_option("--wp", help="Analysis which conforms to the 'white papers' (Alsop et al 2014)", action="store_true", default=False)
    g.add_option("--mc", help="Motion correct data", action="store_true", default=False)
    if oxasl_enable:
        g.add_option("--use-enable", help="Use ENABLE preprocessing step", action="store_true", default=False)
    if oxasl_deblur:
        g.add_option("--deblur", help="Perform deblurring of ASL and calibration data", action="store_true", default=False)
    parser.add_option_group(g)

    g = OptionGroup(parser, "Model fitting options")
    g.add_option("--fixbat", dest="inferbat", help="Fix bolus arrival time", action="store_false", default=True)
    g.add_option("--batsd", help="Bolus arrival time standard deviation (s) - default 1.0 for multi-PLD, 0.1 otherwise", type=float)
    g.add_option("--fixbolus", "--fixtau", dest="infertau", help="Fix bolus duration", action="store_false")
    g.add_option("--art-off", "--artoff", dest="inferart", help="Do not infer arterial component", action="store_false", default=True)
    g.add_option("--spatial-off", "--spatialoff", dest="spatial", help="Do not include adaptive spatial smoothing on CBF", action="store_false", default=True)
    g.add_option("--infertexch", help="Infer exchange time (multi-TE data only)", action="store_true", default=False)
    g.add_option("--infert1", help="Infer T1 value", action="store_true", default=False)
    g.add_option("--infert2", help="Infer T2 value (multi-TE data only)", action="store_true", default=False)
    g.add_option("--t1im", help="Voxelwise T1 tissue estimates", type="image")
    g.add_option("--batim", "--attim", help="Voxelwise BAT (ATT) estimates in seconds", type="image")
    g.add_option("--basil-mask", help="Masking policy to use for Basil model fitting. Does not affect analysis mask used in rest of pipeline. 'dilate' means dilate the default analysis mask. 'none' means use no masking",
                 type="choice", choices=["default", "dilated", "none"])
    g.add_option("--basil-options", "--fit-options", help="File containing additional options for model fitting step", type="optfile", default=None)
    parser.add_option_group(g)
    
    g = OptionGroup(parser, "Physiological parameters (all have default values from literature)")
    g.add_option("--bat", help="Estimated bolus arrival time (s) - default=0.7 (pASL), 1.3 (cASL)", type=float)
    g.add_option("--t1", "--t1t", help="Tissue T1 (s)", type=float, default=1.3)
    g.add_option("--t2", "--t2t", help="Tissue T2 (ms)", type=float, default=50)
    g.add_option("--t2s", help="Tissue T2* (ms)", type=float, default=20)
    g.add_option("--t1b", help="Blood T1 (s)", type=float, default=1.65)
    g.add_option("--t2b", help="Blood T2 (ms) - Lu et al. 2012 MRM 67:42-49, 3T during normoxia", type=float, default=150)
    g.add_option("--t2sb", help="Blood T2* (ms) - Petersen 2006 MRM 55(2):219-232", type=float, default=50)
    parser.add_option_group(g)

def main():
    """
    Entry point for oxasl command line tool
    """
    debug = True
    wsp = None
    try:
        parser = AslOptionParser(usage="oxasl -i <asl_image> [options]", version=__version__)
        parser.add_category(image.Options())
        parser.add_category(struc.Options())
        add_options(parser)
        m0.add_options(parser)
        parser.add_category(reg.Options())
        distcorr.add_options(parser)
        parser.add_category(senscorr.Options())
        parser.add_category(pvc.Options())
        parser.add_category(corrections.Options())
        if oxasl_ve:
            parser.add_category(oxasl_ve.VeaslOptions())
        if oxasl_mp:
            parser.add_category(oxasl_mp.MultiphaseOptions())
        if oxasl_enable:
            parser.add_category(oxasl_enable.EnableOptions())
        if oxasl_multite:
            parser.add_category(oxasl_multite.MultiTEOptions())
        if oxasl_deblur:
            parser.add_category(oxasl_deblur.Options())
        parser.add_category(region_analysis.Options())
        parser.add_category(GenericOptions())
        parser.add_category(output.Options())

        options, _ = parser.parse_args()
        debug = options.debug

        # Some oxasl command-line specific defaults
        if (options.calib is not None or options.calib_first_vol) and options.calib_method is None:
            if options.struc is not None:
                options.calib_method = "refregion"
            else:
                options.calib_method = "voxelwise"

        if options.debug:
            options.save_all = True
        options.output_native = True
        options.output_struc = True
        options.save_mask = True

        if options.asldata is None:
            raise RuntimeError("Input ASL file not specified\n")

        if os.path.exists(options.output) and not options.overwrite:
            raise RuntimeError("Output directory exists - use --overwrite to overwrite it")

        wsp = Workspace(savedir=options.output, auto_asldata=True, **vars(options))
        oxasl(wsp)

    except Exception as e:
        sys.stderr.write("ERROR: " + str(e) + "\n")
        if debug:
            traceback.print_exc()
        if wsp is not None:
            wsp.log.write("ERROR: " + str(e) + "\n")
            wsp.log.write(traceback.format_exc() + "\n")
        sys.exit(1)

def oxasl(wsp):
    """
    Main oxasl pipeline script
    """
    wsp.log.write("OXASL version: %s\n" % __version__)
    for plugin in (oxasl_ve, oxasl_mp, oxasl_deblur, oxasl_enable, oxasl_surfpvc, oxasl_multite):
        if plugin is not None:
            wsp.log.write(" - Found plugin: %s (version %s)\n" % (plugin.__name__, getattr(plugin, "__version__", "unknown")))

    wsp.log.write("\nInput ASL data: %s\n" % wsp.asldata.name)
    wsp.asldata.summary(wsp.log)
    report_asl(wsp)

    # Preprocessing and preparation
    preproc.run(wsp)
    struc.run(wsp)
    moco.run(wsp)
    corrections.run(wsp)
    reg.run(wsp)
    distcorr.run(wsp)
    senscorr.run(wsp)
    corrections.run(wsp)
    mask.run(wsp)

    # Quantification
    filter.run(wsp)
    prequantify.run(wsp)
    quantify.run(wsp)
    pvc.run(wsp)

    # Calibration and output
    m0.run(wsp)
    output.run(wsp)

    # Post processing and reporting
    region_analysis.run(wsp.output)
    if wsp.pvcorr:
        region_analysis.run(wsp.output_pvcorr)

    if wsp.save_report:
        reporting.run(wsp)

    _cleanup(wsp)
    wsp.log.write("\nOutput is %s\n" % wsp.savedir)
    wsp.log.write("OXASL - done\n")

def report_asl(wsp):
    """
    Generate a report page about the input ASL data
    """
    page = wsp.report.page("asl")
    page.heading("ASL input data")
    md_table = [(key, value) for key, value in wsp.asldata.metadata_summary().items()]
    page.table(md_table)
    if wsp.pwi is not None:
        img = wsp.pwi
        img_type = "Perfusion-weighted image"
    else:
        img = wsp.asldata.mean()
        img_type = "Mean ASL data"

    page.heading(img_type, level=1)
    page.image("asldata", LightboxImage(img))

def _cleanup(wsp):
    """
    Remove items from the workspace that are not being output. The
    corresponding files will be deleted.
    """
    wsp.log.write("\nDoing cleanup\n")
    if not wsp.save_all:
        if not wsp.save_preproc:
            wsp.log.write(" - Removing preprocessed data\n")
        if not wsp.save_corrected:
            wsp.log.write(" - Removing corrected data\n")
            wsp.corrected = None
        if not wsp.save_corrections:
            wsp.log.write(" - Removing intermediate correction data\n")
            wsp.senscorr = None
            wsp.distcorr = None
            wsp.moco = None
            wsp.topup = None
            wsp.fieldmap = None
        if not wsp.save_reg:
            wsp.log.write(" - Removing registration data\n")
            wsp.reg = None
        if not wsp.save_basil:
            wsp.log.write(" - Removing model fitting data\n")
            wsp.basil = None
        if not wsp.save_struc:
            wsp.log.write(" - Removing structural segementation data\n")
            wsp.structural = None
        if not wsp.save_calib:
            wsp.log.write(" - Removing calibration data\n")
            wsp.calibration = None
        if not wsp.save_filter:
            wsp.filter = None
        if not wsp.save_input:
            wsp.input = None
        wsp.rois = None
