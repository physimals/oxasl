#!/bin/env python
"""
OXASL: Converts ASL images into perfusion maps
==============================================

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

import numpy as np

from fsl.data.image import Image

try:
    import oxasl_ve
except ImportError:
    oxasl_ve = None

try:
    import oxasl_deblur
except ImportError:
    oxasl_deblur = None

try:
    import oxasl_enable
except ImportError:
    oxasl_enable = None

from oxasl import Workspace, __version__, image, calib, struc, basil, mask, corrections, reg
from oxasl.options import AslOptionParser, GenericOptions, OptionCategory, IgnorableOptionGroup
from oxasl.reporting import LightboxImage

class OxfordAslOptions(OptionCategory):
    """
    OptionCategory which contains options for preprocessing ASL data
    """

    def __init__(self, **kwargs):
        OptionCategory.__init__(self, "oxasl", **kwargs)

    def groups(self, parser):
        ret = []
        g = IgnorableOptionGroup(parser, "Main Options")
        g.add_option("--wp", help="Analysis which conforms to the 'white papers' (Alsop et al 2014)", action="store_true", default=False)
        g.add_option("--mc", help="Motion correct data", action="store_true", default=False)
        g.add_option("--fixbat", dest="inferbat", help="Fix bolus arrival time", action="store_false", default=True)
        g.add_option("--fixbolus", dest="infertau", help="Fix bolus duration", action="store_false")
        g.add_option("--artoff", dest="inferart", help="Do not infer arterial component", action="store_false", default=True)
        g.add_option("--spatial-off", dest="spatial", help="Do not include adaptive spatial smoothing on CBF", action="store_false", default=True)
        if oxasl_enable:
            g.add_option("--use-enable", help="Use ENABLE preprocessing step", action="store_true", default=False)

        ret.append(g)
        g = IgnorableOptionGroup(parser, "Acquisition/Data specific")
        g.add_option("--bat", help="Estimated bolus arrival time (s) - default=0.7 (pASL), 1.3 (cASL)", type=float)
        g.add_option("--batsd", help="Bolus arrival time standard deviation (s)", type=float)
        g.add_option("--t1", help="Tissue T1 (s)", type=float, default=1.3)
        g.add_option("--t1b", help="Blood T1 (s)", type=float, default=1.65)
        ret.append(g)
        g = IgnorableOptionGroup(parser, "Output options")
        g.add_option("--save-corrected", help="Save corrected input data", action="store_true", default=False)
        g.add_option("--save-reg", help="Save registration information (transforms etc)", action="store_true", default=False)
        g.add_option("--save-basil", help="Save Basil modelling output", action="store_true", default=False)
        g.add_option("--save-calib", help="Save calibration output", action="store_true", default=False)
        g.add_option("--save-all", help="Save all output (enabled when --debug specified)", action="store_true", default=False)
        g.add_option("--output-stddev", "--output-std", help="Output standard deviation of estimated variables", action="store_true", default=False)
        g.add_option("--output-var", "--vars", help="Output variance of estimated variables", action="store_true", default=False)
        g.add_option("--output-mni", help="Output in MNI standard space", action="store_true", default=False)
        g.add_option("--no-report", dest="save_report", help="Don't try to generate an HTML report", action="store_false", default=True)
        ret.append(g)
        return ret

def main():
    """
    Entry point for oxasl command line tool
    """
    debug = True
    try:
        parser = AslOptionParser(usage="oxasl -i <asl_image> [options]", version=__version__)
        parser.add_category(image.AslImageOptions())
        parser.add_category(struc.StructuralImageOptions())
        parser.add_category(OxfordAslOptions())
        parser.add_category(calib.CalibOptions(ignore=["perf", "tis"]))
        parser.add_category(reg.RegOptions())
        parser.add_category(corrections.DistcorrOptions())
        if oxasl_ve:
            parser.add_category(oxasl_ve.VeaslOptions())
        if oxasl_enable:
            parser.add_category(oxasl_enable.EnableOptions(ignore=["regfrom",]))
        parser.add_category(GenericOptions())

        options, _ = parser.parse_args()
        debug = options.debug
        if not options.output:
            options.output = "oxasl"
        
        # Some oxasl command-line specific defaults
        if options.calib is not None and options.calib_method is None:
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
        sys.exit(1)

def oxasl(wsp):
    """
    Main oxasl pipeline script
    """
    wsp.log.write("OXASL version: %s\n" % __version__)
    for plugin in (oxasl_ve, oxasl_deblur, oxasl_enable):
        if plugin is not None:
            wsp.log.write(" - Found plugin: %s (version %s)\n" % (plugin.__name__, getattr(plugin, "__version__", "unknown")))

    wsp.log.write("\nInput ASL data: %s\n" % wsp.asldata.name)
    wsp.asldata.summary(wsp.log)

    oxasl_preproc(wsp)
    calib.init(wsp)

    if wsp.asldata.iaf in ("tc", "ct", "diff"):
        model_paired(wsp)
    elif wsp.asldata.iaf == "ve":
        if oxasl_ve is None:
            raise ValueError("Vessel encoded data supplied but oxasl_ve is not installed")
        oxasl_ve.model_ve(wsp)
    else:
        # FIXME support for multiphase data
        raise ValueError("ASL data has format '%s' - not supported by OXASL pipeline" % wsp.asldata.iaf)

    if wsp.save_report:
        do_report(wsp)
    
    do_cleanup(wsp)
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
    try:
        # Not all data can generate a PWI
        img = wsp.asldata.perf_weighted()
        img_type = "Perfusion-weighted image"
    except ValueError:
        img = wsp.asldata.mean()
        img_type = "Mean ASL data"

    page.heading(img_type, level=1)
    page.image("asldata", LightboxImage(img))

def oxasl_preproc(wsp):
    """
    Run standard processing on ASL data

    This method requires wsp to be a Workspace containing certain standard information.
    As a minimum, the attribute ``asldata`` must contain an AslImage object.
    """
    report_asl(wsp)

    struc.init(wsp)
    corrections.apply_corrections(wsp)

    corrections.get_motion_correction(wsp)
    corrections.apply_corrections(wsp)

    reg.reg_asl2calib(wsp)
    reg.reg_asl2struc(wsp, True, False)

    corrections.get_fieldmap_correction(wsp)
    corrections.get_cblip_correction(wsp)
    corrections.get_sensitivity_correction(wsp)
    corrections.apply_corrections(wsp)

    mask.generate_mask(wsp)

    if oxasl_enable and wsp.use_enable:
        wsp.sub("enable")
        oxasl_enable.enable(wsp.enable)
        wsp.corrected.asldata = wsp.enable.asldata_enable

def model_paired(wsp):
    """
    Do model fitting on TC/CT or subtracted data

    Workspace attributes updated
    ----------------------------

     - ``basil``         - Contains model fitting output on data without partial volume correction
     - ``basil_pvcorr``  - Contains model fitting output with partial volume correction if 
                           ``wsp.pvcorr`` is ``True``
     - ``output.native`` - Native (ASL) space output from last Basil modelling output
     - ``output.struc``  - Structural space output
    """
    basil.basil(wsp, output_wsp=wsp.sub("basil"))
    redo_reg(wsp, wsp.basil.finalstep.mean_ftiss)

    wsp.sub("output")

    if wsp.pvcorr:
        # Do partial volume correction fitting
        #
        # FIXME: We could at this point re-apply all corrections derived from structural space?
        # But would need to make sure corrections module re-transforms things like sensitivity map
        #
        # Partial volume correction is very sensitive to the mask, so recreate it
        # if it came from the structural image as this requires accurate ASL->Struc registration
        if wsp.rois.mask_src == "struc":
            wsp.rois.mask_orig = wsp.rois.mask
            wsp.rois.mask = None
            mask.generate_mask(wsp)

        # Generate PVM and PWM maps for Basil
        struc.segment(wsp)
        wsp.structural.wm_pv_asl = reg.struc2asl(wsp, wsp.structural.wm_pv)
        wsp.structural.gm_pv_asl = reg.struc2asl(wsp, wsp.structural.gm_pv)
        wsp.basil_options = {"pwm" : wsp.structural.wm_pv_asl, "pgm" : wsp.structural.gm_pv_asl}
        basil.basil(wsp, output_wsp=wsp.sub("basil_pvcorr"), prefit=False)
        output_native(wsp.output, wsp.basil_pvcorr)
    else:
        output_native(wsp.output, wsp.basil)

    output_trans(wsp.output)

def redo_reg(wsp, pwi):
    """
    Re-do ASL->structural registration using BBR and perfusion image
    """
    wsp.reg.regfrom_orig = wsp.reg.regfrom
    wsp.reg.regto_orig = wsp.reg.regto
    new_regfrom = np.copy(pwi.data)
    new_regfrom[new_regfrom < 0] = 0
    wsp.reg.regfrom = Image(new_regfrom, header=pwi.header)
    wsp.reg.asl2struc_initial = wsp.reg.asl2struc
    wsp.reg.struc2asl_initial = wsp.reg.struc2asl
    reg.reg_asl2struc(wsp, False, True, name="final")

def do_report(wsp):
    """
    Generate HTML report
    """
    report_build_dir = None
    if wsp.debug:
        report_build_dir = os.path.join(wsp.savedir, "report_build")
    wsp.log.write("\nGenerating HTML report\n")
    report_dir = os.path.join(wsp.savedir, "report")
    success = wsp.report.generate_html(report_dir, report_build_dir, log=wsp.log)
    if success:
        wsp.log.write(" - Report generated in %s\n" % report_dir)

OUTPUT_ITEMS = {
    "ftiss" : ("perfusion", 6000, True, "ml/100g/min", "30-50", "10-20"),
    "fblood" : ("aCBV", 100, True, "ml/100g/min", "", ""),
    "delttiss" : ("arrival", 1, False, "s", "", ""),
    "fwm" : ("perfusion_wm", 6000, True, "ml/100g/min", "", "10-20"),
    "deltwm" : ("arrival_wm", 1, False, "s", "", ""),
    "modelfit" : ("modelfit", 1, False, "", "", ""),
    "asldata_diff" : ("asldata_diff", 1, False, "", "", ""),
}
    
def output_native(wsp, basil_wsp, report=None):
    """
    Create native space output images

    :param wsp: Workspace object. Output will be placed in a sub workspace named ``native``
    :param basil_wsp: Workspace in which Basil modelling has been run. The ``finalstep``
                      attribute is expected to point to the final output workspace
    """
    if not wsp.output_native: return

    wsp.sub("native")
    prefixes = ["", "mean"]
    if wsp.output_stddev:
        prefixes.append("std")
    if wsp.output_var:
        prefixes.append("var")
    for fabber_name, oxasl_output in OUTPUT_ITEMS.items():
        for prefix in prefixes:
            is_variance = prefix == "var"
            if is_variance:
                # Variance is not output by Fabber natively so we get it by 
                # squaring the standard deviation. We also pass the flag
                # to the calibration routine so it can square the correction
                # factors
                fabber_output = "std_%s" % fabber_name
            elif prefix:
                fabber_output = "%s_%s" % (prefix, fabber_name)
            else:
                fabber_output = fabber_name
            
            img = basil_wsp.finalstep.ifnone(fabber_output, None)
            if img is not None:
                # Make negative values = 0 and ensure masked value zeroed
                data = np.copy(img.data)
                data[img.data < 0] = 0
                data[wsp.rois.mask.data == 0] = 0
                img = Image(data, header=img.header)
                name, multiplier, calibrate, _, _, _ = oxasl_output
                if prefix and prefix != "mean":
                    name = "%s_%s" % (name, prefix)

                if calibrate:
                    # Anything that needs calibration also requires sensitivity correction
                    img, = corrections.apply_sensitivity_correction(wsp, img)

                if is_variance:
                    img = Image(np.square(img.data), header=img.header)
                setattr(wsp.native, name, img)

                if calibrate and wsp.calib is not None:
                    alpha = wsp.ifnone("calib_alpha", 1.0 if wsp.asldata.iaf == "ve" else 0.85 if wsp.asldata.casl else 0.98)
                    img = calib.calibrate(wsp, img, multiplier=multiplier, alpha=alpha, var=is_variance)
                    name = "%s_calib" % name
                    setattr(wsp.native, name, img)
    
    if wsp.save_mask:
        wsp.native.mask = wsp.rois.mask

    output_report(wsp.native, report=report)
        
def output_report(wsp, report=None):
    """
    Create report pages from output data

    :param wsp: Workspace object containing output
    """
    if report is None:
        report = wsp.report

    roi = wsp.rois.mask.data
    if wsp.structural.struc is not None:
        gm = reg.struc2asl(wsp, wsp.structural.gm_pv).data
        wm = reg.struc2asl(wsp, wsp.structural.wm_pv).data
    
    for oxasl_name, multiplier, calibrate, units, normal_gm, normal_wm in OUTPUT_ITEMS.values():
        name = oxasl_name + "_calib"
        img = getattr(wsp, name)
        if img is None:
            name = oxasl_name
            img = getattr(wsp, name)

        if img is not None and img.ndim == 3:
            page = report.page(name)
            page.heading("Output image: %s" % name)
            if calibrate and name.endswith("_calib"):
                alpha = wsp.ifnone("calib_alpha", 1.0 if wsp.asldata.iaf == "ve" else 0.85 if wsp.asldata.casl else 0.98)
                page.heading("Calibration", level=1)
                page.text("Image was calibrated using supplied M0 image")
                page.text("Inversion efficiency: %f" % alpha)
                page.text("Multiplier for physical units: %f" % multiplier)

            page.heading("Metrics", level=1)
            data = img.data
            table = []
            table.append(["Mean within mask", "%.4g %s" % (np.mean(data[roi > 0.5]), units), ""])
            if wsp.structural.struc is not None:
                table.append(["GM mean", "%.4g %s" % (np.mean(data[gm > 0.5]), units), normal_gm])
                table.append(["Pure GM mean", "%.4g %s" % (np.mean(data[gm > 0.8]), units), normal_gm])
                table.append(["WM mean", "%.4g %s" % (np.mean(data[wm > 0.5]), units), normal_wm])
                table.append(["Pure WM mean", "%.4g %s" % (np.mean(data[wm > 0.9]), units), normal_wm])
            page.table(table, headers=["Metric", "Value", "Typical"])
            
            page.heading("Image", level=1)
            page.image("%s_img" % name, LightboxImage(img, zeromask=False, mask=wsp.rois.mask, colorbar=True))

def output_trans(wsp):
    """
    Create transformed output, i.e. in structural and/or standard space

    FIXME std space not implemented yet
    """
    if not wsp.output_struc or wsp.reg.asl2struc is None:
        return

    wsp.sub("struct")
    for suffix in ("", "_std", "_var", "_calib", "_std_calib", "_var_calib"):
        for output in ("perfusion", "aCBV", "arrival", "perfusion_wm", "arrival_wm", "modelfit", "mask"):
            native_output = getattr(wsp.native, output + suffix)
            # Don't transform 4D output (e.g. modelfit) - too large!
            if native_output is not None and native_output.ndim == 3:
                if wsp.reg.asl2struc is not None:
                    setattr(wsp.struct, output + suffix, reg.asl2struc(wsp, native_output, mask=(output == 'mask')))

def do_cleanup(wsp):
    """
    Remove items from the workspace that are not being output. The
    corresponding files will be deleted.
    """
    wsp.log.write("\nDoing cleanup\n")
    if not wsp.save_all:
        if not wsp.save_corrected:
            wsp.log.write(" - Removing corrected data\n")
            wsp.corrected = None
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
        wsp.input = None
        wsp.rois = None
