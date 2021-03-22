"""
OXASL - Output module

Copyright (c) 2008-2020 Univerisity of Oxford
"""
import itertools

import numpy as np

from fsl.data.image import Image

from oxasl.options import OptionCategory, OptionGroup
from oxasl import reg, corrections, calibration
from oxasl.reporting import LightboxImage

class Options(OptionCategory):
    """
    OptionCategory which contains options for registration of ASL data to structural image
    """

    def __init__(self, **kwargs):
        OptionCategory.__init__(self, "output")

    def groups(self, parser):
        groups = []

        g = OptionGroup(parser, "Output options")
        g.add_option("--save-corrected", help="Save corrected input data", action="store_true", default=False)
        g.add_option("--save-reg", help="Save registration information (transforms etc)", action="store_true", default=False)
        g.add_option("--save-basil", help="Save Basil modelling output", action="store_true", default=False)
        g.add_option("--save-calib", help="Save calibration output", action="store_true", default=False)
        g.add_option("--save-all", help="Save all output (enabled when --debug specified)", action="store_true", default=False)
        g.add_option("--output-stddev", "--output-std", help="Output standard deviation of estimated variables", action="store_true", default=False)
        g.add_option("--output-var", "--vars", help="Output variance of estimated variables", action="store_true", default=False)
        g.add_option("--output-residuals", help="Output residuals (model fit - actual data)", action="store_true", default=False)
        g.add_option("--output-mni", help="Output in MNI standard space", action="store_true", default=False)
        g.add_option("--output-custom", help="Output in custom space (provide path to reference image in space)", type=str)
        g.add_option("--output-custom-mat", help="(Optional) FLIRT transformation from structural space to custom space. " +
                        "If not provided, will FLIRT registration from structural to --output-custom will be used.", type=str)
        g.add_option("--no-report", dest="save_report", help="Don't try to generate an HTML report", action="store_false", default=True)
        groups.append(g)

        return groups

OUTPUT_ITEMS = {
    "ftiss" : ("perfusion", 6000, True, "ml/100g/min", "30-50", "10-20"),
    "fblood" : ("aCBV", 100, True, "ml/100g/min", "", ""),
    "delttiss" : ("arrival", 1, False, "s", "", ""),
    "fwm" : ("perfusion_wm", 6000, True, "ml/100g/min", "", "10-20"),
    "deltwm" : ("arrival_wm", 1, False, "s", "", ""),
    "modelfit" : ("modelfit", 1, False, "", "", ""),
    "modelfit_mean" : ("modelfit_mean", 1, False, "", "", ""),
    "residuals" : ("residuals", 1, False, "", "", ""),
    "asldata_diff" : ("asldata_diff", 1, False, "", "", ""),
    "T_exch" : ("texch", 1, False, "", "", ""),
}

def run(wsp):
    for basildir in wsp.basildirs:
        if basildir:
            basil_wsp = "basil_%s" % basildir
            output_wsp = "output_%s" % basildir
        else:
            basil_wsp = "basil"
            output_wsp = "output"
        _output_basil(getattr(wsp, basil_wsp), wsp.sub(output_wsp), basildir)

def _output_basil(basil_wsp, output_wsp, basildir):
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
    if basil_wsp.image_space is None:
        _output_native(output_wsp.sub("native"), basil_wsp, basildir)
        _output_trans(output_wsp)
    else:
        _output_native(output_wsp, basil_wsp, basildir)

def _output_native(wsp, basil_wsp, basildir, report=None):
    """
    Create output images from a Basil run

    This includes basic sanity-processing (removing negatives and NaNs) plus
    calibration using existing M0

    :param wsp: Workspace object for output
    :param basil_wsp: Workspace in which Basil modelling has been run. The ``finalstep``
                      attribute is expected to point to the final output workspace
    """
    # Output the differenced data averaged across repeats for kinetic curve comparison
    # with the model
    if wsp.asldata.iaf in ("tc", "ct", "diff"):
        wsp.diffdata_mean = wsp.asldata.diff().mean_across_repeats()

    # Output model fitting results
    prefixes = ["", "mean"]
    if wsp.output_stddev:
        prefixes.append("std")
    if wsp.output_var or wsp.region_analysis:
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
                # Make negative/nan values = 0 and ensure masked value zeroed
                # since Basil may use a different fitting mask to the pipeline mask
                data = np.copy(img.data)
                data[~np.isfinite(data)] = 0
                data[img.data < 0] = 0
                mask = reg.change_space(wsp, wsp.rois.mask, img)
                data[mask.data == 0] = 0
                img = Image(data, header=img.header)
                name, multiplier, calibrate, _, _, _ = oxasl_output
                if prefix and prefix != "mean":
                    name = "%s_%s" % (name, prefix)

                #if calibrate:
                #    # Anything that needs calibration also requires sensitivity correction
                #    img, = corrections.apply_sensitivity_correction(wsp, img)

                if is_variance:
                    img = Image(np.square(img.data), header=img.header)
                setattr(wsp, name, img)

                if calibrate and wsp.calib is not None:
                    img = calibration.run(wsp, img, multiplier=multiplier, var=is_variance)
                    name = "%s_calib" % name
                    setattr(wsp, name, img)

    if wsp.save_mask:
        wsp.mask = wsp.rois.mask

    output_report(wsp, basildir, report=report)

def output_report(wsp, basildir, report=None):
    """
    Create report pages from output data

    :param wsp: Workspace object containing output
    """
    if report is None:
        report = wsp.report

    roi, gm, wm = None, None, None
    for oxasl_name, multiplier, calibrate, units, normal_gm, normal_wm in OUTPUT_ITEMS.values():
        name = oxasl_name + "_calib"
        img = getattr(wsp, name)
        if img is None:
            name = oxasl_name
            img = getattr(wsp, name)

        if img is not None and img.ndim == 3:
            if basildir:
                page_name = "%s_%s" % (name, basildir)
            else:
                page_name = name
            page = report.page(page_name)
            page.heading("Output image: %s" % name)
            if calibrate and name.endswith("_calib"):
                alpha = wsp.ifnone("calib_alpha", 1.0 if wsp.asldata.iaf in ("ve", "vediff") else 0.85 if wsp.asldata.casl else 0.98)
                page.heading("Calibration", level=1)
                page.text("Image was calibrated using supplied M0 image")
                page.text("Inversion efficiency: %f" % alpha)
                page.text("Multiplier for physical units: %f" % multiplier)

            page.heading("Metrics", level=1)
            data = img.data
            if roi is None:
                roi = reg.change_space(wsp, wsp.rois.mask, img).data
            table = []
            table.append(["Mean within mask", "%.4g %s" % (np.mean(data[roi > 0.5]), units), ""])
            if wsp.structural.struc is not None:
                if gm is None:
                    gm = reg.change_space(wsp, wsp.structural.gm_pv, img).data
                    wm = reg.change_space(wsp, wsp.structural.wm_pv, img).data
                table.append(["GM mean", "%.4g %s" % (np.mean(data[gm > 0.5]), units), normal_gm])
                table.append(["Pure GM mean", "%.4g %s" % (np.mean(data[gm > 0.8]), units), normal_gm])
                table.append(["WM mean", "%.4g %s" % (np.mean(data[wm > 0.5]), units), normal_wm])
                table.append(["Pure WM mean", "%.4g %s" % (np.mean(data[wm > 0.9]), units), normal_wm])
            page.table(table, headers=["Metric", "Value", "Typical"])

            page.heading("Image", level=1)
            page.image("%s_img" % name, LightboxImage(img, zeromask=False, mask=wsp.rois.mask, colorbar=True))

def __output_trans_helper(wsp):
    """
    Generator to provide all the combinations of variables for the output_
    trans() function. Note that 4D data and upwards will be skipped. 
    """

    # We loop over these combinations of variables for each output space
    # (structural, standard, custom)
    suffixes = ("", "_std", "_var", "_calib", "_std_calib", "_var_calib")
    outputs = ("perfusion", "aCBV", "arrival", "perfusion_wm", 
        "arrival_wm", "modelfit", "modelfit_mean", "residuals", "texch", "mask")

    for suff, out in itertools.product(suffixes, outputs):
        data = getattr(wsp.native, out + suff)
        # Don't transform 4D output (e.g. modelfit) - too large!
        if (data is not None) and data.ndim == 3:
            yield suff, out, data

def _output_trans(wsp):
    """
    Create transformed output, i.e. in structural and/or standard space
    """
    output_spaces = []
    if wsp.output_struc and wsp.reg.asl2struc is not None:
        output_spaces.append(("struc", "structural"))
    
    if wsp.output_mni:
        if wsp.reg.struc2asl is None:
            wsp.log.write(" - WARNING: No structural registration - cannot output in standard space\n")
        else:
            output_spaces.append(("std", "standard (MNI)"))
    
    if wsp.output_custom:
        output_spaces.append(("custom", "user-defined custom"))

    for space, name in output_spaces:
        wsp.log.write("\nGenerating output in %s space\n" % name)
        output_wsp = wsp.sub(space)
        for suffix, output, native_output in __output_trans_helper(wsp): 
            setattr(output_wsp, output + suffix, reg.change_space(wsp, native_output, space, mask=(output == 'mask')))
        wsp.log.write(" - DONE\n")
