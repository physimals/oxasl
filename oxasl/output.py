"""
OXASL - Output module

Copyright (c) 2008-2020 Univerisity of Oxford
"""
import itertools

import numpy as np

from fsl.data.image import Image
from fsl.data.atlases import AtlasRegistry

from oxasl.options import OptionCategory, OptionGroup
from oxasl import reg, calibration
from oxasl.reporting import LightboxImage

class Options(OptionCategory):
    """
    OptionCategory which contains options for output files to keep
    """

    def __init__(self, **kwargs):
        OptionCategory.__init__(self, "output")

    def groups(self, parser):
        groups = []

        g = OptionGroup(parser, "Output options")
        g.add_option("--save-input", help="Save unprocessed input data", action="store_true", default=False)
        g.add_option("--save-preproc", help="Save preprocessed input data", action="store_true", default=False)
        g.add_option("--save-corrected", help="Save corrected input data (distortion/moco/sensitivity corrected)", action="store_true", default=False)
        g.add_option("--save-corrections", help="Save intermediate correction data (e.g. fieldmaps, moco)", action="store_true", default=False)
        g.add_option("--save-reg", help="Save registration information (transforms etc)", action="store_true", default=False)
        g.add_option("--save-filter", help="Save data filtering output", action="store_true", default=False)
        g.add_option("--save-quantification", help="Save quantification output", action="store_true", default=False)
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
    "noise_means" : ("noise", 1, False, "", "", ""),
    "modelfit" : ("modelfit", 1, False, "", "", ""),
    "modelfit_mean" : ("modelfit_mean", 1, False, "", "", ""),
    "residuals" : ("residuals", 1, False, "", "", ""),
    "asldata_diff" : ("asldata_diff", 1, False, "", "", ""),
    "T_exch" : ("texch", 1, False, "", "", ""),
}

def run(wsp):
    wsp.log.write("\nGenerating output images\n")
    for quantify_wsp in wsp.quantify_wsps:
        try:
            quantify_name = quantify_wsp[quantify_wsp.index("_"):]
        except ValueError:
            quantify_name = ""
        output_wsp = "output%s" % quantify_name
        _output_from_quantification(getattr(wsp, quantify_wsp), wsp.sub(output_wsp))

def _output_from_quantification(quantify_wsp, output_wsp):
    """
    """
    if quantify_wsp.image_space is None:
        _output_native(output_wsp.sub("native"), quantify_wsp)
        _output_trans(output_wsp)
    else:
        _output_native(output_wsp, quantify_wsp)

def _output_native(wsp, quantify_wsp):
    """
    Create output images from quantification

    This includes basic sanity-processing (removing negatives and NaNs) plus
    calibration using existing M0

    :param wsp: Workspace object for output
    :param quantify_wsp: Workspace in which quantification has been run. The ``finalstep``
                      attribute is expected to point to the final output workspace
    """
    wsp.log.write(" - Generating native (ASL) space output\n")
    # Output the differenced data averaged across repeats for kinetic curve comparison
    # with the model
    if wsp.asldata.iaf in ("tc", "ct", "diff"):
        wsp.diffdata_mean = wsp.asldata.diff().mean_across_repeats()

    # Save the analysis mask
    wsp.mask = quantify_wsp.analysis_mask

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

            img = quantify_wsp.finalstep.ifnone(fabber_output, None)
            if img is not None:
                # Make negative/nan values = 0 and ensure masked value zeroed
                # since quantification may use a different fitting mask to the pipeline mask
                data = np.copy(img.data)
                data[~np.isfinite(data)] = 0
                data[img.data < 0] = 0
                data[wsp.mask.data == 0] = 0
                img = Image(data, header=img.header)
                name, multiplier, calibrate, units, normal_gm, normal_wm = oxasl_output
                if prefix and prefix != "mean":
                    name = "%s_%s" % (name, prefix)

                # Variance is not output by Fabber natively so we get it by
                # squaring the standard deviation.
                if is_variance:
                    img = Image(np.square(img.data), header=img.header)
                setattr(wsp, name, img)

                if calibrate:
                    for method in wsp.calibration.calib_method:
                        if method != "prequantified":
                            calib_wsp = getattr(wsp.calibration, method)
                            img_calib = calibration.run(calib_wsp, img, multiplier=multiplier, var=is_variance)
                            sub_wsp_name = "calib_%s" % method
                            calib_output_wsp = getattr(wsp, sub_wsp_name)
                            if calib_output_wsp is None:
                                calib_output_wsp = wsp.sub(sub_wsp_name)
                            setattr(calib_output_wsp, name, img_calib)
                            output_report(calib_output_wsp, name, units, normal_gm, normal_wm, method)
                        else:
                            output_report(wsp, name, units, normal_gm, normal_wm, method)
                else:
                    output_report(wsp, name, units, normal_gm, normal_wm)

def output_report(wsp, name, units, normal_gm, normal_wm, calib_method="none"):
    """
    Create report pages from output data

    :param wsp: Workspace object containing output
    """
    report = wsp.report

    img = getattr(wsp, name)
    if img is not None and img.ndim == 3:
        page = report.page("%s_%s" % (name, calib_method))
        page.heading("Output image: %s (calibration: %s)" % (name, calib_method))
        if calib_method != "none":
            alpha = wsp.ifnone("calib_alpha", 1.0 if wsp.asldata.iaf in ("ve", "vediff") else 0.85 if wsp.asldata.casl else 0.98)
            page.heading("Calibration", level=1)
            page.text("Calibration method: %s" % calib_method)
            page.text("Inversion efficiency: %f" % alpha)

        page.heading("Metrics", level=1)
        data = img.data
        roi = reg.change_space(wsp, wsp.mask, img).data
        table = []
        table.append(["Mean within mask", "%.4g %s" % (np.mean(data[roi > 0.5]), units), ""])
        if wsp.structural.struc is not None:
            if wsp.structural.gm_pv_asl is None:
                wsp.structural.gm_pv_asl = reg.change_space(wsp, wsp.structural.gm_pv, img).data
                wsp.structural.wm_pv_asl = reg.change_space(wsp, wsp.structural.wm_pv, img).data
            gm = np.asarray(wsp.structural.gm_pv_asl.data)
            wm = np.asarray(wsp.structural.wm_pv_asl.data)
            if wsp.structural.cortex_asl is None:
                wsp.structural.cortex_asl = _get_cortex(wsp)
            cortex = wsp.structural.cortex_asl
            some_gm, some_wm = gm > 0.5, wm > 0.5
            pure_gm, pure_wm = gm > wsp.ifnone("gm_thresh", 0.8), wm > wsp.ifnone("wm_thresh", 0.9)
            table.append(["GM mean", "%.4g %s" % (np.mean(data[some_gm]), units), normal_gm])
            table.append(["Pure GM mean", "%.4g %s" % (np.mean(data[pure_gm]), units), normal_gm])
            table.append(["Cortical GM mean", "%.4g %s" % (np.mean(data[np.logical_and(pure_gm, cortex)]), units), normal_gm])
            table.append(["WM mean", "%.4g %s" % (np.mean(data[some_wm]), units), normal_wm])
            table.append(["Pure WM mean", "%.4g %s" % (np.mean(data[pure_wm]), units), normal_wm])
            table.append(["Cerebral wM mean", "%.4g %s" % (np.mean(data[np.logical_and(pure_wm, cortex)]), units), normal_wm])

        page.table(table, headers=["Metric", "Value", "Typical"])

        page.heading("Image", level=1)
        page.image("%s_%s_img" % (name, calib_method), LightboxImage(img, zeromask=False, mask=wsp.mask, colorbar=True))

def _get_cortex(wsp):
    atlases = AtlasRegistry()
    atlases.rescanAtlases()
    atlas = atlases.loadAtlas("harvardoxford-subcortical", loadSummary=False, resolution=2)
    cort_std = Image(np.mean(atlas.data[..., 0:2] + atlas.data[..., 11:13], axis=-1) * 2, header=atlas.header)
    cort_asl = reg.change_space(wsp, cort_std, "asl").data
    return cort_asl > 50

def __output_trans_helper(wsp):
    """
    Generator to provide all the combinations of variables for the output_
    trans() function. Note that 4D data and upwards will be skipped. 
    """

    # We loop over these combinations of variables for each output space
    # (structural, standard, custom)
    suffixes = ("", "_std", "_var")
    outputs = (
        "perfusion", "aCBV", "arrival", "perfusion_wm", 
        "arrival_wm", "modelfit", "modelfit_mean", 
        "residuals", "texch", "mask"
    )

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
        wsp.log.write(" - Generating output in %s space\n" % name)
        output_wsp = wsp.sub(space)
        for suffix, output, native_data in __output_trans_helper(wsp): 
            setattr(output_wsp, output + suffix, reg.change_space(wsp, native_data, space, mask=(output == 'mask')))
            for method in wsp.calibration.calib_method:
                if method != "prequantified":
                    sub_wsp_name = "calib_%s" % method
                    native_calib_output_wsp = getattr(wsp.native, sub_wsp_name)
                    native_calib_data = getattr(native_calib_output_wsp, output + suffix)
                    calib_output_wsp = output_wsp.sub(sub_wsp_name)
                    setattr(calib_output_wsp, output + suffix, reg.change_space(wsp, native_calib_data, space, mask=(output == 'mask')))
                else:
                    native_data = getattr(wsp.native, output + suffix)
                    setattr(output_wsp, output + suffix, reg.change_space(wsp, native_data, space, mask=(output == 'mask')))
