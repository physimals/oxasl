"""
Volumetric quantification using SVB
"""
import os.path as op

import numpy as np
from fsl.data.image import Image

try:
    from vaby import InferenceMethod, DataModel, NP_DTYPE
    from vaby.utils import setup_logging
except ImportError:
    DataModel = None

def _set_defaults(wsp):
    # Single or Multi TI setup
    if wsp.asldata.ntis == 1:
        # Single TI data - don't try to infer arterial component of bolus duration, we don't have enough info
        wsp.log.write(" - Operating in Single TI mode - no arterial component, fixed bolus duration\n")
        wsp.inferart = False
        wsp.infertau = False
        batsd_default = 0.1
    else:
        # For multi TI/PLD data, set a more liberal prior for tissue ATT since we should be able to
        # determine this from the data. NB this leaves the arterial BAT alone.
        batsd_default = 1

    if wsp.wp:
        # White paper mode - this overrides defaults, but can be overwritten by command line
        # specification of individual parameters
        wsp.log.write(" - Analysis in white paper mode: T1 default=1.65, BAT default=0, voxelwise calibration\n")
        t1_default = 1.65
        bat_default = 0.0
    else:
        t1_default = 1.3
        if wsp.asldata.casl:
            bat_default = 1.3
        else:
            bat_default = 0.7

    wsp.t1 = wsp.ifnone("t1", t1_default)
    wsp.t1b = wsp.ifnone("t1b", 1.65)
    wsp.bat = wsp.ifnone("bat", bat_default)
    wsp.batsd = wsp.ifnone("batsd", batsd_default)
    wsp.infertiss = wsp.ifnone("infertiss", True)
    wsp.inferbat = wsp.ifnone("inferbat", len(wsp.asldata.tis) > 1)
    wsp.inferart = wsp.ifnone("inferart", len(wsp.asldata.tis) > 1)

    # if we are doing CASL then fix the bolus duration, unless explicitly told us otherwise
    wsp.infertau = wsp.ifnone("infertau", not wsp.asldata.casl)

    # Take a copy of extra user-specified BASIL options
    wsp.svb_options = dict(wsp.ifnone("svb_options", {}))

def _define_mask(wsp):
    # Save a copy of the *analysis* mask to be used for this data as this might change later
    # e.g. for PVC. Fitting may occur in a different mask. Two locations for compatibility...
    if wsp.rois is not None and wsp.rois.mask is not None:
        wsp.analysis_mask = wsp.rois.mask
    else:
        wsp.analysis_mask = wsp.mask

def run(wsp):
    if DataModel is None:
        raise ValueError("VABY is not installed or could not be imported")

    wsp.log.write("\nRunning SVB modelling on ASL data\n")
    wsp = wsp.sub("vaby")
    _set_defaults(wsp)
    _define_mask(wsp)
    wsp.asldata_diff = wsp.asldata.diff().reorder("tr")
    wsp.log.write(" - Scaling ASL data\n")
    dgmp = np.percentile(wsp.asldata_diff.data, 75)
    wsp.log.write(f" - Data 75th percentile: {dgmp}\n")
    wsp.scaling_factor = 10.0 ** (np.floor(np.log10(dgmp)) + 1)
    wsp.log.write(f" - Scaling factor: {wsp.scaling_factor}\n")
    scaled_data = wsp.asldata_diff.data / wsp.scaling_factor
    wsp.asldata_diff_scaled = wsp.asldata_diff.derived(scaled_data)

    gm_att_prior = 1.3
    wm_att_prior = 1.6
    data_config = {
        "data" : wsp.asldata_diff_scaled.nibImage,
        "mask" : wsp.analysis_mask.nibImage
    }

    model_config = {
        "name" : "aslrest",
        "casl" : True,
        "tis" : wsp.asldata_diff.tis,
        "taus" : wsp.asldata_diff.taus,
        "repeats" : wsp.asldata_diff.rpts,
        "inferart" : False,
        "att_init" : "max",
        "post_dist" : "FoldedNormal",
        "params" : {
            "cbf" : {},
            "att" : {},
            "abv" : {
                "prior_dist" : "ARD",
            },
            "artt" : {
                "post_dist" : "LogNormal",
            },
        },
    }

    if wsp.pvcorr:
        structure_config = {
            "gm" : {
                "type" : "WeightedVolume",
                "weights" : wsp.structural.gm_pv_asl.data,
                "model" : {
                    "att" : gm_att_prior,
                    "t1" : 1.3,
                    "pc" : 0.98,
                    "fcalib" : 0.01,
                    "params" : {}
                }
            },
            "wm" : {
                "type" : "WeightedVolume",
                "weights" : wsp.structural.wm_pv_asl.data,
                "model" : {
                    "att" : wm_att_prior,
                    "t1" : 1.1,
                    "pc" : 0.8,
                    "fcalib" : 0.003,
                    "params" : {}
                }
            }
        }
    else:
        structure_config = {
            "voxels" : {
                "type" : "WeightedVolume",
                "model" : {
                    "att" : gm_att_prior,
                    "params" : {}
                },
            },
        }

    infer_config = {
        "method" : "svb",
        "max_epochs" : None,
        "display_step" : wsp.ifnone("svb_display_step", 50),
        "learning_rate" : wsp.ifnone("svb_learning_rate", 0.1),
        "batch_size" : len(wsp.asldata_diff.tis),
        "sample_size" : wsp.ifnone("svb_sample_size", 4),
        "sample_size_increase_rate": wsp.ifnone("svb_ss_increase_rate", 2),
        "log_avg" : "mean",
        "debug" : False,
    }

    # Step 1 - no arterial component
    outdir = op.join(wsp.savedir, "step1")
    setup_logging(clear=True, log_level="INFO", log_stream=wsp.log, save_log=True, outdir=outdir)
    data_model = DataModel(data_config, model_config, structure_config)
    inference_method = InferenceMethod.get_inference_method(data_model, **infer_config)
    state = inference_method.run()
    inference_method.save(state, outdir=outdir, save_model_fit=True, save_mean=True, save_input_data=True, save_pv=True, save_var=True, save_src_space=True)

    # Step 2 - tissue + arterial
    if wsp.inferart:
        outdir = op.join(wsp.savedir, "step2")
        setup_logging(clear=True, log_level="INFO", log_stream=wsp.log, save_log=True, outdir=outdir)
        for struc_name in structure_config:
            mean = state[struc_name]["mean"]
            var = state[struc_name]["var"]
            struc_model_config = structure_config[struc_name]["model"]
            for idx, param_name in enumerate(["cbf", "att"]):
                if param_name not in struc_model_config["params"]:
                    struc_model_config["params"][param_name] = {}
                struc_model_config["params"][param_name]["post_mean"] = mean[idx]
                struc_model_config["params"][param_name]["post_var"] = var[idx]

        weights = np.zeros(data_model.srcdata.shape, dtype=NP_DTYPE)
        weights[data_model.srcdata.mask > 0] = data_model.net_pv
        structure_config["art"] = {
            "type" : "WeightedVolume",
            "weights" : weights,
            "model" : {
                "artonly" : True,
                "inferart" : True,
                "params" : {},
            }
        }
        data_model = DataModel(data_config, model_config, structure_config)
        inference_method = InferenceMethod.get_inference_method(data_model, **infer_config)
        state = inference_method.run()
        inference_method.save(state, outdir=outdir, save_model_fit=True, save_mean=True, save_input_data=True, save_pv=True, save_var=True, save_src_space=True)

    # Step 3 - spatial
    outdir = op.join(wsp.savedir, "step3")
    setup_logging(clear=True, log_level="INFO", log_stream=wsp.log, save_log=True, outdir=outdir)
    for struc_name in structure_config:
        mean = state[struc_name]["mean"]
        var = state[struc_name]["var"]
        struc_model_config = structure_config[struc_name]["model"]
        if struc_name != "art":
            for idx, param_name in enumerate(["cbf", "att"]):
                if param_name not in struc_model_config["params"]:
                    struc_model_config["params"][param_name] = {}
                struc_model_config["params"][param_name]["post_mean"] = mean[idx]
                struc_model_config["params"][param_name]["post_var"] = var[idx]
        else:
            for idx, param_name in enumerate(["abv", "artt"]):
                if param_name not in struc_model_config["params"]:
                    struc_model_config["params"][param_name] = {}
                struc_model_config["params"][param_name]["post_mean"] = mean[idx]
                struc_model_config["params"][param_name]["post_var"] = var[idx]

    model_config["params"]["cbf"]["prior_dist"] = "MRFSpatial"
    model_config["params"]["att"]["prior_dist"] = "MRFCombined"
    model_config["params"]["artt"]["prior_dist"] = "MRFCombined"
    data_model = DataModel(data_config, model_config, structure_config)
    inference_method = InferenceMethod.get_inference_method(data_model, **infer_config)
    state = inference_method.run()
    inference_method.save(state, outdir=outdir, save_model_fit=True, save_mean=True, save_input_data=True, save_pv=True, save_var=True, save_src_space=True)

    wsp_output = wsp.sub("finalstep")
    ftiss = Image(op.join(outdir, "mean_cbf_combined.nii.gz"))
    wsp_output.mean_ftiss = ftiss.data * wsp.scaling_factor
    fblood = Image(op.join(outdir, "mean_abv_combined.nii.gz"))
    wsp_output.mean_fblood = fblood.data * wsp.scaling_factor
    modelfit = Image(op.join(outdir, "modelfit.nii.gz"))
    wsp_output.modelfit = modelfit.data * wsp.scaling_factor
    wsp_output.mean_delttiss = Image(op.join(outdir, "mean_att_combined.nii.gz"))
    wsp_output.mean_deltblood = Image(op.join(outdir, "mean_artt_combined.nii.gz"))
    wsp.quantify_wsps.append("vaby")
