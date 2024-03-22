"""
Volumetric quantification using SVB
"""
import os.path as op

import numpy as np
from fsl.data.image import Image

try:
    import regtricks as rt
    from ssvb import SSVBFit
    from ssvb.data import DataModel
    from ssvb.main import run_inference
    from ssvb.structure import Cortex, Volumetric, VoxelROI
except:
    SSVBFit = None

TISSUE_PROPERTIES = {
    "GM": {"att": 1.3, "pc": 0.98, "fcalib": 0.01, "t1": 1.3},
    "WM": {"att": 1.6, "pc": 0.8, "fcalib": 0.003, "t1": 1.1},
    "mixed": {"artatt": 0.9, "att": 1.3, "pc": 0.9, "fcalib": 0.01, "t1": 1.3},
    "art": {"att": 0.9, "pc": 0.9, "fcalib": 0.01, "t1": 1.3},
}

def _get_default_tissue_properties(data_model, mode="hybrid"):
    overrides = {}
    for s in data_model.structures:
        if mode == "hybrid":
            if type(s) is Volumetric:
                # WM volumetric structure
                overrides[s.name] = TISSUE_PROPERTIES["WM"]
            else:
                # Cortex and single-node subcortical ROIs
                overrides[s.name] = TISSUE_PROPERTIES["GM"]
        elif mode == "volume":
            if isinstance(s, Volumetric) and (len(data_model.structures) == 1):
                overrides[s.name] = TISSUE_PROPERTIES["mixed"]
        elif mode == "volume_pvec":
            overrides[s.name] = TISSUE_PROPERTIES[s.name]
        else:
            raise RuntimeError("unrecognised mode")
    return overrides

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
    if SSVBFit is None:
        raise ValueError("SSVB is not installed or could not be imported")

    wsp.log.write("\nRunning SVB modelling on ASL data\n")

    wsp = wsp.sub("ssvb")
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

    data = wsp.asldata_diff_scaled.data
    if wsp.pvcorr:
        wsp.log.write(" - Using separate GM/WM modelling\n")
        data_model = DataModel.volumetric_pvec(data, {
            "GM": wsp.structural.gm_pv_asl.data,
            "WM": wsp.structural.wm_pv_asl.data,
        })
        mode = "volume_pvec"
    else:
        wsp.log.write(" - Combined GM/WM modelling\n")
        data_model = DataModel.volumetric(data, wsp.analysis_mask.data)
        mode = "volume"

    # Hack
    data_model.nii = wsp.asldata_diff_scaled.nibImage

    # Base options for use in all fits
    overrides = _get_default_tissue_properties(data_model, mode)

    # FIXME this is really ugly. Because we want ATT to be an inference parameter,
    # we have to re-create its entry in the overrides dict to define the distribution
    # type and mean. It would be nicer if the some later code could pick up a the
    # numerical value and use that automatically as the mean of the corresponding post / prior ?
    step1_overrides = {}
    for k in overrides.keys():
        new_att = {
            "att": {
                "prior_mean": overrides[k].get("att"),
            },
        }
        step1_overrides[k] = {**overrides[k], **new_att}

    print(step1_overrides)
    model_options = dict(
        tau=wsp.asldata.taus,
        tis=wsp.asldata.tis,
        taus=wsp.asldata_diff.taus,
        repeats=wsp.asldata.rpts,
        casl=wsp.asldata.casl,
        prior_dist="N",
        post_dist="F",
    )

    # Step 1: CBF only
    data_model.set_fwd_model(
        "aslrest",
        overrides=step1_overrides,
        infer_cbf=True,
        infer_att=wsp.inferbat,
        **model_options
    )

    fit = SSVBFit(data_model)
    fit_options = {
        "debug": True,
        "display_step": wsp.ifnone("svb_display_step", 50),
        "batch_size": wsp.asldata.ntis,
        "learning_rate": wsp.ifnone("svb_learning_rate", 0.1),
        "epochs" : wsp.ifnone("svb_epochs", -5),
        "lr_decay_rate": wsp.ifnone("svb_lr_decay_rate", 1.0),
        "sample_size": wsp.ifnone("svb_sample_size", 4),
        "ss_increase_rate": wsp.ifnone("svb_ss_increase_rate", 2),
    }
    outdir = "svb_123_out"
    runtime, fit, training_history = run_inference(
        fit, op.join(outdir, "step1"), **fit_options
    )

    step2_overrides = {}
    for k in step1_overrides.keys():
        new_att = {
            "artcbf": {
                "prior_dist": "A",
            },
            "artatt": {
                "prior_mean": overrides[k].get("artatt"),
            }
        }
        step2_overrides[k] = {**step1_overrides[k], **new_att}

    if wsp.inferart:
        # Step 2: Arterial
        data_model.set_fwd_model(
            "aslrest",
            overrides=step2_overrides,
            infer_cbf=True,
            infer_att=wsp.inferbat,
            infer_artcbf=wsp.inferart,
            infer_artatt=wsp.inferart,
            **model_options
        )
        fit2 = SSVBFit(data_model)
        fit2.continue_from(fit)
        runtime, fit, training_history = run_inference(
            fit2, op.join(outdir, "step2"), **fit_options
        )
        fit = fit2

    # Step 3: spatial
    step3_overrides = {}
    for k in step2_overrides.keys():
        new_att = {
            "cbf": {
                "prior_dist": "M",
            },
            "att": {
                "prior_dist": "NM",
            },
            "artcbf": {
                "prior_dist": "A",
            },
            "artatt": {
                "prior_dist": "NM",
            }
        }
        step3_overrides[k] = {**step2_overrides[k], **new_att}

    data_model.set_fwd_model(
        "aslrest",
        overrides=step3_overrides,
        infer_cbf=True,
        infer_att=wsp.inferbat,
        infer_artcbf=wsp.inferart,
        infer_artatt=wsp.inferart,
        **model_options
    )
    fit3 = SSVBFit(data_model)
    fit3.continue_from(fit)
    runtime, fit, training_history = run_inference(
        fit3, op.join(outdir, "step3"), **fit_options
    )

    wsp_output = wsp.sub("finalstep")
    ftiss = Image(op.join(outdir, "step3", "mean_cbf_voxels.nii.gz"))
    wsp_output.mean_ftiss = Image(ftiss.data * wsp.scaling_factor, header=ftiss.header)
    fblood = Image(op.join(outdir, "step3", "mean_artcbf_voxels.nii.gz"))
    wsp_output.mean_fblood = Image(fblood.data * wsp.scaling_factor, header=ftiss.header)
    modelfit = Image(op.join(outdir, "step3", "modelfit.nii.gz"))
    wsp_output.modelfit = Image(modelfit.data * wsp.scaling_factor, header=ftiss.header)
    wsp_output.mean_delttiss = Image(op.join(outdir, "step3", "mean_att_voxels.nii.gz"))
    wsp_output.mean_deltblood = Image(op.join(outdir, "step3", "mean_artatt_voxels.nii.gz"))
    wsp.quantify_wsps.append("ssvb")
