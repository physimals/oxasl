"""
OXASL: Multi-step Bayesian model fitting for ASL data

Based on the BASIL shell command from FSL
"""
import sys
import math

import numpy as np
import scipy.ndimage

from fsl.wrappers import LOAD
from fsl.data.image import Image

from oxasl import __version__, __timestamp__, reg

def run(wsp, prefit=True, **kwargs):
    """
    Run multi-step model fitting on ASL data in a workspace

    :param wsp: Workspace object
    :param prefit: If True, run a pre-fitting step using the mean over repeats of the ASL data

    Required workspace attributes
    -----------------------------

     - ``asldata`` : AslImage object

    Optional workspace attributes
    -----------------------------

     - ``mask`` : Brain mask (fsl.Image)
     - ``wp`` : If True, use 'white paper' mode (Alsop et al) - modifies some defaults and infers tissue component only
     - ``infertiss`` : If True, infer tissue component (default: True)
     - ``inferbat`` : If True, infer bolus arrival time (default: False)
     - ``infertau`` : If True, infer bolus duration (default: False)
     - ``inferart`` : If True, infer arterial component (default: False)
     - ``infert1`` : If True, infer T1 (default: False)
     - ``inferpc`` : If True, infer PC (default: False)
     - ``t1``: Assumed/initial estimate for tissue T1 (default: 1.65 in white paper mode, 1.3 otherwise)
     - ``t1b``: Assumed/initial estimate for blood T1 (default: 1.65)
     - ``bat``: Assumed/initial estimate for bolus arrival time (s) (default 0 in white paper mode, 1.3 for CASL, 0.7 otherwise)
     - ``t1im`` : T1 map as Image
     - ``pgm`` :  Grey matter partial volume map as Image
     - ``pwm`` : White matter partial volume map as Image
     - ``initmvn`` : MVN structure to use as initialization as Image
     - ``spatial`` : If True, include final spatial VB step (default: False)
     - ``onestep`` : If True, do all inference in a single step (default: False)
     - ``basil_options`` : Optional dictionary of additional options for underlying model
    """
    wsp.log.write("\nRunning BASIL Bayesian modelling on ASL data in '%s' data space\n" % wsp.ifnone("image_space", "asl"))

    # Set default options based on data characteristics
    _set_defaults(wsp)
    _define_mask(wsp)

    # Do pre-fit on averaged data if requested
    if prefit and max(wsp.asldata.rpts) > 1:
        wsp.log.write(" - Doing initial fit on mean at each TI\n\n")
        init_wsp = wsp.sub("init")
        main_wsp = wsp.sub("main")
        run_multistep_fitting(init_wsp, wsp.asldata.mean_across_repeats())
        wsp.basil_options["continue-from-mvn"] = wsp.init.finalstep.finalMVN
        main_wsp.initmvn = wsp.basil_options["continue-from-mvn"]
    else:
        main_wsp = wsp

    # Then do main run on full ASL data
    wsp.log.write("\n - Doing fit on full ASL data\n\n")
    run_multistep_fitting(main_wsp, wsp.asldata)
    wsp.finalstep = main_wsp.finalstep

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

    if wsp.t1 is None:
        wsp.t1 = t1_default
    if wsp.t1b is None:
        wsp.t1b = 1.65
    if wsp.bat is None:
        wsp.bat = bat_default
    if wsp.batsd is None:
        wsp.batsd = batsd_default
    if wsp.infertiss is None:
        wsp.infertiss = True

    # if we are doing CASL then fix the bolus duration, unless explicitly told us otherwise
    if wsp.infertau is None:
        wsp.infertau = not wsp.asldata.casl

    # Take a copy of extra user-specified BASIL options
    wsp.basil_options = dict(wsp.ifnone("basil_options", {}))

    # If we only have one volume, set a nominal noise prior as it is not possible to
    # estimate from the data
    if wsp.asldata.nvols / wsp.asldata.ntc == 1:
        wsp.log.write(" - Restricting noise prior as only one ASL volume\n")
        wsp.basil_options["prior-noise-stddev"] = 1.0

def _define_mask(wsp):
    # Save a copy of the *analysis* mask to be used for this data as this might change later
    # e.g. for PVC. Fitting may occur in a different mask. Two locations for compatibility...
    if wsp.rois is not None and wsp.rois.mask is not None:
        wsp.analysis_mask = wsp.rois.mask
    else:
        wsp.analysis_mask = wsp.mask

    # Now determine the *fitting* mask
    mask_policy = wsp.ifnone("basil_mask", "default")
    if mask_policy in ("default", "dilated"):
        wsp.log.write(" - Using pipeline analysis mask\n")
        # Two possible locations for compatibility
        if wsp.rois is not None and wsp.rois.mask is not None:
            wsp.basil_mask = wsp.rois.mask
        else:
            wsp.basil_mask = wsp.mask
        if mask_policy == "dilated":
            # Use 3x3x3 kernel for compatibility with fslmaths default
            wsp.log.write(" - Dilating mask for Basil analysis\n")
            struct = scipy.ndimage.generate_binary_structure(3, 3)
            wsp.basil_mask = Image(scipy.ndimage.binary_dilation(wsp.basil_mask.data, structure=struct).astype(np.int32), header=wsp.basil_mask.header)
    elif mask_policy == "none":
            wsp.log.write(" - Not using mask for Basil - will fit every voxel\n")
            wsp.basil_mask = Image(np.ones(wsp.asldata.data.shape[:3]), header=wsp.asldata.header)
    else:
        raise ValueError("Unrecognized mask policy: %s" % mask_policy)

def run_multistep_fitting(wsp, asldata):
    """
    Run multi-step Bayesian model fitting on ASL data

    See ``run`` for details of workspace attributes used

    :param wsp: Workspace object
    :param asldata: AslImage object to use as input data
    """
    if len(asldata.tes) > 1:
        steps = fitting_steps_multite(wsp, asldata)
    else:
        steps = fitting_steps(wsp, asldata)

    prev_result = None
    wsp.asldata_diff = asldata.diff().reorder("rt")

    for idx, step in enumerate(steps):
        step_wsp = wsp.sub("step%i" % (idx+1))
        desc = "Step %i of %i: %s" % (idx+1, len(steps), step.desc)
        if prev_result is not None:
            desc += " - Initialise with step %i" % idx
        step_wsp.log.write(desc + "     ")
        result = step.run(prev_result, log=wsp.log, fsllog=wsp.fsllog, savedir=step_wsp.savedir)
        for key, value in result.items():
            if key == "modelfit":
                # Treat model fit specially - make it an AslImage and also output a mean
                # across repeats version for comparison
                value = wsp.asldata_diff.derived(value.data, header=value.header)
                modelfit_mean = value.mean_across_repeats()
                setattr(step_wsp, "modelfit_mean", modelfit_mean)
            setattr(step_wsp, key, value)

        if step_wsp.logfile is not None and step_wsp.savedir is not None:
            step_wsp.set_item("logfile", step_wsp.logfile, save_fn=str)

        prev_result = result
    wsp.finalstep = step_wsp
    wsp.log.write("\nEnd\n")

def _handle_slicedt(wsp, options):
    """
    Handle the slicedt parameter if we are quantifying in
    a space other than the native ASL space

    Since slicedt is defined relative to the acquisition (ASL) space we
    generate a 3D volume in that space where every voxel contains the
    effective slice time offset. We then transform this into the 
    quantification space we are using and pass it to model fitting using the 'tiimg'
    option.
    """
    img_space = wsp.ifnone("image_space", "asl")
    if img_space != "asl":
        asldata = options["data"]
        _x, _y, z, _t = np.indices(list(asldata.data.shape[:3]) + [asldata.ntis,])
        tis_arr = np.array(asldata.tis) + (z.astype(np.float32) * options["slicedt"])
        tis_img = Image(tis_arr, header=options["data"].header)
        wsp.tiimg = reg.change_space(wsp, tis_img, wsp.ifnone("image_space", "asl"))

        del options["slicedt"]
        ti_idx = 1
        while "ti%i" % ti_idx in options:
            del options["ti%i" % ti_idx]
            ti_idx += 1

        options["tiimg"] = wsp.tiimg

def fitting_steps(wsp, asldata):
    """
    Get the steps required for a BASIL run

    This is separated for the case where an alternative process wants to run
    the actual modelling, or so that the steps can be checked prior to doing
    an actual run.

    Arguments are the same as the ``basil`` function. No workspace is required.
    """
    if asldata is None:
        raise ValueError("Input ASL data is None")

    wsp.log.write("BASIL v%s\n" % __version__)
    asldata.summary(log=wsp.log)
    asldata = asldata.diff().reorder("rt")

    # Default model fitting options for VB runs and spatial steps. Note that values
    # could are None (e.g. sliceband) - these will not be passed to the model fitting step
    options = {
        "data" : asldata,
        "model" : "aslrest",
        "disp" : "none",
        "exch" : "mix",
        "method" : "vb",
        "noise" : "white",
        "allow-bad-voxels" : True,
        "max-iterations" : 20,
        "convergence" : "trialmode",
        "max-trials" : 10,
        "save-mean" : True,
        "save-mvn" : True,
        "save-std" : True,
        "save-noise-mean" : True,
        "save-model-fit" : True,
        "save-residuals" : wsp.ifnone("output_residuals", False),
    }

    if wsp.basil_mask is not None:
        options["mask"] = wsp.basil_mask

    # We choose to pass TIs (not PLDs). The asldata object ensures that
    # TIs are correctly derived from PLDs, when these are specified, by adding
    # the bolus duration.
    options["tis"] = list(asldata.tis)
    options["rpts"] = list(asldata.rpts)

    # Bolus duration - use a single value where possible as cannot infer otherwise
    taus = getattr(asldata, "taus", [1.8,])
    if min(taus) == max(taus):
        options["tau"] = taus[0]
    else:
        options["taus"] = list(taus)

    # Other asl data parameters
    for attr in ("casl", "slicedt", "sliceband"):
        if getattr(asldata, attr, None) is not None:
            options[attr] = getattr(asldata, attr)

    _handle_slicedt(wsp, options)

    if wsp.noiseprior:
        # Use an informative noise prior
        if wsp.noisesd is None:
            snr = wsp.ifnone("snr", 10)
            wsp.log.write(" - Using SNR of %f to set noise std dev\n" % snr)

            # Estimate signal magntiude FIXME diffdata_mean is always 3D?
            if wsp.diffdata_mean.ndim > 3:
                datamax = np.amax(wsp.diffdata_mean.data, 3)
            else:
                datamax = wsp.diffdata_mean.data
            brain_mag = np.mean(datamax.data[wsp.basil_mask.data != 0])
            # this will correspond to whole brain CBF (roughly) - about 0.5 of GM
            noisesd = math.sqrt(brain_mag * 2 / snr)
        else:
            noisesd = wsp.noisesd
        wsp.log.write(" - Using a prior noise sd of: %f\n" % noisesd)
        options["prior-noise-stddev"] = noisesd

    # Add Basil-specific options defined on the workspace
    options.update(wsp.ifnone("basil_options", {}))

    # Additional optional workspace arguments
    for attr in ("t1", "t1b", "bat", "FA", "pwm", "pgm", "batsd"):
        value = getattr(wsp, attr)
        if value is not None:
            options[attr] = value

    # Options for final spatial step
    prior_type_spatial = "M"
    prior_type_mvs = "A"
    options_svb = {
        "method" : "spatialvb",
        "param-spatial-priors" : "N+",
        "convergence" : "maxits",
        "max-iterations": 20,
    }

    wsp.log.write("Model: %s\n" % options["model"])
    wsp.log.write("Dispersion model: %s\n" % options["disp"])
    wsp.log.write("Compartment exchange model: %s\n" % options["exch"])
    inferdisp = options["disp"] != "none"
    inferexch = options["exch"] != "mix"

    # Partial volume correction
    pvcorr = "pgm" in options or "pwm" in options
    if pvcorr:
        if not wsp.infertiss:
            raise ValueError("ERROR: PV correction is not compatible with --artonly option (there is no tissue component)")

        options["incpve"] = True
        if "pgm" not in options or "pwm" not in options:
            raise ValueError("Only one partial volume map (GM / WM) was supplied for PV correctioN")
        # Need a spatial step with more iterations for the PV correction
        wsp.spatial = True
        options_svb["max-iterations"] = 200
        # Ignore partial volumes below 0.1
        pgm_img = options.pop("pgm")
        pwm_img = options.pop("pwm")
        pgm = np.copy(pgm_img.data)
        pwm = np.copy(pwm_img.data)
        pgm[pgm < 0.1] = 0
        pgm[pgm > 1] = 1
        pwm[pwm < 0.1] = 0
        pwm[pwm > 1] = 1
        pgm = Image(pgm, header=pgm_img.header)
        pwm = Image(pwm, header=pwm_img.header)

    # Set general parameter inference and inclusion
    if wsp.infertiss:
        options["inctiss"] = True
    if wsp.inferbat:
        options["incbat"] = True
        options["inferbat"] = True # Infer in first step
    if wsp.inferart:
        options["incart"] = True
    if wsp.inferpc:
        options["incpc"] = True
    if wsp.infertau:
        options["inctau"] = True
    if wsp.infert1:
        options["inct1"] = True

    # Keep track of the number of spatial priors specified by name
    spriors = 1

    if wsp.initmvn:
        # we are being supplied with an initial MVN
        wsp.log.write("Initial MVN being loaded %s\n" % wsp.initmvn.name)
        options["continue-from-mvn"] = wsp.initmvn

    # T1 image prior
    if wsp.t1im is not None:
        spriors = _add_prior(options, spriors, "T_1", type="I", image=wsp.t1im)

    # BAT image prior
    if wsp.batim is not None:
        # With a BAT image prior we must include BAT even if we are not inferring it
        # (in this case the image prior will be treated as ground truth)
        spriors = _add_prior(options, spriors, "delttiss", type="I", image=wsp.batim)
        options["incbat"] = True

    steps = []
    components = ""

    ### --- TISSUE MODULE ---
    if wsp.infertiss:
        components += " Tissue "
        options["infertiss"] = True
        step_desc = "VB - %s" % components
        if not wsp.onestep:
            steps.append(FittingStep(wsp, options, step_desc))

        # setup spatial priors ready
        spriors = _add_prior(options_svb, spriors, "ftiss", type=prior_type_spatial)

    ### --- ARTERIAL MODULE ---
    if wsp.inferart:
        components += " Arterial "
        options["inferart"] = True
        step_desc = "VB - %s" % components
        if not wsp.onestep:
            steps.append(FittingStep(wsp, options, step_desc))

        # setup spatial priors ready
        spriors = _add_prior(options_svb, spriors, "fblood", type=prior_type_mvs)

    ### --- BOLUS DURATION MODULE ---
    if wsp.infertau:
        components += " Bolus duration "
        options["infertau"] = True
        step_desc = "VB - %s" % components
        if not wsp.onestep:
            steps.append(FittingStep(wsp, options, step_desc))

    ### --- MODEL EXTENSIONS MODULE ---
    # Add variable dispersion and/or exchange parameters and/or pre-capiliary
    if inferdisp or inferexch or wsp.inferpc:
        if inferdisp:
            components += " dispersion"
            options["inferdisp"] = True
        if inferexch:
            components += " exchange"
            options["inferexch"] = True
        if wsp.inferpc:
            components += " pre-capiliary"
            options["inferpc"] = True

        step_desc = "VB - %s" % components
        if not wsp.onestep:
            steps.append(FittingStep(wsp, options, step_desc))

    ### --- T1 MODULE ---
    if wsp.infert1:
        components += " T1 "
        options["infert1"] = True
        step_desc = "VB - %s" % components
        if not wsp.onestep:
            steps.append(FittingStep(wsp, options, step_desc))

    ### --- PV CORRECTION MODULE ---
    if pvcorr:
        # setup ready for PV correction, which has to be done with spatial priors
        components += " PVE"
        options["pvcorr"] = True

        # set the image priors for the PV maps
        spriors = _add_prior(options, spriors, "pvgm", type="I", image=pgm)
        spriors = _add_prior(options, spriors, "pvwm", type="I", image=pwm)
        spriors = _add_prior(options, spriors, "fwm", type="M")

        if steps:
            # Add initialisaiton step for PV correction - ONLY if we have something to init from
            steps.append(PvcInitStep(wsp, {"data" : asldata, "mask" : wsp.basil_mask, "pgm" : pgm, "pwm" : pwm}, "PVC initialisation"))

    ### --- SPATIAL MODULE ---
    if wsp.spatial:
        step_desc = "Spatial VB - %s" % components
        options.update(options_svb)
        del options["max-trials"]

        if not wsp.onestep:
            steps.append(FittingStep(wsp, options, step_desc))

    ### --- SINGLE-STEP OPTION ---
    if wsp.onestep:
        steps.append(FittingStep(wsp, options, step_desc))

    if not steps:
        raise ValueError("No steps were generated - no parameters were set to be inferred")

    return steps

def fitting_steps_multite(wsp, asldata, **kwargs):
    """
    Get the steps required for a BASIL run on multi-TE data

    This is separated for the case where an alternative process wants to run
    the actual modelling, or so that the steps can be checked prior to doing
    an actual run.

    Arguments are the same as the ``fitting_steps`` function.
    """
    if asldata is None:
        raise ValueError("Input ASL data is None")

    wsp.log.write("BASIL v%s\n" % __version__)
    asldata.summary(log=wsp.log)
    asldata = asldata.diff().reorder("rt")

    # Default model fitting options for VB runs and spatial steps.
    options = {
        "data" : asldata,
        "model" : "asl_multite",
        "method" : "vb",
        "noise" : "white",
        "allow-bad-voxels" : True,
        "max-iterations" : 20,
        "convergence" : "trialmode",
        "max-trials" : 10,
        "save-mean" : True,
        "save-mvn" : True,
        "save-std" : True,
        "save-model-fit" : True,
        "mask" : wsp.basil_mask,
    }

    # We choose to pass TIs (not PLDs). The asldata object ensures that
    # TIs are correctly derived from PLDs, when these are specified, by adding
    # the bolus duration.
    _list_option(options, asldata.tis, "ti")

    # Pass multiple TEs
    _list_option(options, asldata.tes, "te")

    # Bolus duration must be constant for multi-TE model
    if min(asldata.taus) != max(asldata.taus):
        raise ValueError("Multi-TE model does not support variable bolus durations")
    else:
        options["tau"] = asldata.taus[0]

    # Repeats must be constant for multi-TE model
    if min(asldata.rpts) != max(asldata.rpts):
        raise ValueError("Multi-TE model does not support variable repeats")
    else:
        options["repeats"] = asldata.rpts[0]

    # Other asl data parameters
    for attr in ("casl", "slicedt", "sliceband"):
        if getattr(asldata, attr, None) is not None:
            options[attr] = getattr(asldata, attr)

    # Keyword arguments override options
    options.update(kwargs)

    # Additional optional workspace arguments
    for attr in ("t1", "t1b", "t2", "t2b"):
        value = getattr(wsp, attr)
        if value is not None:
            if attr.startswith("t2"):
                # Model expects T2 in seconds not ms
                options[attr] = float(value) / 1000
            else:
                options[attr] = value

    # Options for final spatial step
    prior_type_spatial = "M"
    prior_type_mvs = "A"
    options_svb = {
        "method" : "spatialvb",
        "param-spatial-priors" : "N+",
        "convergence" : "maxits",
        "max-iterations": 20,
    }

    wsp.log.write("Model: %s\n" % options["model"])
    
    # Set general parameter inference and inclusion
    if not wsp.infertiss:
        wsp.log.write("WARNING: infertiss=False but ftiss is always inferred in multi-TE model\n")
    if not wsp.inferbat:
        wsp.log.write("WARNING: inferbat=False but BAT is always inferred in multi-TE model\n")
    if wsp.inferart:
        wsp.log.write("WARNING: inferart=True but multi-TE model does not support arterial component\n")
    if wsp.infertau:
        options["infertau"] = True
    if wsp.infert1:
        options["infert1"] = True
    if wsp.infert2:
        options["infert2"] = True

    # Keep track of the number of spatial priors specified by name
    spriors = 1

    if wsp.initmvn:
        # we are being supplied with an initial MVN
        wsp.log.write("Initial MVN being loaded %s\n" % wsp.initmvn.name)
        options["continue-from-mvn"] = wsp.initmvn

    # T1 image prior
    if wsp.t1im:
        spriors = _add_prior(options, spriors, "T_1", type="I", image=wsp.t1im)

    # BAT image prior
    if wsp.batim is not None:
        # With a BAT image prior we must include BAT even if we are not inferring it
        # (in this case the image prior will be treated as ground truth)
        spriors = _add_prior(options, spriors, "delttiss", type="I", image=wsp.batim)
        options["incbat"] = True

    steps = []
    components = ""

    ### --- TISSUE MODULE ---
    #if wsp.infertiss:
    if True:
        components += " Tissue"

        ### Inference options
        if wsp.infertau:
            components += " Bolus duration"
            options["infertau"] = True
        if wsp.infert1:
            components += " T1"
            options["infert1"] = True
        if wsp.infertexch:
            components += " Exchange time"
            options["infertexch"] = True

        step_desc = "VB - %s" % components
        if not wsp.onestep:
            steps.append(FittingStep(wsp, options, step_desc))

        # Setup spatial priors ready
        spriors = _add_prior(options_svb, spriors, "ftiss", type=prior_type_spatial)

    ### --- SPATIAL MODULE ---
    if wsp.spatial:
        step_desc = "Spatial VB - %s" % components
        options.update(options_svb)
        del options["max-trials"]

        if not wsp.onestep:
            steps.append(FittingStep(wsp, options, step_desc))

    ### --- SINGLE-STEP OPTION ---
    if wsp.onestep:
        steps.append(FittingStep(wsp, options, step_desc))

    if not steps:
        raise ValueError("No steps were generated - no parameters were set to be inferred")

    return steps

def _list_option(options, values, name):
    for idx, value in enumerate(values):
        options["%s%i" % (name, idx+1)] = value

def _add_prior(options, prior_idx, param, **kwargs):
    options["PSP_byname%i" % prior_idx] = param
    for key, value in kwargs.items():
        options["PSP_byname%i_%s" % (prior_idx, key)] = value
    return prior_idx + 1

class Step(object):
    """
    A step in the fitting process
    """
    def __init__(self, wsp, options, desc):
        self.wsp = wsp
        self.options = dict(options)
        self.desc = desc
        # Need to convert all images to quantification image space
        for key in list(options.keys()):
            poss_img = self.options[key]
            if isinstance(poss_img, Image):
                image_space = wsp.ifnone("image_space", "asl")
                self.options[key] = reg.change_space(wsp, poss_img, image_space, mask=(key == 'mask'))

class FittingStep(Step):
    """
    A step which fits a model to the data
    """
    def __init__(self, wsp, options, desc):
        Step.__init__(self, wsp, options, desc)
        self._impl = wsp.ifnone("basil_method", "fabber")

    def run(self, prev_output, log=sys.stdout, fsllog=None, savedir=None, **kwargs):
        """
        Run model fitting (fabber or vaby), initialising it from the output of a previous step
        """
        print("steprun: ", prev_output)
        if prev_output is not None:
            self.options["continue-from-mvn"] = prev_output["finalMVN"]
        if self._impl == "fabber":
            from oxasl.wrappers import fabber
            fabber_options = dict(self.options)
            for k, v in self.options.items():
                if isinstance(v, list):
                    fabber_options.pop(k)
                    _list_option(fabber_options, v, k[:-1])
            ret = fabber(fabber_options, output=LOAD, progress_log=log, log=fsllog, **kwargs)
        elif self._impl == "vaby":
            print("vaby")
            import vaby
            vaby_options = {}
            for k, v in self.options.items():
                k = k.replace("-", "_")
                if isinstance(v, Image):
                    v = v.nibImage
                vaby_options[k] = v
            vaby_options["model_name"] = vaby_options.pop("model")
            vaby_options["initial_posterior"] = vaby_options.pop("continue-from-mvn", None)
            vaby_options["repeats"] = vaby_options.pop("rpts", None)
            vaby_options["method"] = "avb"
            vaby_options["save_posterior"] = True
            vaby_options["output"] = savedir
            vaby_options["save_log"] = True
            output = {}
            print("running vaby")
            _runtime, _state = vaby.run(outdict=output, **vaby_options)
            print(_state)
            ret = {}
            for k, v in output.items():
                v = Image(v)
                if k == "posterior":
                    ret["finalMVN"] = v
                else:
                    ret[k.replace("_native", "")] = v
        else:
            raise ValueError("Unrecognized model fitting implementation: %s", self._impl)
        log.write("\n")
        return ret

class PvcInitStep(Step):
    """
    A step which initialises partial volume correction
    """
    def run(self, prev_output, log=sys.stdout, fsllog=None, **kwargs):
        """
        Update the MVN from a previous step to include initial estimates
        for PVC parameters
        """
        log.write("Initialising partial volume correction...\n")
        # set the inital GM amd WM values using a simple PV correction
        wm_cbf_ratio = 0.4

        # Modified pvgm map
        temp_pgm = np.copy(self.options["pgm"].data)
        temp_pgm[temp_pgm < 0.2] = 0.2

        # First part of correction psuedo WM CBF term
        prev_ftiss = prev_output["mean_ftiss"].data
        wm_cbf_term = (prev_ftiss * wm_cbf_ratio) * self.options["pwm"].data

        gmcbf_init = (prev_ftiss - wm_cbf_term) / temp_pgm
        wmcbf_init = gmcbf_init * wm_cbf_ratio

        mvn = prev_output["finalMVN"]
        gmcbf_init = Image(gmcbf_init, header=mvn.header)
        wmcbf_init = Image(wmcbf_init, header=mvn.header)

        # HACK: This seems to be required to get the fslpy decorators to write
        # the temporary file correctly
        mask = Image(self.options["mask"].data, header=self.options["mask"].header)

        # load these into the MVN
        mvn = prev_output["finalMVN"]
        from oxasl.wrappers import mvntool
        params = prev_output["paramnames"]
        mvn = mvntool(mvn, params.index("ftiss")+1, output=LOAD, mask=mask, write=True, valim=gmcbf_init, var=0.1, log=fsllog)["output"]
        mvn = mvntool(mvn, params.index("fwm")+1, output=LOAD, mask=mask, write=True, valim=wmcbf_init, var=0.1, log=fsllog)["output"]
        log.write("DONE\n")
        return {"finalMVN" : mvn, "gmcbf_init" : gmcbf_init, "wmcbf_init" : wmcbf_init}
