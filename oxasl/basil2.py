#!/usr/bin/env python
"""
BASIL - Bayesian model fitting for ASL

The BASIL module is a little more complex than the other Workspace based
modules because of the number of options available and the need for flexibility
in how the modelling steps are run.

The main function is ``basil`` which performs model fitting on ASL data 
in the Workspace ``asldata`` attribute. 

    wsp = Workspace()
    wsp.asldata = AslImage("asldata.nii.gz", tis=[1.6,])
    wsp.infertiss = True
    basil(wsp, output_wsp=wsp.sub("basil"))
    basil.finalstep.mean_ftiss.save("mean_ftiss.nii.gz")

Because of the number of options possible for the modelling process, the 
workspace attribute ``basil_options`` can be set as a dictionary of extra
options relevant only to Basil:

    wsp = Workspace()
    wsp.asldata = AslImage("asldata.nii.gz", tis=[1.6,])
    wsp.basil_options = {"infertiss" : True, "spatial" : True}
    basil(wsp, output_wsp=wsp.sub("basil"))
    basil.finalstep.mean_ftiss.save("mean_ftiss.nii.gz")

All options specified in basil_options are either consumed by Basil, or
if not passed directly to the model.

Copyright (c) 2008-2018 University of Oxford
"""

import sys
import math

import numpy as np

from fsl.wrappers import LOAD
from fsl.data.image import Image

from oxasl import __version__, __timestamp__, Workspace, image
from oxasl.options import AslOptionParser, OptionCategory, IgnorableOptionGroup, GenericOptions

def basil(wsp, prefit=True):
    wsp.log.write("\nRunning BASIL Bayesian modelling on ASL data\n")

    if wsp.asldata is None:
        raise ValueError("Input ASL data is None")

    wsp.log.write("  - BASIL v%s\n" % __version__)
    wsp.asldata.summary(log=wsp.log)

    if prefit and max(wsp.asldata.rpts) > 1:
        # Initial BASIL run on mean data
        wsp.log.write(" - Doing initial fit on mean at each TI\n\n")
        init_wsp = wsp.sub("init")
        init_wsp.asldata = wsp.asldata.mean_across_repeats()
        runbasil(init_wsp)
        main_wsp = wsp.sub("main")
        main_wsp.initmvn = wsp.init.finalstep.finalMVN
    else:
        main_wsp = wsp

    # Main run on full ASL data
    wsp.log.write("\n - Doing fit on full ASL data\n\n")
    runbasil(main_wsp)

    # Provide a reference to the final step
    wsp.finalstep = main_wsp.finalstep

def runbasil(wsp):
    """
    Run Bayesian model fitting on ASL data

    See ``basil`` for details of workspace attributes used

    :param wsp: Workspace object
    :param asldata: AslImage object to use as input data
    :param output_wsp: Optional Workspace object for storing output files. If not specified
                       ``wsp`` is used instead
    """
    steps = basil_steps(wsp)

    prev_result = None
    wsp.asldata_diff = wsp.asldata.diff().reorder("rt")

    for idx, step in enumerate(steps):
        step_wsp = wsp.sub("step%i" % (idx+1))
        desc = "Step %i of %i: %s" % (idx+1, len(steps), step.desc)
        if prev_result is not None:
            desc += " - Initialise with step %i" % idx
        step_wsp.log.write(desc + "     ")
        result = step.run(prev_result, log=wsp.log)
        for key, value in result.items():
            setattr(step_wsp, key, value)

        if step_wsp.logfile is not None and step_wsp.savedir is not None:
            step_wsp.set_item("logfile", step_wsp.logfile, save_fn=str)

        prev_result = result
    wsp.finalstep = step_wsp
    wsp.log.write("\nEnd\n")

def basil_steps(wsp):
    # Get data in the format expected by BASIL
    asldata = wsp.asldata.diff().reorder("rt")

    # Default Fabber options for VB runs and spatial steps. Note that attributes
    # which are None (e.g. sliceband) are not passed to Fabber
    options = {
        "data" : asldata,
        "mask" : wsp.rois.mask,
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
        "save-model-fit" : True,
        "t1b" : 1.65,
    }

    # Default is to infer tissue parameters only
    infertiss = wsp.ifnone("infertiss", True)
    inferbat = wsp.ifnone("inferbat", False)
    inferart = wsp.ifnone("inferart", False)
    inferpc = wsp.ifnone("inferpc", False)
    infert1 = wsp.ifnone("infert1", False)

    # Default to inferring exch/disp parameters if a relevant model is used
    inferdisp = wsp.ifnone("inferdisp", options["disp"] != "none")
    inferexch = wsp.ifnone("inferexch", options["exch"] != "mix")

    # if we are doing CASL then fix the bolus duration, unless explicitly told us otherwise
    infertau = wsp.ifnone("infertau", not asldata.casl)
    
    # Single or Multi TI setup
    if asldata.ntis == 1:
        # Single TI data - don't try to infer arterial component of bolus duration, we don't have enough info
        wsp.log.write(" - Operating in Single TI mode - no arterial component, fixed bolus duration\n")
        inferart = False
        infertau = False
        options["batsd"] = 0.1
    else:
        # For multi TI/PLD data, set a more liberal prior for tissue ATT since we should be able to 
        # determine this from the data. NB this leaves the arterial BAT alone.
        options["batsd"] = 1

    if wsp.wp:
        # White paper mode - this overrides defaults, but can be overwritten by command line 
        # specification of individual parameters
        wsp.log.write(" - Analysis in white paper mode: T1 default=1.65, BAT default=0, voxelwise calibration\n")
        options["t1"] = 1.65
        options["bat"] = 0.0
    else:
        options["t1"] = 1.3
        if asldata.casl:
            options["bat"] = 1.3
        else:
            options["bat"] = 0.7

    # Override defaults using workspace arguments
    for attr in ("t1", "t1b", "FA", "bat", "batsd", "exch", "disp"):
        value = getattr(wsp, attr)
        if value is not None:
            options[attr] = value

    # We choose to pass TIs (not PLDs). The asldata object ensures that
    # TIs are correctly derived from PLDs, when these are specified, by adding
    # the bolus duration.
    for idx, ti in enumerate(asldata.tis):
        options["ti%i" % (idx+1)] = ti
        options["rpt%i" % (idx+1)] = asldata.rpts[idx]

    # Bolus duration - use a single value where possible as cannot infer otherwise
    taus = getattr(asldata, "taus", [1.8,])
    if min(taus) == max(taus):
        options["tau"] = taus[0]
    else:
        for idx, tau in enumerate(taus):
            options["tau%i" % (idx+1)] = tau

    # Other asl data parameters
    for attr in ("casl", "slicedt", "sliceband"):
        if getattr(asldata, attr, None) is not None:
            options[attr] = getattr(asldata, attr)

    if wsp.noiseprior:
        # Use an informative noise prior
        if wsp.noisesd is None:
            snr = wsp.ifnone("snr", 10)
            wsp.log.write(" - Using SNR of %f to set noise std dev\n" % snr)

            # Estimate signal magntiude FIXME diffdata_mean is always 3D?
            diffdata_mean = asldata.mean()
            if diffdata_mean.ndim > 3:
                datamax = np.amax(diffdata_mean.data, 3)
            else:
                datamax = diffdata_mean.data
            brain_mag = np.mean(datamax.data[wsp.rois.mask.data != 0])
            # this will correspond to whole brain CBF (roughly) - about 0.5 of GM
            noisesd = math.sqrt(brain_mag * 2 / snr)
        else:
            noisesd = wsp.noisesd
        wsp.log.write(" - Using a prior noise sd of: %f\n" % noisesd)
        options["prior-noise-stddev"] = noisesd

    # Options for final spatial step
    prior_type_spatial = "M"
    prior_type_mvs = "A"
    options_svb = {
        "method" : "spatialvb",
        "param-spatial-priors" : "N+",
        "convergence" : "maxits",
        "max-iterations": 20,
    }

    wsp.log.write("Model (in fabber) is : %s\n" % options["model"])
    wsp.log.write("Dispersion model option is %s\n" % options["disp"])
    wsp.log.write("Compartment exchange model option is %s\n" % options["exch"])

    # Partial volume correction
    pvcorr = wsp.pwm is not None or wsp.pgm is not None
    if pvcorr:
        if wsp.pwm is None or wsp.pgm is None:
            raise ValueError("Only one partial volume map (GM / WM) was supplied for PV correctioN")
        # Need a spatial step with more iterations for the PV correction
        options_svb["max-iterations"] = 200
        # Ignore partial volumes below 0.1
        pgm = np.copy(wsp.pgm.data)
        pwm = np.copy(wsp.pwm.data)
        pgm[pgm < 0.1] = 0
        pgm[pgm > 1] = 1
        pwm[pwm < 0.1] = 0
        pwm[pwm > 1] = 1
        pgm = Image(pgm, header=wsp.pgm.header)
        pwm = Image(pwm, header=wsp.pwm.header)
        
    if pvcorr and not infertiss:
        raise ValueError("ERROR: PV correction is not compatible with --artonly option (there is no tissue component)")

    # Set general parameter inference and inclusion
    if infertiss:
        options["inctiss"] = True
    if inferbat:
        options["incbat"] = True
        options["inferbat"] = True # Infer in first step
    if inferart:
        options["incart"] = True
    if inferpc:
        options["incpc"] = True
    if infertau:
        options["inctau"] = True
    if infert1:
        options["inct1"] = True
    if pvcorr:
        options["incpve"] = True

    # Keep track of the number of spatial priors specified by name
    spriors = 1 

    if wsp.initmvn:
        # we are being supplied with an initial MVN
        wsp.log.write("Initial MVN being loaded %s\n" % wsp.initmvn.name)
        options["continue-from-mvn"] = wsp.initmvn
    
    # T1 image prior
    if wsp.t1im:
        spriors = _add_prior(options, spriors, "T_1", type="I", image=wsp.t1im)

    steps = []
    components = ""

    ### --- TISSUE MODULE ---
    if infertiss:
        components += " Tissue "
        options["infertiss"] = True
        step_desc = "VB - %s" % components
        if not wsp.onestep:
            steps.append(FabberStep(options, step_desc))

        # setup spatial priors ready
        spriors = _add_prior(options_svb, spriors, "ftiss", type=prior_type_spatial)

    ### --- ARTERIAL MODULE ---
    if inferart:
        components += " Arterial "
        options["inferart"] = True
        step_desc = "VB - %s" % components
        if not wsp.onestep:
            steps.append(FabberStep(options, step_desc))

        # setup spatial priors ready
        spriors = _add_prior(options_svb, spriors, "fblood", type=prior_type_mvs)

    ### --- BOLUS DURATION MODULE ---
    if infertau:
        components += " Bolus duration "
        options["infertau"] = True
        step_desc = "VB - %s" % components
        if not wsp.onestep:
            steps.append(FabberStep(options, step_desc))

    ### --- MODEL EXTENSIONS MODULE ---
    # Add variable dispersion and/or exchange parameters and/or pre-capiliary
    if inferdisp or inferexch or wsp.inferpc:
        if inferdisp:
            components += " dispersion"
            options["inferdisp"] = True
        if inferexch:
            components += " exchange"
            options["inferexch"] = True
        if inferpc:
            components += " pre-capiliary"
            options["inferpc"] = True

        step_desc = "VB - %s" % components
        if not wsp.onestep:
            steps.append(FabberStep(options, step_desc))

    ### --- T1 MODULE ---
    if infert1:
        components += " T1 "
        options["infert1"] = True
        step_desc = "VB - %s" % components
        if not wsp.onestep:
            steps.append(FabberStep(options, step_desc))

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
            steps.append(PvcInitStep({"data" : asldata, "mask" : wsp.rois.mask, "pgm" : pgm, "pwm" : pwm}, "PVC initialisation"))

    ### --- SPATIAL MODULE ---
    if wsp.spatial:
        step_desc = "Spatial VB - %s" % components
        options.update(options_svb)
        del options["max-trials"]

        if not wsp.onestep:
            steps.append(FabberStep(options, step_desc))

    ### --- SINGLE-STEP OPTION ---
    if wsp.onestep:
        steps.append(FabberStep(options, step_desc))
        
    if not steps:
        raise ValueError("No steps were generated - no parameters were set to be inferred")
        
    return steps

def _add_prior(options, prior_idx, param, **kwargs):
    options["PSP_byname%i" % prior_idx] = param
    for key, value in kwargs.items():
        options["PSP_byname%i_%s" % (prior_idx, key)] = value
    return prior_idx + 1

class Step(object):
    """
    A step in the Basil modelling process
    """
    def __init__(self, options, desc):
        self.options = dict(options)
        self.desc = desc

class FabberStep(Step):
    """
    A Basil step which involves running Fabber
    """
    def run(self, prev_output, log=sys.stdout):
        """
        Run Fabber, initialising it from the output of a previous step
        """
        if prev_output is not None:
            self.options["continue-from-mvn"] = prev_output["finalMVN"]
        from .wrappers import fabber
        ret = fabber(self.options, output=LOAD, progress=log)
        log.write("\n")
        return ret

class PvcInitStep(Step):
    """
    A Basil step which initialises partial volume correction
    """
    def run(self, prev_output, log=sys.stdout):
        """
        Update the MVN from a previous step to include initial estimates
        for PVC parameters
        """
        log.write("Initialising partial volume correction...\n")
        mask = self.options["mask"]
        # set the inital GM amd WM values using a simple PV correction
        wm_cbf_ratio = 0.4

        # Modified pvgm map
        #fsl.maths(pgm, " -sub 0.2 -thr 0 -add 0.2 temp_pgm")
        temp_pgm = np.copy(self.options["pgm"].data)
        temp_pgm[temp_pgm < 0.2] = 0.2

        # First part of correction psuedo WM CBF term
        #fsl.run("mvntool", "--input=temp --output=temp_ftiss --mask=%s --param=ftiss --param-list=step%i/paramnames.txt --val" % (mask.name, prev_step))
        #fsl.maths("temp_ftiss", "-mul %f -mul %s wmcbfterm" % (wm_cbf_ratio, pwm.name))
        prev_ftiss = prev_output["mean_ftiss"].data
        wm_cbf_term = (prev_ftiss * wm_cbf_ratio) * self.options["pwm"].data

        #fsl.maths("temp_ftiss", "-sub wmcbfterm -div temp_pgm gmcbf_init")
        #fsl.maths("gmcbf_init -mul %f wmcbf_init" % wm_cbf_ratio)
        gmcbf_init = (prev_ftiss - wm_cbf_term) / temp_pgm
        wmcbf_init = gmcbf_init * wm_cbf_ratio

        mvn = prev_output["finalMVN"]
        gmcbf_init = Image(gmcbf_init, header=mvn.header)
        wmcbf_init = Image(wmcbf_init, header=mvn.header)

        # load these into the MVN
        mvn = prev_output["finalMVN"]
        from .wrappers import mvntool
        params = prev_output["paramnames"]
        mvn = mvntool(mvn, params.index("ftiss")+1, output=LOAD, mask=mask, write=True, valim=gmcbf_init, var=0.1)["output"]
        mvn = mvntool(mvn, params.index("fwm")+1, output=LOAD, mask=mask, write=True, valim=wmcbf_init, var=0.1)["output"]
        log.write("DONE\n")
        return {"finalMVN" : mvn, "gmcbf_init" : gmcbf_init, "wmcbf_init" : wmcbf_init}

class BasilOptions(OptionCategory):
    """
    BASIL option category
    """

    def __init__(self, **kwargs):
        OptionCategory.__init__(self, "basil", **kwargs)

    def groups(self, parser):
        groups = []
        
        group = IgnorableOptionGroup(parser, "BASIL options", ignore=self.ignore)
        group.add_option("--infertau", help="Infer bolus duration", action="store_true", default=False)
        group.add_option("--inferart", help="Infer macro vascular (arterial) signal component", action="store_true", default=False)
        group.add_option("--inferpc", help="Infer pre-capillary signal component", action="store_true", default=False)
        group.add_option("--infert1", help="Include uncertainty in T1 values", action="store_true", default=False)
        group.add_option("--artonly", help="Remove tissue component and infer only arterial component", action="store_true", default=False)
        group.add_option("--fixbat", dest="inferbat", help="Fix bolus arrival time", action="store_false", default=True)
        group.add_option("--spatial", help="Add step that implements adaptive spatial smoothing on CBF", action="store_true", default=False)
        group.add_option("--fast", help="Faster analysis (1=faster, 2=single step", type=int, default=0)
        group.add_option("--noiseprior", help="Use an informative prior for the noise estimation", action="store_true", default=False)
        group.add_option("--noisesd", help="Set a custom noise std. dev. for the nosie prior", type=float)
        group.add_option("--model-options", help="File containing additional model-specific options")
        groups.append(group)

        group = IgnorableOptionGroup(parser, "Model options", ignore=self.ignore)
        group.add_option("--disp", help="Model for label dispersion", default="none")
        group.add_option("--exch", help="Model for tissue exchange (residue function)", default="mix")
        groups.append(group)

        group = IgnorableOptionGroup(parser, "Partial volume correction / CBF estimation (enforces --spatial)", ignore=self.ignore)
        group.add_option("--pgm", help="Gray matter PV map", type="image")
        group.add_option("--pwm", help="White matter PV map", type="image")
        groups.append(group)

        group = IgnorableOptionGroup(parser, "Special options", ignore=self.ignore)
        group.add_option("--t1im", help="Voxelwise T1 tissue estimates", type="image")
        groups.append(group)

        return groups

def main():
    """
    Entry point for BASIL command line application
    """
    try:
        parser = AslOptionParser(usage="basil -i <ASL input file> [options...]", version=__version__)
        parser.add_category(image.AslImageOptions())
        parser.add_category(BasilOptions())
        parser.add_category(GenericOptions())
        
        options, _ = parser.parse_args()
        if not options.output:
            options.output = "basil"

        if not options.asldata:
            sys.stderr.write("Input file not specified\n")
            parser.print_help()
            sys.exit(1)
        
        wsp = Workspace(savedir=options.output, separate_input=True, auto_asldata=True, **vars(options))
        
        # Deal with --artonly
        if wsp.artonly:
            wsp.infertiss = False
            wsp.inferart = True

        # Adjust number of iterations based on fast option
        if not wsp.fast:
            num_iter, num_trials, onestep = 20, 10, False
        elif wsp.fast == 1:
            num_iter, num_trials, onestep = 5, 2, False
        elif wsp.fast == 2:
            num_iter, num_trials, onestep = 10, 5, True
        else:
            raise ValueError("Not a valid option for fast: %s" % str(wsp.fast))

        wsp.max_iterations = num_iter
        wsp.max_trials = num_trials
        wsp.onestep = onestep

        # Read in additional model options from a file
        wsp.model_options = {}
        if wsp.model_options:
            with open(wsp.model_options) as optfile:
                for line in optfile:
                    keyval = line.strip().rstrip("\n").lstrip("--").split("=", 1)
                    key = keyval[0].strip().replace("-", "_")
                    if key != "":
                        if len(keyval) == 1:
                            wsp.model_options[key] = True
                        else:
                            wsp.model_options[key] = keyval[1].strip()

        # Run BASIL processing
        basil(wsp)
        
    except ValueError as exc:
        sys.stderr.write("\nERROR: " + str(exc) + "\n")
        sys.stderr.write("Use --help for usage information\n")
        sys.exit(1)

if __name__ == "__main__":
    main()
