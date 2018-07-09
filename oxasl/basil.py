#!/usr/bin/env python
"""
BASIL - Bayesian model fitting for ASL

Copyright (c) 2008-2018 University of Oxford
"""

import sys

from fsl.wrappers import LOAD

from ._version import __version__, __timestamp__
from .image import AslImage, AslImageOptions
from .workspace import Workspace
from .options import AslOptionParser, OptionCategory, IgnorableOptionGroup, GenericOptions
from .wrappers import fabber, mvntool

def basil(asldata, wsp=None, log=sys.stdout, **kwargs):
    """
    Run Bayesian model fitting on ASL data

    :param asldata: ASL data (can be AslImage or fsl.Image)
    :param wsp: Optional Workspace object for storing output files. If not specified
                and in-memory workspace is created and returned
    :param log: File-like stream for logging output. Not required if Workspace
                is provided as in this case the workspace log is used
    :param mask: Brain mask (fsl.Image)
    :param infertiss: If True, infer tissue component
    :param inferbat: If True, infer bolus arrival time
    :param infertau: If True, infer bolus duration
    :param inferart: If True, infer arterial component
    :param infert1: If True, infer T1 
    :param inferpc: If True, infer PC 
    :param spatial: If True, include final spatial VB step
    :param onestep: If True, do all inference in a single step
    :param t1im: T1 map as Image
    :param pgm:  Grey matter partial volume map as Image
    :param pwm: White matter partial volume map as Image
    :param initmvn: MVN structure to use as initialization as Image   
    :param kwargs: Additional model options can be passed as keyword arguments,
                   e.g. ``ti1=1.8``
    :return: Workspace containing output. This is an object whose attributes
             are the output data from the run. If a workspace was provided which
             is associated with a physical directory, output data is saved 
             there too.
    """   
    if not wsp:
        wsp = Workspace(log=log)

    steps = get_steps(asldata, log=wsp.log, **kwargs)
    prev_result = None
    for idx, step in enumerate(steps):
        step_wsp = wsp.sub("step%i" % (idx+1))
        desc = "Step %i of %i: %s" % (idx+1, len(steps), step.desc)
        if prev_result:
            desc += " - Initialise with step %i" % idx
        step_wsp.log.write(desc + "\n")
        result = step.run(prev_result, log=log)
        for key, value in result.items():
            setattr(step_wsp, key, value)
        prev_result = result
    wsp.log.write("\nEnd\n")
    return wsp

def get_steps(asldata, mask=None, 
              infertiss=True, inferbat=True, 
              infertau=False, inferart=False, infert1=False, inferpc=False,
              spatial=False, onestep=False,
              t1im=None, pgm=None, pwm=None,
              initmvn=None,
              log=sys.stdout, **kwargs):
    """
    Get the steps required for a BASIL run

    This is separated for the case where an alternative process wants to run
    the actual modelling, or so that the steps can be checked prior to doing
    an actual run.

    Arguments are the same as the ``basil`` function. No workspace is required.
    """
    if not asldata:
        raise ValueError("Input ASL data is None")

    log.write("BASIL v%s\n" % __version__)
    asldata.summary(log=log)
    asldata = asldata.diff()

    # Spatial prior types
    prior_type_spatial = kwargs.pop("spatial_prior", "M")
    prior_type_mvs = kwargs.pop("ard_prior", "A") 

    # Fabber command line options for all runs, VB runs and spatial runs
    options = {
        "model" : "aslrest",
        "method" : "vb",
        "noise" : "white",
        "allow-bad-voxels" : True,
        "max-iterations" : 20,
        "convergence" : "trialmode",
        "max-trials" : 10,
        "disp" : "none",
        "exch" : "mix",
        "save-mean" : True,
        "save-mvn" : True,
        "save-std" : True,
        "data" : asldata,
        "mask" : mask,
    }
    for idx, ti in enumerate(asldata.tis):
        options["ti%i" % (idx+1)] = ti
        options["rpt%i" % (idx+1)] = asldata.rpts[idx]

    # Options to override for final spatial step
    options_svb = {
        "method" : "spatialvb",
        "param-spatial-priors" : "N+",
        "max-iterations": kwargs.pop("svb-max-iterations", 20),
    }

    # Options can be overridden using keyword arguments
    options.update(kwargs)

    log.write("Model (in fabber) is : %s\n" % options["model"])
    log.write("Dispersion model option is %s\n" % options["disp"])
    log.write("Compartment exchange model option is %s\n" % options["exch"])
    inferdisp = options["disp"] != "none"
    inferexch = options["exch"] != "mix"

    # Partial volume correction
    pvcorr = pgm or pwm
    if pvcorr:
        if pgm is None or pwm is None:
            raise ValueError("Only one partial volume map (GM / WM) was supplied for PV correctioN")
        # Need a spatial step with more iterations for the PV correction
        spatial = True
        options_svb["max-iterations"] = 200
        
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

    if initmvn:
        # we are being supplied with an initial MVN
        log.write("Initial MVN being loaded %s\n" % initmvn.name)
        options["continue-from-mvn"] = initmvn
    
    # T1 image prior
    if t1im:
        spriors = _add_prior(options, spriors, "T_1", type="I", image=t1im)

    steps = []
    components = ""

    ### --- TISSUE MODULE ---
    if infertiss:
        components += " Tissue "
        options["infertiss"] = True
        step_desc = "VB - %s" % components
        if not onestep:
            steps.append(FabberStep(options, step_desc))

        # setup spatial priors ready
        spriors = _add_prior(options_svb, spriors, "ftiss", type=prior_type_spatial)

    ### --- ARTERIAL MODULE ---
    if inferart:
        components += " Arterial "
        options["inferart"] = True
        step_desc = "VB - %s" % components
        if not onestep:
            steps.append(FabberStep(options, step_desc))

        # setup spatial priors ready
        spriors = _add_prior(options_svb, spriors, "fblood", type=prior_type_mvs)

    ### --- BOLUS DURATION MODULE ---
    if infertau:
        components += " Bolus duration "
        options["infertau"] = True
        step_desc = "VB - %s" % components
        if not onestep:
            steps.append(FabberStep(options, step_desc))

    ### --- MODEL EXTENSIONS MODULE ---
    # Add variable dispersion and/or exchange parameters and/or pre-capiliary
    if inferdisp or inferexch or inferpc:
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
        if not onestep:
            steps.append(FabberStep(options, step_desc))

    ### --- T1 MODULE ---
    if infert1:
        components += " T1 "
        options["infert1"] = True
        step_desc = "VB - %s" % components
        if not onestep:
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
            steps.append(PvcInitStep({"data" : asldata, "mask" : mask, "pgm" : pgm, "pwm" : pwm}, "PVC initialisation"))

    ### --- SPATIAL MODULE ---
    if spatial:
        step_desc = "Spatial VB - %s" % components
        options.update(options_svb)
        del options["max-trials"]

        if not onestep:
            steps.append(FabberStep(options, step_desc))

    ### --- SINGLE-STEP OPTION ---
    if onestep:
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

        return fabber(self.options, output=LOAD, progress=log)

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
        temp_pgm = self.options["pgm"].data
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

        # load these into the MVN, GM cbf is always param 1
        mvn = prev_output["finalMVN"]
        mvn = mvntool(mvn, "ftiss", output=LOAD, mask=mask, param_list="FIXME", write=True, valim=gmcbf_init, var=0.1)
        mvn = mvntool(mvn, "fwm", output=LOAD, mask=mask, param_list="FIXME", write=True, valim=wmcbf_init, var=0.1)
        log.write("DONE\n")
        return {"finalMVN" : mvn}

class BasilOptions(OptionCategory):
    """
    BASIL option category
    """

    def __init__(self, **kwargs):
        OptionCategory.__init__(self, "basil", **kwargs)

    def groups(self, parser):
        groups = []
        
        group = IgnorableOptionGroup(parser, "BASIL options", ignore=self.ignore)
        group.add_option("--optfile", "-@", dest="optfile", help="If specified, file containing additional Fabber options (e.g. --ti1=1.8)")
        groups.append(group)

        group = IgnorableOptionGroup(parser, "Extended options", ignore=self.ignore)
        group.add_option("--infertau", dest="infertau", help="Infer bolus duration", action="store_true", default=False)
        group.add_option("--inferart", dest="inferart", help="Infer macro vascular (arterial) signal component", action="store_true", default=False)
        group.add_option("--inferpc", dest="inferpc", help="Infer pre-capillary signal component", action="store_true", default=False)
        group.add_option("--infert1", dest="infert1", help="Include uncertainty in T1 values", action="store_true", default=False)
        group.add_option("--artonly", dest="artonly", help="Remove tissue component and infer only arterial component", action="store_true", default=False)
        group.add_option("--fixbat", dest="fixbat", help="Fix bolus duration", action="store_true", default=False)
        group.add_option("--spatial", dest="spatial", help="Add step that implements adaptive spatial smoothing on CBF", action="store_true", default=False)
        group.add_option("--fast", dest="fast", help="Faster analysis (1=faster, 2=single step", type=int, default=0)
        groups.append(group)

        group = IgnorableOptionGroup(parser, "Model options", ignore=self.ignore)
        group.add_option("--disp", dest="disp", help="Model for label dispersion", default="none")
        group.add_option("--exch", dest="exch", help="Model for tissue exchange (residue function)", default="mix")
        groups.append(group)

        group = IgnorableOptionGroup(parser, "Partial volume correction / CBF estimation (enforces --spatial)", ignore=self.ignore)
        group.add_option("--pgm", dest="pgm", help="Gray matter PV map", type="image")
        group.add_option("--pwm", dest="pwm", help="White matter PV map", type="image")
        groups.append(group)

        group = IgnorableOptionGroup(parser, "Special options", ignore=self.ignore)
        group.add_option("--t1im", dest="t1im", help="Voxelwise T1 tissue estimates", type="image")
        groups.append(group)

        return groups

def main():
    """
    Entry point for BASIL command line application
    """
    try:
        parser = AslOptionParser(usage="basil -i <ASL input file> [options...]", version=__version__)
        parser.add_category(AslImageOptions())
        parser.add_category(BasilOptions())
        parser.add_category(GenericOptions())
        
        options, _ = parser.parse_args(sys.argv)
        if not options.output:
            options.output = "basil"

        if not options.asldata:
            sys.stderr.write("Input file not specified\n")
            parser.print_help()
            sys.exit(1)

        options.asldata = AslImage(options.asldata, **parser.filter(vars(options), "image"))

        # Deal with --artonly
        if options.artonly:
            options.infertiss = False
            options.inferart = True

        # Deal with --fixbat
        options.inferbat = not options.fixbat

        # Adjust number of iterations based on fast option
        if not options.fast:
            num_iter, num_trials, onestep = 20, 10, False
        elif options.fast == 1:
            num_iter, num_trials, onestep = 5, 2, False
        elif options.fast == 2:
            num_iter, num_trials, onestep = 10, 5, True
        else:
            raise ValueError("Not a valid option for fast: %s" % str(options.fast))
        options.max_iterations = num_iter
        options.max_trials = num_trials
        options.onestep = onestep

        # Read in additional model options from a file
        if options.optfile:
            for line in open(options.optfile):
                keyval = line.strip().rstrip("\n").lstrip("--").split("=", 1)
                key = keyval[0].strip().replace("-", "_")
                if key != "":
                    if len(keyval) == 1:
                        setattr(options, key, True)
                    else:
                        setattr(options, key, keyval[1].strip())

        # Run BASIL processing, passing options as keyword arguments using **
        wsp = Workspace(savedir=options.output)
        basil(wsp=wsp, **vars(options))
        
    except ValueError as exc:
        sys.stderr.write("\nERROR: " + str(exc) + "\n")
        sys.stderr.write("Use --help for usage information\n")
        sys.exit(1)

if __name__ == "__main__":
    main()
