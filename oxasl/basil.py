#!/usr/bin/env python
"""
BASIL - Bayesian model fitting for ASL

Copyright (c) 2008-2018 University of Oxford
"""

import sys
from optparse import OptionParser

from fsl.data.image import Image

from . import __version__, AslImage, AslOptionGroup, fslwrap as fsl
from .image import add_data_options

__timestamp__ = "TEMP"

def _add_prior(options, prior_idx, param, **kwargs):
    options["PSP_byname%i" % prior_idx] = param
    for key, value in kwargs.items():
        options["PSP_byname%i_%s" % (prior_idx, key)] = value
    return prior_idx + 1

def _add_step(steps, step, step_desc, infile, mask, options, prev_step=None):
    steps.append((step, step_desc, infile, mask, dict(options), prev_step))
    return len(steps)+1, len(steps)

def _do_step(wsp, step, step_desc, infile, mask, options, prev_step=None, log=sys.stdout):
    """
    Perform a single model fitting step

    :param wsp: FSL workspace
    :param step: Step number
    :param step_desc: Description of step. The special value "PVC" indicates partial volume
                      correction preprocessing step
    :param infile: Input data as Image
    :param mask: Optional mask as Image (or None)
    :param options: Dictionary of Fabber options, not including the ``--`` command line prefix
    :param prev_step: Optional number of previous step to initialize from
    :param log: File stream for log output
    """
    if prev_step == "PVC":
        _do_pvc_init(wsp, mask, prev_step, options)
    else:
        if prev_step is not None:
            step_desc += " - init with STEP %i" % prev_step
            options["continue-from-mvn"] = "step%i/finalMVN" % prev_step

        log.write("\n%s\n" % step_desc)
        wsp.fabber("asl", infile, mask, options, output_name="step%i" % step)

def _do_pvc_init(wsp, mask, prev_step, options):
    # set the inital GM amd WM values using a simple PV correction
    wm_cbf_ratio = 0.4

    # Modified pvgm map
    #fsl.maths(pgm, " -sub 0.2 -thr 0 -add 0.2 temp_pgm")
    temp_pgm = options["pgm"].nibImage.get_data()
    temp_pgm[temp_pgm < 0.2] = 0.2

    # First part of correction psuedo WM CBF term
    #fsl.run("mvntool", "--input=temp --output=temp_ftiss --mask=%s --param=ftiss --param-list=step%i/paramnames.txt --val" % (mask.name, prev_step))
    #fsl.maths("temp_ftiss", "-mul %f -mul %s wmcbfterm" % (wm_cbf_ratio, pwm.name))
    prev_ftiss = Image("step%i/mean_ftiss" % prev_step).nibImage.get_data()
    wm_cbf_term = (prev_ftiss * wm_cbf_ratio) * options["pwm"].nibImage.get_data()

    #fsl.maths("temp_ftiss", "-sub wmcbfterm -div temp_pgm gmcbf_init")
    #fsl.maths("gmcbf_init -mul %f wmcbf_init" % wm_cbf_ratio)
    gmcbf_init = (prev_ftiss - wm_cbf_term) / temp_pgm
    wmcbf_init = gmcbf_init * wm_cbf_ratio

    # load these into the MVN, GM cbf is always param 1
    mvnfile = "step%i/finalMVN" % prev_step
    fsl.run("mvntool", "--input=%s --output=%s --mask=%s --param=ftiss --param-list=step%i/paramnames.txt --write --valim=gmcbf_init --var=0.1" % (mask.name, mvnfile, mvnfile, prev_step))
    fsl.run("mvntool", "--input=%s --output=%s --mask=%s --param=fwm --param-list=step%i/paramnames.txt --write --valim=wmcbf_init --var=0.1" % (mask.name, mvnfile, mvnfile, prev_step))

def run_steps(wsp, steps, log=sys.stdout):
    for step in steps:
        _do_step(wsp, *step, log=sys.stdout)
    log.write("End\n")

def get_steps(asldata, mask=None, 
              infertiss=True, inferbat=True, 
              infertau=False, inferart=False, infert1=False, inferpc=False,
              spatial=False, onestep=False,
              t1im=None, pgm=None, pwm=None,
              initmvn=None,
              log=sys.stdout, **kwargs):
    """
    Run Bayesian ASL model fitting

    :param wsp: FSL workspace for output files
    :param asldata: ASL data (can be AslImage or fslwrap.Image)
    :param mask: Brain mask (fslwrap.Image)
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
    :param t1im: White matter partial volume map as Image
    :param initmvn: MVN structure to use as initialization as Image   
    :param kwargs: Additional model options can be passed as keyword arguments,
                   e.g. ``ti1=1.8``
    """
    if not asldata:
        raise ValueError("basil: input ASL data is None")

    log.write("BASIL v%s\n" % __version__)
    #log.write("Working directory: %s\n" % wsp.workdir)
    asldata.summary(log=log)

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
        log.write("Initial MVN being loaded %s" % initmvn.name)
        options["continue-from-mvn"] = initmvn.name
    
    # T1 image prior
    if t1im:
        spriors = _add_prior(options, spriors, "T_1", type="I", image=t1im.name)

    step = 1
    steps = []
    prev_step = None
    step_params = ""

    ### --- TISSUE MODULE ---
    if infertiss:
        step_params += " Tissue "
        options["infertiss"] = True
        step_desc = "STEP %i: VB - %s" % (step, step_params)
        if not onestep:
            step, prev_step = _add_step(steps, step, step_desc, asldata, mask, options, prev_step)

        # setup spatial priors ready
        spriors = _add_prior(options_svb, spriors, "ftiss", type=prior_type_spatial)

    ### --- ARTERIAL MODULE ---
    if inferart:
        step_params += " Arterial "
        options["inferart"] = True
        step_desc = "STEP %i: VB - %s" % (step, step_params)
        if not onestep:
            step, prev_step = _add_step(steps, step, step_desc, asldata, mask, options, prev_step)

        # setup spatial priors ready
        spriors = _add_prior(options_svb, spriors, "fblood", type=prior_type_mvs)

    ### --- BOLUS DURATION MODULE ---
    if infertau:
        step_params += " Bolus duration "
        options["infertau"] = True
        step_desc = "STEP %i: VB - %s" % (step, step_params)
        if not onestep:
            step, prev_step = _add_step(steps, step, step_desc, asldata, mask, options, prev_step)

    ### --- MODEL EXTENSIONS MODULE ---
    # Add variable dispersion and/or exchange parameters and/or pre-capiliary
    if inferdisp or inferexch or inferpc:
        if inferdisp:
            step_params += " dispersion"
            options["inferdisp"] = True
        if inferexch:
            step_params += " exchange"
            options["inferexch"] = True
        if inferpc:
            step_params += " pre-capiliary"
            options["inferpc"] = True

        step_desc = "STEP %i: VB - %s" % (step, step_params)	
        if not onestep:
            step, prev_step = _add_step(steps, step, step_desc, asldata, mask, options, prev_step)

    ### --- T1 MODULE ---
    if infert1:
        step_params += " T1 "
        options["infert1"] = True
        step_desc = "STEP %i: VB - %s" % (step, step_params)
        if not onestep:
            step, prev_step = _add_step(steps, step, step_desc, asldata, mask, options, prev_step)

    ### --- PV CORRECTION MODULE ---
    if pvcorr:
        # setup ready for PV correction, which has to be done with spatial priors
        step_params += " PVE"
        options["pvcorr"] = True

        # set the image priors for the PV maps
        spriors = _add_prior(options, spriors, "pvgm", type="I", image=pgm.name)
        spriors = _add_prior(options, spriors, "pvwm", type="I", image=pwm.name)
        spriors = _add_prior(options, spriors, "fwm", type="M")

        if prev_step:
            # Add initialisaiton step for PV correction - ONLY if we have something to init from 
            # (either step greater than 1 or initmvn set)
            step, prev_step = _add_step(steps, step, "PVC", asldata, mask, {"pgm" : pgm.name, "pwm" : pwm.name}, prev_step)

    ### --- SPATIAL MODULE ---
    if spatial:
        step_desc = "STEP %i: Spatial VB %s" % (step, step_params)
        options.update(options_svb)
        del options["max-trials"]

        if not onestep:
            step, prev_step = _add_step(steps, step, step_desc, asldata, mask, options, prev_step)

    ### --- SINGLE-STEP OPTION ---
    if onestep:
        step, prev_step = _add_step(steps, step, step_desc, asldata, mask, options, prev_step)
        
    if not steps:
        raise ValueError("No steps were generated - no parameters were set to be inferred")
        
    return steps

def add_basil_options(p, ignore=()):
    g = AslOptionGroup(p, "BASIL options", ignore=ignore)
    g.add_option("-m", dest="mask", help="Mask file")
    g.add_option("--optfile", "-@", dest="optfile", help="If specified, file containing additional Fabber options (e.g. --ti1=1.8)")
    p.add_option_group(g)
    g = AslOptionGroup(p, "Extended options", ignore=ignore)
    g.add_option("--infertau", dest="infertau", help="Infer bolus duration", action="store_true", default=False)
    g.add_option("--inferart", dest="inferart", help="Infer macro vascular (arterial) signal component", action="store_true", default=False)
    g.add_option("--inferpc", dest="inferpc", help="Infer pre-capillary signal component", action="store_true", default=False)
    g.add_option("--infert1", dest="infert1", help="Include uncertainty in T1 values", action="store_true", default=False)
    g.add_option("--artonly", dest="artonly", help="Remove tissue component and infer only arterial component", action="store_true", default=False)
    g.add_option("--fixbat", dest="fixbat", help="Fix bolus duration", action="store_true", default=False)
    g.add_option("--spatial", dest="spatial", help="Add step that implements adaptive spatial smoothing on CBF", action="store_true", default=False)
    g.add_option("--fast", dest="fast", help="Faster analysis (1=faster, 2=single step", type=int, default=0)
    p.add_option_group(g)
    g = AslOptionGroup(p, "Model options", ignore=ignore)
    g.add_option("--disp", dest="disp", help="Model for label dispersion", default="none")
    g.add_option("--exch", dest="exch", help="Model for tissue exchange (residue function)", default="mix")
    p.add_option_group(g)
    g = AslOptionGroup(p, "Partial volume correction / CBF estimation (enforces --spatial)", ignore=ignore)
    g.add_option("--pgm", dest="pgm", help="Gray matter PV map")
    g.add_option("--pwm", dest="pwm", help="White matter PV map")
    p.add_option_group(g)
    g = AslOptionGroup(p, "Special options", ignore=ignore)
    g.add_option("--t1im", dest="t1im", help="Voxelwise T1 tissue estimates")
    p.add_option_group(g)

def main():
    """
    Entry point for BASIL command line application
    """
    usage = """BASIL

    basil -i <ASL input file> -m <mask image file> -o <output dir>"""

    try:
        p = OptionParser(usage=usage, version=__version__)
        add_data_options(p)
        add_basil_options(p)
        options, _ = p.parse_args(sys.argv)
        
        # Convert options into a dictionary. This is a bit easier to iterate over
        # and turn into keyword arguments to the basil function
        options = vars(options)

        # Create workspace which is the equivalent of the output directory
        wsp = fsl.Workspace(options["output"], echo=True, debug=options.pop("debug", False))

        # Names and descriptions of options which are images
        images = {
            "mask" : "Brain mask", 
            "t1im": "T1 map", 
            "pwm" : "White matter PV map", 
            "pgm" : "Grey matter PV map",
        }

        # Convert image options into Image objects and copy into workspace - this 
        # also checks they exist and can be loaded
        if not options.get("asldata", None):
            raise ValueError("ASL input data not specified")
        options["asldata"] = AslImage(options["asldata"], **options).diff().reorder("rt")
        wsp.add_img(options["asldata"])
        for opt, role in images.items():
            if options[opt]:
                options[opt] = Image(image=options[opt], name=role)
                wsp.add_img(options[opt])

        # Remove options consumed by AslImage
        for opt in ["order", "ntis", "tis", "nplds", "plds", "nrpts", "rpts", "nphases", "phases"]:
            options.pop(opt, None)
            
        # Deal with --artonly
        artonly = options.pop("artonly")
        if artonly:
            options["infertiss"] = False
            options["inferart"] = True

        # Deal with --fixbat
        options["inferbat"] = not options.pop("fixbat", False)

        # Adjust number of iterations based on fast option
        fast = options.pop("fast", 0)
        if fast == 0:
            num_iter, num_trials, onestep = 20, 10, False
        elif fast == 1:
            num_iter, num_trials, onestep = 5, 2, False
        elif fast == 2:
            num_iter, num_trials, onestep = 10, 5, True
        else:
            raise ValueError("Not a valid option for fast: %s" % str(fast))
        options["max-iterations"] = num_iter
        options["max-trials"] = num_trials
        options["onestep"] = onestep

        # Read in additional model options from a file
        optfile = options.pop("optfile", None)
        if optfile:
            for line in open(optfile):
                keyval = line.strip().rstrip("\n").lstrip("--").split("=", 1)
                key = keyval[0].strip()
                if key != "":
                    if len(keyval) == 1:
                        options[key] = True
                    else:
                        options[key] = keyval[1].strip()

        # Run BASIL processing, passing options as keyword arguments using **
        steps = get_steps(**options)
        run_steps(wsp, steps)
        
    except ValueError as e:
        sys.stderr.write("\nERROR: " + str(e) + "\n")
        sys.stderr.write("Use --help for usage information\n")
        sys.exit(1)

if __name__ == "__main__":
    main()
