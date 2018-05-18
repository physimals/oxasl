#!/bin/env python
"""
OXFORD_ASL: Converts ASL images in to perfusion maps

Michael Chappell, FMRIB Image Analysis Group & IBME

Copyright (c) 2008-2013 Univerisity of Oxford
"""

import sys
import os
import traceback
from optparse import OptionParser

import numpy as np

from . import __version__, AslOptionGroup, AslImage, fslwrap as fsl, basil
from .image import add_data_options
from .calib import add_calib_options
from .basil import add_basil_options

def add_oxasl_options(parser, ignore=()):
    g = AslOptionGroup(parser, "Main Options", ignore=ignore)
    g.add_option("-m", dest="mask", help="Brain mask (in native space of ASL data)")
    g.add_option("--wp", dest="wp", help="Analysis which conforms to the 'white papers' (Alsop et al 2014)", action="store_true", default=False)
    g.add_option("--mc", dest="mc", help="Motion correct data", action="store_true", default=False)
    parser.add_option_group(g)
    g = AslOptionGroup(parser, "Acquisition/Data specific", ignore=ignore)
    g.add_option("--casl", dest="casl", help="ASL acquisition is  pseudo cASL (pcASL) rather than pASL", action="store_true", default=False)
    g.add_option("--bolus", dest="bolus", help="Bolus duration", type=float, default=1.8)
    g.add_option("--bat", dest="bat", help="Bolus arrival time (default=0.7 (pASL), 1.3 (cASL)", type=float)
    g.add_option("--t1", dest="t1", help="Tissue T1 value", type=float, default=1.3)
    g.add_option("--t1b", dest="t1b", help="Blood T1 value", type=float, default=1.65)
    g.add_option("--slicedt", dest="slicedt", help="Timing difference between slices", type=float, default=0.0)
    g.add_option("--sliceband", dest="sliceband", help="Number of slices per pand in multi-band setup", type=int)
    g.add_option("--artsupp", dest="artsupp", help="Arterial suppression (vascular crushing) was used", action="store_true", default=False)
    parser.add_option_group(g)
    g = AslOptionGroup(parser, "Structural image (optional) (see also Registration)", ignore=ignore)
    g.add_option("--fslanat", dest="fslanat", help=" An fsl_anat directory from structural image")
    g.add_option("-s", dest="struct", help="Structural image (whole head)")
    g.add_option("--sbrain", dest="sbrain", help="Structural image (Brain extracted)")
    g.add_option("--fastsrc", dest="fastsrc", help="Images from a FAST segmentation - if not set FAST will be run on structural image")
    g.add_option("--senscorr", dest="senscorr", help="Use bias field (from segmentation) for sensitivity correction", action="store_true", default=False)
    parser.add_option_group(g)

"""
Usage_extended() {
    echo "Extended options (all optional):"
    echo " Analysis"
    echo " --artoff     : Do not infer arterial signal - (same as --artsupp)"
    echo " --fixbolus   : Bolus duration is fixed, e.g. by QUIPSSII or CASL (otheriwse it will be estimated)"
    echo "                 {default: on (cASL); off (pASL)"
    echo " --fixbat     : Fix the bolus arrival time (to value specified by --bat)"
    echo " --fulldata   : Never average multiple measurements at each TI"
    echo " --noiseprior : Use an informative prior for the noise estimation"
    echo "   --noisesd  : Set a custom noise std. dev. for the nosie prior"

    echo " Registration (requires structural image)"
    echo " --asl2struc : transformation matrix from asl data to structural image"
    echo "                (skips registration)"
    echo " --regfrom   : image to use a basis for inital registration (already BETed)"
    echo "               (must be same resolution as ASL data)"
#    echo " -t          : structural to standard space transformation matrix"
#    echo "               (requires structural image)"
#    echo " --structout : (also) save the maps in structural space (for use with -t)"
#    echo " -S          : standard brain image - {default: MNI152_T1_2mm}"
    echo " -r          : low resolution structural image (already BETed)"

    echo " Distortion correction using fieldmap (see epi_reg):"
    echo "  requires structural image to be provided"
    echo " --fmap=<image>         : fieldmap image (in rad/s)"
    echo " --fmapmag=<image>      : fieldmap magnitude image - wholehead extracted"
    echo " --fmapmagbrain=<image> : fieldmap magnitude image - brain extracted"
    echo " --echospacing=<val>    : Effective EPI echo spacing (sometimes called dwell time) - in seconds"
    echo " --pedir=<dir>          : phase encoding direction, dir = x/y/z/-x/-y/-z"
    echo " {--nofmapreg}          : do not perform registration of fmap to T1 (use if fmap already in T1-space) "
    echo " Distortion correction using phase-encode-reversed calibration image (TOPUP):"
    echo " --cblip=<image>        : phase-encode-reversed (blipped) calibration image"
    echo " --echospacing=<val>    : Effective EPI echo spacing (sometimes called dwell time) - in seconds"
    echo " --pedir=<dir>          : phase encoding direction, dir = x/y/z/-x/-y/-z"
    echo ""
 
    echo "Partial volume correction"
    echo " --pvcorr    : Do partial volume correction"
    echo "  PV estimates will be taken from:"
    echo "   fsl_anat dir (--fslanat), if supplied"
    echo "   exising fast segmentation (--fastsrc), if supplied"
    echo "   FAST segmenation of structural (if using -s and --sbet)"
    echo "   User supplied PV estimates (--pvgm, --pvwm)"   
    echo "   --pvgm    : Partial volume estimates for GM"
    echo "   --pvwm    : Partial volume estimates for WM"

    echo " Epochs "
    echo " --elen      : Length of each epoch in TIs"
    echo " --eol       : Overlap of each epoch in TIs- {deafult: 0}"
    echo ""
    echo " Miscellaneous "
    echo " --model-options=<file>  : File containing additional model options to be passed to BASIL/Fabber"
    echo "                           This is an advanced setting for model-specific options ONLY"
    echo "                           Do NOT include any standard oxford_asl options here" 
    echo ""
    echo " Notes:"
    echo " Input data - Any combination of label (tag) control images is accepted, see asl_file for more detailed"
    echo "              usage of the input options"
    echo " Analysis  -  If only one TI is specified then the following options are set:"
    echo "              --fixbat --fixbolus --artoff --fulldata"
    echo "              White Paper mode - to conform to the 'white paper' the follow options are also set:"
    echo "              --t1=1.65 --exch=simple --cmethod=voxel"
    echo " Registration - Results are saved in native space by default in case you need to revisit this."
    echo "                By default the perfusion image is used with BBR for the final registration"
    echo "                An inital registration is performed using mean of the ASL data"
    echo "                An alternative base image for (inital) registration can be supplied with --regfrom"
    echo "                Final (BBR) registration can be turned off using --finalonly"
    echo "                For custom registration: use asl_reg with native space results."
    echo " Calibration - Performed using asl_calib using CSF as a reference ('longtr' mode) where structural image is supplied"
    echo "               Voxelwise calibration is performed in 'white paper' mode or in the absence of structural image"
    echo "               For custom calibration: do not set M0 image and then use asl_calib separately"
    echo " Masking - for processing purposes a brain mask is applied to the data, this will be"
    echo "            derrived from (in order of preference):"
    echo "            > Brain extracted structural and (inital) registration"
    echo "            > mean ASL data"
    echo "            > 'regfrom' image if supplied"

}
"""

def main():
    usage = """OXFORD_ASL

    oxford_asl -i <perfusion image> -c <calibration image> --method <voxelwise|refregion> -o <output filename> [options]
    """
    debug = True
    try:
        p = OptionParser(usage=usage, version=__version__)
        add_data_options(p)
        add_oxasl_options(p)
        add_calib_options(p, ignore=["perf", "brain_mask", "output", "debug", "tis"])
        add_basil_options(p, ignore=["mask"])
        options, _ = p.parse_args()
        
        # Convert to dictionary for easier handling
        options = vars(options)
        debug = options.pop("debug", False)

        outdir = options.pop("output")
        if outdir is None:
            outdir = "oxasl"
        wsp_output = fsl.Workspace(workdir=outdir, debug=debug)
        print(wsp_output.workdir)

        if debug:
            options["wsp"] = wsp_output.sub("temp")

        options["asl_img"] = AslImage(options.pop("asldata", None), role="Input", **options)
        images = {
            "calib" : "Calibration",
        }
        for opt, role in images.items():
            if options[opt]:
                options[opt] = fsl.Image(options[opt], role=role)

        # Remove options consumed by AslImage
        for opt in ["order", "ntis", "tis", "nplds", "plds", "nrpts", "rpts", "nphases", "phases"]:
            options.pop(opt, None)

        oxasl(wsp_output, **options)
        # save command line to logfile
        #log=$outdir/logfile
        #echo $# > $log

    except Exception as e:
        sys.stderr.write("ERROR: " + str(e) + "\n")
        if debug:
            traceback.print_exc()
        sys.exit(1)

def oxasl(wsp_out, asl_img, mask=None,
          fslanat=None, struc_img=None, struc_bet_img=None, 
          calib=None, calib_ref_img=None, calib_blip_img=None,
          wsp=None, log=sys.stdout, **kwargs):
 
    # Workspace for temporary data
    if wsp is None:
        wsp = fsl.Workspace()

    # What spaces results will be output in. Always include native space
    output_spaces = {
        "native_space" : {} 
    }

    # if kwargs.get("std_trans", False):
    #     if not struc_img:
    #         raise ValueError("Structural image is required along with transformation matrix to output results in standard space")
    #     stdout = True

    log.write("OXFORD_ASL - running\n")
    log.write("Version: %s\n" % __version__)

    # Check required inputs are present
    if not asl_img:
        raise ValueError("ASL input data must be provided")
    
    log.write("Input ASL data: %s\n" % asl_img.fpath)
    asldata = wsp.add_img(asl_img, name="asldata")
    asldata.summary(log)

    diffdata = asldata.diff()
    wsp.add_img(diffdata, "diffdata")

    diffdata_mean = diffdata.mean_across_repeats()
    wsp.add_img(diffdata_mean, "diffdata_mean")

    perf_weighted = asldata.perf_weighted()
    wsp.add_img(perf_weighted, "pwi")

    if fslanat:
        # Copy over the structural and brain extracted structural
        log.write("Using FSL_ANAT output directory: %s" % fslanat)
        if (os.path.isfile(os.path.join(fslanat, "T1_biascorr")) and
            os.path.isfile(os.path.join(fslanat, "T1_biascorr_brain"))):
            log.write("Using bias-corrected structural images")
            struc_img = wsp.add_img(fsl.Image(os.path.join(fslanat, "T1_biascorr")), name="struc")
            struc_bet_img = wsp.add_img(fsl.Image(os.path.join(fslanat, "T1_biascorr_brain")), name="struc_bet")
        else:
            log.write("Using non bias-corrected structural images")
            struc_img = wsp.add_img(fsl.Image(os.path.join(fslanat, "T1")), name="struc")
            struc_bet_img = wsp.add_img(fsl.Image(os.path.join(fslanat, "T1_brain")), name="struc_bet")
        
        output_spaces["struct_space"] = {
            "xfm" : None # FIXME
        }
        
        if os.path.isfile(os.path.join(fslanat, "T1_to_MNI_nonlin_coeff")):
            output_spaces["std_space"] = {
                "warp" : os.path.join(fslanat, "T1_to_MNI_nonlin_coeff")
            }
        elif os.path.isfile(os.path.join(fslanat, "T1_to_MNI_lin.mat")):
            output_spaces["std_space"] = {
                "xfm" : os.path.join(fslanat, "T1_to_MNI_lin")
            }

    elif struc_img:
        output_spaces["struct_space"] = None
        log.write("Using structural image: %s\n" % struc_img.fpath)
        wsp.add_img(struc_img, "struc")
        if not struc_bet_img:
            log.write("Running BET on structural image\n")
            struc_bet_img = wsp.bet(struc_img, "struc_bet")
        else:
            struc_bet_img = wsp.add_img(struc_bet_img, "struc_bet")

        output_spaces["struct_space"] = {
            "xfm" : None # FIXME
        }

    elif kwargs.get("struc_img_lowres", None):
        struc_img_lowres = kwargs.pop("struc_img_lowres")
        log.write("Low-resolution tructural image: %s\n" % struc_img_lowres.fpath)     
        struc_img_lowres = wsp.add_img(struc_img_lowres, "lowstruc")

    elif "std_xfm" in kwargs or "std_warp" in kwargs:
        output_spaces["std_space"] = {
            "xfm" : kwargs.pop("std_xfm", None),
            "warp" : kwargs.pop("std_warp", None),
            "brain" : kwargs.pop("std_brain", None),
        }
        
        if not output_spaces["std_space"]["std_brain"]:
            output_spaces["std_space"]["std_brain"] = os.path.join(fsldir, "data", "standard", "MNI152_T1_2mm")
        log.write("Standard brain is: %s" % output_spaces["std_space"]["std_brain"])
    
    else:
        log.write("No structural data found - output will be in native space only\n")

    log.write("\nPre-processing\n")
       
    # Pre-process calibration images
    calib = preproc_calib(wsp, calib, "calib", "Calibration")
    calib_ref_img = preproc_calib(wsp, calib_ref_img, "cref", "Calibration reference")
    calib_blip_img = preproc_calib(wsp, calib_blip_img, "cblip", "Calibration BLIP")

    if not mask:
        mask = generate_mask(wsp, meanasl=perf_weighted, struc_bet_img=struc_bet_img, calib_img=calib, log=log)
    wsp.add_img(mask, name="mask")

    if kwargs.pop("mc", False): 
        motion_correct(wsp, asldata, ref=calib, log=log)

    # Single or Multi TI setup
    if asldata.ntis == 1:
        # Single TI data - don't try to infer arterial component of bolus duration, we don't have enough info
        log.write("Operating in Single TI mode\n")
        kwargs["artoff"] = True
        kwargs["fixbolus"] = True
        singleti = True
    
    casl = kwargs.pop("casl", False)

    if kwargs.get("wp", False):
        # White paper mode - this overrides defaults, but can be overwritten by command line 
        # specification of individual parameters
        log.write("Analysis in white paper mode\n")
        t1_default = 1.65
        bat_default = 0.0
        calib_method = "voxel"
    else:
        t1_default = 1.3
        if casl:
            bat_default = 1.3
        else:
            bat_default = 0.7

    kwargs["t1"] = kwargs.get("t1", t1_default)
    kwargs["t1b"] = kwargs.get("t1b", 1.65)
    kwargs["bat"] = kwargs.get("bat", bat_default)
    kwargs["bolusdur"] = kwargs.get("bolusdur", 1.8)

    # if we are doing CASL then fix the bolus duration, except where the user has 
    # explicitly told us otherwise
    kwargs["fixbolus"] = kwargs.get("fixbolus", casl)
        
    # Initial BASIL run on mean data
    log.write("\nDoing initial fit on mean at each TI\n")
    run_basil(wsp, "init", diffdata_mean, mask=mask, log=log, **kwargs)

    log.write("Output is %s\n" % wsp_out.workdir)
    log.write("OXFORD_ASL - done\n")

def preproc_calib(wsp, img, output_name, desc, log=sys.stdout):
    """
    Standard pre-processing of calibration images
    
    assumed to be in native space but not necessarily motion corrected or 
    aligned to ASL data
    """
    if img:
        log.write("%s image from %s\n" % (desc, img.fpath))
        data = img.data()
        if img.ndim == 4:
            if img.shape[3] > 1:
                log.write("Removing first volume to ensure data is in steady state" % desc)
                data = data[..., :-1]
                    
                #if [ ! -z $moco ]; then
                #    #motion correction
                #    mcflirt -in $tempdir/calib -o $tempdir/calib
                #fi
            
            # Take the mean across time axis
            data = np.mean(data, axis=-1)

        return wsp.add_img(img.derived(data=data, name=output_name))

def generate_mask(wsp, meanasl, struc_bet_img=None, regfrom_img=None, calib_img=None, struc_lowres_bet_img=None, log=sys.stdout):
    """
    Generate mask for ASL data
    """
    log.write("Automatic mask generation\n")
    if struc_bet_img:
        # Preferred option is to use brain extracted structural
        struc_bet_mask = wsp.maths(struc_bet_img, args="-bin")
        mask = wsp.apply_xfm(struc_bet_mask, regfrom_img, "struct2asl.mat", args="-interp trilinear")
        mask = wsp.maths(mask, args="-thr 0.25 -bin -fillh")
        #wsp.fslcpgeom(regfrom_img, mask) FIXME
        log.write("Mask generated from brain extracted structural image: %s\n" % struc_bet_img.fpath)
    elif regfrom_img:
        # Otherwise use the regfrom image (should already be BETed) - note that regfrom may have been set in the Registration section.
        mask = wsp.maths(regfrom_img, args="-bin")
        log.write("Mask generated from regfrom image: %s" % regfrom_img.fpath)
    elif calib_img:
        # Next option is to use betted version of mean M0 calib image as mask
        mask = wsp.bet(calib_img, mask=True, brain=False)
        log.write("Mask generated from calibration image (post BET)\n")
    #elif [ ! -z $lowstrucflag ]; then
    #    # Use the low resolution strucutral image to create mask (ahould already be BETed)
    #    flirt -in $tempdir/lowstruc -applyxfm -init $FSLDIR/etc/flirtsch/ident.mat -out $tempdir/mask -paddingsize 0.0 -interp trilinear -ref $tempdir/asldata
    #    fslmaths $tempdir/mask -bin $tempdir/mask
    #    log.write("Mask generated from low res. structural\n")
    else:
        # otherwise just use mean time series - use a fairly low fraction value to avoid erosion
        mask = wsp.bet(meanasl, mask=True, brain=False, args="-f 0.2")
        log.write("Mask generated from mean time series\n")

    return mask

def run_basil(wsp, subdir, asldata, mask, log=sys.stdout, **kwargs):
    steps = basil.get_steps(asldata, mask=mask, log=log, **kwargs)
    basil_wsp = wsp.sub(subdir)
    basil.run_steps(basil_wsp, steps, log)
    final_step = 1
    while 1:
        if not os.path.isdir(os.path.join(basil_wsp.workdir, "step%i" % (final_step+1))):
            break
        final_step += 1
    basil_wsp.imcp(os.path.join(basil_wsp.workdir, "step%i/finalMVN" % final_step), os.path.join(basil_wsp.workdir, "finalMVN"))

# # extract images from BASIL results (and throw away values below zero)
# fslmaths ${finalstep}/mean_ftiss -thr 0 $2/ftiss
# if [ ! -z $senscorr ]; then
# # sensitivity correction
#     fslmaths $2/ftiss -div $outdir/native_space/sensitivity $2/ftiss
# fi

# if [ -z $fixbat ]; then
#     fslmaths ${finalstep}/mean_delttiss -thr 0 $2/delttiss
# fi
# if [ -z $artoff ]; then
#     fslmaths ${finalstep}/mean_fblood -thr 0 $2/fblood
#     if [ ! -z $senscorr ]; then
# # sensitivity correction
# 	fslmaths $2/fblood -div $outdir/native_space/sensitivity $2/fblood
#     fi
# fi

# #Partial volume correction - sort out basil results when PV corrected
# if [ ! -z $pvcorr ]; then
#     fslmaths ${finalstep}/mean_fwm -thr 0 $2/ftisswm
#     if [ ! -z $senscorr ]; then
#         # sensitivity correction
# 	fslmaths $2/ftisswm -div $outdir/native_space/sensitivity $2/ftisswm
#     fi
#     if [ -z $fixbat ]; then
# 	fslmaths ${finalstep}/mean_deltwm -thr 0 $2/deltwm
#     fi

# fi

# if [ ! -z $varout ]; then
# #get varainces out of finalMVN
#     fabber_var -d ${finalstep} -m $tempdir/mask
# # do correction of negative values
#     fslmaths ${finalstep}/var_ftiss -bin -add 1 -uthr 1 -mul 1e12 -add ${finalstep}/var_ftiss $2/var_ftiss
#     if [ ! -z $senscorr ]; then
# # sensitivity correction
#     fslmaths $2/var_ftiss -div $outdir/native_space/sensitivity -div $outdir/native_space/sensitivity $2/var_ftiss
# fi
#     if [ -z $fixbat ]; then
# 	fslmaths ${finalstep}/var_delttiss -bin -add 1 -uthr 1 -mul 1e12 -add ${finalstep}/var_delttiss $2/var_delttiss
#     fi
# fi

# #copy the final MVN to the temp directory for future use
# imcp ${finalstep}/finalMVN $2/finalMVN
# cp ${finalstep}/paramnames.txt $2/paramnames.txt


# }

def output(wsp, basil_name, output_name, spaces):
    pass
    # for sp, xfm in spaces:
    #     output_wsp = wsp.sub(sp)
    #     if xfm:
    #         reg_img = wsp.apply_xfm(basil_name, struc, xfm)
    #         wsp.add_img(reg_img, name=output_name)
    #     else:
    #         wsp.add_img(basil_img, name=output_name)

#def output_masked(output_name, subdir, mask):
#    # Function to output images having been masked
#    # currently we only do this in native space
#    if nativeout:
#    	fslmaths $outdir/native_space/$subdir/$output_name -mas $tempdir/${maskname}mask $outdir/native_space/$subdir/${output_name}_masked

#def calibrate(basil_name, output_name, m0val, multiplier, subdir):
#    fslmaths  $tempdir/$subdir/$param -div $m0val -mul $multiplier $tempdir/$subdir/${param}_calib 
#    Output ${param}_calib ${parname}_calib $subdir

#def report(output_name, subdir, masktype):
#    if pvexist:
#    	# generate text reports on parameters - the parameter must have been output first for this to work#
#	    #NB we only do this is the PVE are available (and thus the reqd masks will exist)
#
#        repval=`fslstats $outdir/native_space/$subdir/$parname -k $tempdir/${masktype}mask_pure -m`
#        echo $repval > $outdir/native_space/$subdir/${parname}_${masktype}_mean.txt
#        Log "Mean $parname in $masktype is $repval"

#def normalise(output_name, subdir, masktype)
#    if pvexist:
#        # also output the perfusion images normalised by the mean mask value - the parameter must have been output first for this to work
#        #NB we only do this is the PVE are available (and thus the reqd masks will exist)
#        
#        # get normalization from reported value in the output directory
#        normval=`cat $outdir/native_space/$subdir/${parname}_${masktype}_mean.txt`
#        
#        if stdout:
#            fslmaths $outdir/std_space/$subdir/$parname -div $normval $outdir/std_space/$subdir/${parname}_norm 
#        if nativeout:
#            fslmaths $outdir/native_space/$subdir/$parname -div $normval $outdir/native_space/$subdir/${parname}_norm
#        if structout:
#            fslmaths $outdir/struct_space/$subdir/$parname -div $normval $outdir/struct_space/$subdir/${parname}_norm 

# def normalise_var(param, subdir, masktype):
#     if pvexist:
#         # normalisaiton for a variance image

#         # get normalization from reported value in the output directory
#         normval=`cat $outdir/native_space/$subdir/${parname}_${masktype}_mean.txt`
#         #need to square the value as we are outputting a variance
#         normval=`echo "$normval * $normval" | bc`
        
#         if stdout:
#             fslmaths $outdir/std_space/$subdir/${parname}_var -div $normval $outdir/std_space/$subdir/${parname}_var_norm 
#         if nativeout:
#             fslmaths $outdir/native_space/$subdir/${parname}_var -div $normval $outdir/native_space/$subdir/${parname}_var_norm
#         if structout:
#             fslmaths $outdir/struct_space/$subdir/${parname}_var -div $normval $outdir/struct_space/$subdir/${parname}_var_norm 

# def registration(regbase, transopt, distout):
#     echo "Performing registration"
#     regbase=$1 #the i/p to the function is the image to use for registration
#     transopt=$2 # other options to pass to asl_reg
#     distout=$3 # we want to do distortion correction and save in the subdir distout
    
#     extraoptions=" "
#     if [ ! -z $lowstrucflag ]; then
# 	extraoptions=$extraoptions"-r $tempdir/lowstruc"
#     fi
#     if [ ! -z $debug ]; then
# 	extraoptions=$extraoptions" --debug"
#     fi
    
#     #if [ ! -z $reginit ]; then
#     #    extraoptions=$extraoptions" -c $reginit"
#     #fi
    
#     if [ -z $distout ]; then
# 	# normal registration
# 	$asl_reg -i $regbase -o $tempdir -s $tempdir/struc --sbet $tempdir/struc_bet $transopt $extraoptions

# 	if [ ! -z $trans ]; then
# 	    # compute the transformation needed to standard space if we have the relvant structural to standard transform
# 	    convert_xfm -omat $tempdir/asl2std.mat -concat $trans $tempdir/asl2struc.mat
# 	fi
	    
#     else
# 	# registration for distortion correction
# 	fmapregstr=""
# 	if [ ! -z $nofmapreg ]; then
# 	    fmapregstr="--nofmapreg"
# 	fi
# 	$asl_reg -i $regbase -o $tempdir/$distout -s $tempdir/struc --sbet $tempdir/struc_bet $transopt $extraoptions --fmap=$tempdir/fmap --fmapmag=$tempdir/fmapmag --fmapmagbrain=$tempdir/fmapmagbrain --pedir=$pedir --echospacing=$echospacing $fmapregstr
#     fi


# def calibration() {
#  echo "Calculating M0a - calling ASL_CALIB"
#     extraoptions=""
#     if [ ! -z $debug ]; then
# 	extraoptions=$extraoptions" --debug"
#     fi

#     #if [ ! -z $cref ]; then
# 	# pass calibration reference image to asl_calib
# 	#extraoptions=$extraoptions" --cref $tempdir/cref"
#     if [ ! -z $senscorr ]; then
# 	    # use a sensitivity iamge from elsewhere
# 	    Log "Sensitivity image $outdir/native_space/sensitivity being loaded into asl_calib"
# 	    extraoptions=$extraoptions" --isen $outdir/native_space/sensitivity"
#     fi

#     if [ -z $te ]; then
# 	#by default assume TE is zero
# 	te=0
#     fi

#     if [ ! -z $t2star ]; then
# 	# tell asl_calib to correct for T2* rather than T2
# 	extraoptions=$extraoptions" --t2star"
#     fi

#     if [ ! -z $tissref ]; then
# 	# Specify reference tissue type
# 	extraoptions=$extraoptions" --tissref $tissref"
#     fi

#     if [ ! -z $t1csf ]; then
# 	# supply the T1 of csf
# 	extraoptions=$extraoptions" --t1r $t1csf"
#     fi

#     if [ ! -z $t2csf ]; then
# 	# Supply the T2(*) of CSF
# 	extraoptions=$extraoptions" --t2r $t2csf"
#     fi

#     if [ ! -z $t2bl ]; then
# 	# Supply the T2(*) of blood
# 	extraoptions=$extraoptions" --t2b $t2bl"
#     fi

#     if [ ! -z $debug ]; then
# 	#run asl_calib in debug mode
# 	extraoptions=$extraoptions" --debug"
#     fi

#     # setup the main options that we will pass to aslcalib regardless of whether we are auot generating reference mask
#     maincaliboptions="--cgain $cgain --te $te --tr $tr"

#     if [ -z $csfflag ]; then
#     # call asl_calib in normal (auto csf) mode

#     # use low res structural for auto generation of csf mask if availible
#     # otherwise just use structural image
# #	if [ -z $lowstrucflag ]; then
# #	    usestruc=$tempdir/struc_bet
# #	    usetrans=$tempdir/asl2struct.mat
# #	else
# #	    usestruc=$tempdir/lowstruc_bet
# #	    usetrans=$tempdir/asl2lowstruct.mat
# #	fi

# 	usestruc=$tempdir/struc_bet
# 	usetrans=$tempdir/asl2struct.mat

# 	if [ ! -z $fasthasrun ]; then
# 	    # we have already run FAST so we can pass the PVE for CSF to asl_calib (to save running FAST again)
# 	    extraoptions=$extraoptions" --refpve $tempdir/pvcsf_struct"
# 	fi
   
# 	$asl_calib -c $calib -s $usestruc -t $usetrans -o $outdir/calib --bmask $tempdir/mask --osen $outdir/native_space/sensitivity $maincaliboptions $extraoptions 

#     else
#     # a manual csf mask has been supplied
# 	$asl_calib -c $calib -m $csf -o $outdir/calib --bmask $tempdir/mask --osen $outdir/native_space/sensitivity $maincaliboptions $extraoptions
#     fi
# }

# Dooutput() {
# # Do all the outputs - using the supplied subdirectiory of the results

# if [ -z $1 ]; then
#     subdir=/ #need a default 'empty' value for this
# else
#     subdir=$1
# fi

# # perfusion
# Output ftiss perfusion $subdir
# Report perfusion $subdir gm
# Normalise perfusion $subdir gm

# # arrival
# if [ -z $fixbat ]; then
#     Output delttiss arrival $subdir
#     Report arrival $subdir gm
# fi
# # aBV
# if [ -z $artoff ]; then
#     Output fblood aCBV $subdir
# fi

# # white matter values
# if [ $subdir = "pvcorr" ]; then
#     Output ftisswm perfusion_wm $subdir
#     Report perfusion_wm $subdir wm
#     Normalise perfusion_wm $subdir wm
#     if [ -z $fixbat ]; then
# 	Output deltwm arrival_wm $subdir
# 	Report arrival_wm $subdir wm
#     fi
    
# else 
#     Report perfusion $subdir wm
#     if [ -z $fixbat ]; then
# 	Report arrival $subdir wm
#     fi
# fi

# # Masked results (PVcorr)
# if [ $subdir = "pvcorr" ]; then
#     OutputMasked perfusion $subdir gm
#     OutputMasked perfusion_wm $subdir wm
#     if [ -z $fixbat ]; then
# 	OutputMasked arrival $subdir gm
# 	OutputMasked deltwm $subdir wm
#     fi
# fi

# #Optionally provide variance results
# if [ ! -z $varout ]; then
#     Output var_ftiss perfusion_var $subdir
#     Normalise_var perfusion $subdir gm
#     if [ -z $fixbat ]; then
#     Output var_delttiss arrival_var $subdir
#     fi
    
# fi

# # calibrated results
# if [ ! -z $calibflag ]; then
#     if [ $cmethod = 'single' ]; then
# 	malpha=`echo "$Mo * $alpha" | bc` #include the inversion efficiency when we do the final calibration
#     elif [ $cmethod = 'voxel' ]; then
# 	fslmaths $outdir/calib/M0 -mul $alpha $tempdir/malpha
# 	malpha=$tempdir/malpha
#     fi

#     Calibrate ftiss perfusion $malpha 6000 $subdir
#     Report perfusion_calib $subdir gm
   
#     if [ $subdir = "pvcorr" ]; then
# 	OutputMasked perfusion_calib $subdir gm
# 	Calibrate ftisswm perfusion_wm $malpha 6000 $subdir
# 	Report perfusion_wm_calib $subdir wm
# 	OutputMasked perfusion_wm_calib $subdir wm
#     else
# 	Report perfusion_calib $subdir wm
#     fi

#     if [ ! -z $varout ]; then
# 	if [ $cmethod = 'single' ]; then
# 	    Mosq=`echo "$Mo * $Mo * $alpha * $alpha" | bc` #include the inversion efficiency when we do the final calibration
# 	elif [ $cmethod = 'voxel' ]; then
# 	    fslmaths $outdir/calib/M0 -mul $outdir/calib/M0 -mul $alpha -mul $alpha $tempdir/mosq
#         Mosq=$tempdir/mosq
# 	fi
	
# 	Calibrate var_ftiss perfusion_var $Mosq 36000000 $subdir
#     fi

#     if [ -z $artoff ];then
#         # output aCBV as a percentage
# 	    Calibrate fblood aCBV $malpha 100 $subdir
#     fi
# fi

# # advanced output
# if [ ! -z $advout ]; then
#     if [ ! -d  $outdir/advanced/$subdir ]; then mkdir $outdir/advanced/$subdir; fi
#     imcp $tempdir/$subdir/finalMVN $outdir/advanced/$subdir/finalMVN
#    cp $tempdir/$subdir/paramnames.txt $outdir/advanced/$subdir/paramnames.txt
# fi

# }

def motion_correct(wsp, asldata, calib_img=None, cref_img=None, cblip_img=None, log=sys.stdout):
    """
    Motion Correction (main)
    
    Note motion correction within calibration data is done above
    """
    log.write("Motion Correction\n")
    if calib_img:
        # use supplied image as our reference for motion correction
        # Normally the calibration image since this will be most consistent if the data has a range of different TIs and background suppression etc
        # this also removes motion effects between asldata and calibration image
        log.write("Using calibration image as reference\n")
        asldata_mc, mats = wsp.mcflirt(asldata, ref=calib_img, mats=True)

        # To reduce interpolation of the ASL data change the transformations so that we end up in the space of the central volume of asldata
        asl2calib = mats[len(mats)/2+1]
        calib2asl = np.linalg.inv(middlemat)
        log.write("middle volume->calib:\n%s\n" % str(asl2calib))
        log.write("calib->middle volume:\n%s\n" % str(calib2asl))

        # Convert all the volumes to this space
        asl_moco_mats = [np.dot(mat, calib2asl) for mat in mats]
        wsp.apply_xfm_4d(asldata, asldata_mc, asl_moco_mats, output_name=asldata.iname + "_mc", args="-fourdigit")

        # Convert all calibration images to align with asldata
        wsp.apply_xfm(calib_img, asldata, xfm=calib2asl)
        if cref_img:
            wsp.apply_xfm(cref_img, asldata, xfm=calib2asl)
        if cblip_img:
            wsp.apply_xfm(cblip_img, asldata, xfm=calib2asl)

    else:
         log.write("Using ASL data middle volume as reference\n")
         mcdata, mats = wsp.mcflirt(asldata, mats=True)
         
    #cat $tempdir/asldata.mat/MAT* > $tempdir/asldata.cat # save the motion matrices for distortion correction if reqd

def get_asl2struc(wsp, asldata, regfrom):
    """
    Registration (1/2)
    # Make sure we have some form of transformation between the ASL data and the structural (if that has been supplied)
    # only 'initial' step in asl_reg is used here
    """
    pass
    #register=0
    #if [ ! -z $strucflag ]; then # if structural image has not been suppled then skip the registration
    #    register=1
    #    if [ ! -z $asl2struc ]; then # we have been supplied with a transformation matrix - we do not need registration, but we do want to transform the results
    #    register=0
    #    Log "Using existing asl to structural transform: $asl2struc"
    #    cp $asl2struc $tempdir/asl2struct.mat
    #    convert_xfm -omat $tempdir/struct2asl.mat -inverse $tempdir/asl2struct.mat
    #    fi
    #fi

    # if not regfrom:
    #     # No regfrom iamge supplied so we will use the mean of the asl timeseries - unless it was diff data
    #     # NB in the case of really good background suppresion this might not be best option even if raw ASL data has 
    #     # been supplied, in which case a calibration image could be provided to regfrom (or failing that the PWI). 
    #     # We wont force that here as we dont know.
    #     if calib_img and "p" not in asldata.order and "P" not in asldata.order:
    #         #if available use calibration image
    #         regfrom_img = wsp.bet(calib_img, "calib_brain")
    #     else:
    #         regfrom_img = wsp.bet(meanasl_img, "meanasl_brain", "-f 0.2")

    # #if [ $register -eq 1 ]; then
    # # registration here using asl_reg (inital only)
    # log.write("Performing registration\n")
    # log.write("Using %s as base for inital registration" % regfrom_img.iname)
    # extraoptions="--mainonly " # to ensure we only do the initil flirt part
    # if low_struc:
    #     extraoptions += "-r $tempdir/lowstruc"

    # Registration(regfrom_img, extraoptions)
    # #convert_xfm -omat $tempdir/struct2asl.mat -inverse $tempdir/asl2struct.mat

def segment():
    pass
    ### Segmentation of structural image - if we have a structural image we ALWAYS ensure we have a segmentation
    # if [ ! -z $fslanat ]; then
    #     # we are being supplied with an fslanat directory
    #     fasthasrun=1 #this means that we have PVE for calibration & PVC purposes
        
    #     # copy over the things we need and place them using the names used elsewhere
    #     imcp $fslanat/T1_fast_pve_0 $tempdir/pvcsf_struct #indicate that it is in structural space!
    #     imcp $fslanat/T1_fast_pve_1 $tempdir/pvgm_struct
    #     imcp $fslanat/T1_fast_pve_2 $tempdir/pvwm_struct

    #     if [ ! -z $fslanat/T1_fast_bias ]; then # test to check that there is a bias field here
    #     Log "Bias field extracted from $fslanat sucessfully"
    #     imcp $fslanat/T1_fast_bias $tempdir/biasfield_struct
    #     else
    #     Log "No Bias field found in $fslanat"
    #     fi
        
    # elif [ ! -z $struc ]; then
    #     # do we have the results from FAST already? If not run it
    #     if [ -z $fastsrc ]; then
    #     echo "Segmenting the structural image"
    #     Log "Segmenting the structural image"
    #     fast -B -b -o $tempdir/seg -p $tempdir/struc_bet
    #     fastsrc=$tempdir/seg
    #     else
    #     # FAST has been run externally
    #     Log "Using FAST outputs at: $fastsrc"
    #     fi
        
    #     # we are now sure we have FAST outputs
    #     fasthasrun=1

    #     # copy over the things we need and place them using the names used elsewhere
    #     imcp ${fastsrc}_pve_0 $tempdir/pvcsf_struct #indicate that it is in structural space!
    #     imcp ${fastsrc}_pve_1 $tempdir/pvgm_struct
    #     imcp ${fastsrc}_pve_2 $tempdir/pvwm_struct

    #     if [ ! -z ${fastsrc}_bias ]; then # test to see if there is a bias field in the FAST output
    #     Log "Bias field extracted from ${fastsrc} sucessfully"
    #     imcp ${fastsrc}_bias $tempdir/biasfield_struct
    #     else
    #     Log "No Bias field found with ${fastsrc}"
    #     fi
        
    # fi

    # # some useful preproc to do with FAST outputs
    # if [ ! -z $fasthasrun ]; then
    #     # create a tissseg (wmseg) image for BBR in asl_reg
    #     fslmaths $tempdir/pvwm_struct -thr 0.5 -bin ${tempdir}/tissseg

    #     if [ ! -z $tempdir/biasfield_struct ]; then
    #     # transform the bias field and invert to use for sensitivity correction in calibration
    #     applywarp --ref=$tempdir/asldata --in=$tempdir/biasfield_struct --out=$tempdir/biasfield --premat=$tempdir/struct2asl.mat --super --interp=spline --superlevel=4
    #     fi

    #     if [ ! -z $senscorr ]; then
    #     if [ ! -z $tempdir/biasfield ]; then #make sure we have the biasfield (from above) before attempting this
    #         Log "Creating sensitivity map from biasfield"
    #         fslmaths $tempdir/biasfield -recip $outdir/native_space/sensitivity
    #     fi	    
    #     fi

    # fi

def distcorr():
    pass
    ### Distortion Correction
    # Do TOPUP if applicable
    # if [ ! -z $cblip ]; then
    #     if [ -z $calib ]; then
    #     echo "WARNING: Cannot do TOPUP on blip-reversed calibration image ($cblip) without correpsonding calibration image"
    #     elif [ -z $echospacing ] || [ -z $pedir ]; then
    #     echo "WARNING: Cannot do TOPUP on blip-reversed calibration image without echospacing (dwell time) and pahse encode direction"
    #     else
    #     echo "Distortion correction: running topup"
        
    #         #create topup params
    #     case $pedir in
    #         x)
    #         echo "1 0 0 $echospacing" > $tempdir/topup_params.txt
    #         echo "-1 0 0 $echospacing" >> $tempdir/topup_params.txt
    #         ;;
    #             -x)
    #         echo "-1 0 0 $echospacing" > $tempdir/topup_params.txt
    #         echo "1 0 0 $echospacing" >> $tempdir/topup_params.txt
    #         ;;
    #         y)
    #         echo "0 1 0 $echospacing" > $tempdir/topup_params.txt
    #         echo "0 -1 0 $echospacing" >> $tempdir/topup_params.txt
    #         ;;
    #         -y)
    #         echo "0 -1 0 $echospacing" > $tempdir/topup_params.txt
    #         echo "0 1 0 $echospacing" >> $tempdir/topup_params.txt
    #         ;;
    #         z)
    #         echo "0 0 1 $echospacing" > $tempdir/topup_params.txt
    #         echo "0 0 -1 $echospacing" >> $tempdir/topup_params.txt
    #         ;;
    #             -z)
    #         echo "0 0 -1 $echospacing" > $tempdir/topup_params.txt
    #         echo "0 0 1 $echospacing" >> $tempdir/topup_params.txt
    #         ;;
    #     esac
        
    #     # do topup
    #     fslmerge -t $tempdir/calib_blipped $tempdir/calib $cblip 
    #     topup --imain=$tempdir/calib_blipped --datain=$tempdir/topup_params.txt --config=b02b0.cnf --out=$tempdir/topupresult --fout=$tempdir/topupresult_fmap
    #     topupresult=$tempdir/topupresult
    #     fi
    # fi

    # #Fieldmaps
    # if [ ! -z $topupresult ]; then
    # #    if [ -e $tempdir/struc.* ]; then
    # #	# currently DISABLED and applytopup is used with topup results
    # #	#we will do the distorition correction using epi_reg so that we can merge with motion correction matrices and also get the jacobian
    # #	# the fieldmap provided is from topup and will be in ASL space
    # #	# convert ot the correct units of rad/s from Hz
    # #	fslmaths ${topupresult}_fmap -mul 3.1459 -mul 2 $tempdir/fmap
    # #	# use the existing registration to get the fieldmap into T1 space
    # #	flirt -in $tempdir/fmap -out $tempdir/fmap -ref $tempdir/struc -applyxfm -init $tempdir/asl2struct.mat
    # #	# asl_reg/epi_reg will expect a fieldmap magnitude image (although doesn't really need it in this case) - just use the structural
    # #	imcp $tempdir/struc $tempdir/fmapmag
    # #	imcp $tempdir/struc_bet $tempdir/fmapmagbrain
    # #	nofmapreg=1
    # #    else
    #     echo "Distortion Correction using TOPUP"
    #     # we will use apply topup - this does not do the jacboian magntiude correction - therefore strictly only okay if using voxelwise calibration
    #     applytopup --imain=$tempdir/calib,$cblip --inindex=1,2 --datain=$tempdir/topup_params.txt --topup=${topupresult} --out=$tempdir/calib --method=jac
    #     repeatsubtract=1;
    #     applytopup --imain=$tempdir/asldata --datain=$tempdir/topup_params.txt --inindex=1 --topup=${topupresult} --out=$tempdir/asldata --method=jac #ND using asldata as this has been motion corrected by this point (if requested)
    #     if [ ! $cmethod="voxel" ]; then
    #         echo "WARNING: Using apply_topup this does not correct for magntiude using the jocbian in distortion correction - this is not optimal when not using voxelwise calibration, to avoid this supply structural image(s)"
    #     fi
        
    #     if [ ! -z $cref ]; then
    #         applytopup --imain=$tempdir/cref --datain=$tempdir/topup_params.txt --inindex=1 --topup=${topupresult} --out=$tempdir/cref --method=jac
    # #	fi
    #     fi
    # elif [ ! -z $fmap ]; then
    #     # a fieldmap has been provided that needs registration - copy images over
    #     imcp $fmap $tempdir/fmap
    #     imcp $fmapmag $tempdir/fmapmag
    #     imcp $fmapmagbrain $tempdir/fmapmagbrain
    # fi

    # if [ -e $tempdir/fmap.* ]; then
    #     echo "Distortion Correction using asl_reg"
        
    #     # Do registration to T1 to get distortion correction warp
    #     # use whatever registration matrix we already have to initialise here 
    #     distbase=$tempdir/pwi # use the perfusion-weighted image (mean over all TIs) as the best basis we have for registration at this point
    #     if [ -z $finalreg ]; then
    #     distbase=$regfrom # use whatever image we have been using for (inital) registration
    #     fi

    #     Registration $distbase "-m $tempdir/mask --tissseg $tempdir/tissseg --imat $tempdir/asl2struct.mat --finalonly" distcorr
        
    #     # generate the correction warp
    #     convertwarp -r $tempdir/meanasl -o $tempdir/asldist_warp -w $tempdir/distcorr/asl2struct_warp.nii.gz --postmat=$tempdir/distcorr/struct2asl.mat --rel -j $tempdir/distcorr/jacobian_parts
    #     fslmaths $tempdir/distcorr/jacobian_parts -Tmean $tempdir/distcorr/jacobian
        
    #     # Now apply the correction to the data.
    #     # note that we use the orignal data here and apply the motion correction as part of the process
    #     appremat=""
    #     if [ ! -z $moco ]; then
    #     appremat="--premat=$tempdir/asldata.cat"
    #     fi
    #     applywarp -i $tempdir/asldata_orig -r $tempdir/meanasl -o $tempdir/asldata $appremat -w $tempdir/asldist_warp --rel --interp=spline --paddingsize=1
    #     fslmaths $tempdir/asldata -mul $tempdir/distcorr/jacobian $tempdir/asldata
    #     repeatsubtract=1;
            
    #     # Now apply the correction to the calibration image
    #     applywarp -i $tempdir/calib -r $tempdir/calib -o $tempdir/calib -w $tempdir/asldist_warp --rel --interp=spline --paddingsize=1
    #     fslmaths $tempdir/calib -mul $tempdir/distcorr/jacobian $tempdir/calib
    
    #     if [ ! -z $cref ]; then
    #     applywarp -i $tempdir/cref -r $tempdir/calib -o $tempdir/cref -w $tempdir/asldist_warp --rel --interp=spline --paddingsize=1
    #     fslmaths $tempdir/cref -mul $tempdir/distcorr/jacobian $tempdir/cref
    #     fi
    # fi

    # # Repeat the label-control subtraction on the corrected data
    # if [ ! -z $repeatsubtract ]; then
    #     if [ $iaf = 'diff' ]; then
    #     # make sure the block format is correct for BASIL
    #     asl_file --data=$tempdir/asldata --ntis=$ntis --ibf=$ibf --iaf=$iaf --obf=tis --out=$tempdir/diffdata --mean=$tempdir/diffdata_mean
    #     else
    #     # create label-control difference data using asl_file - this gets it into the right block form for BASIL (blocks of TIs)
    #     asl_file --data=$tempdir/asldata --ntis=$ntis --ibf=$ibf --iaf=$iaf --obf=tis --diff --out=$tempdir/diffdata --mean=$tempdir/diffdata_mean
    #     fi
    # fi
    # ### End of: Distortion Correction

    # Mask the calibration images, saving the wholehead images (although currently unused)
    #if [ ! -z $calib ]; then
    #    imcp $tempdir/calib $tempdir/calib_wholehead
    #    fslmaths $tempdir/calib -mas $tempdir/mask $tempdir/calib
    #fi
    #if [ ! -z $cref ]; then
    #    imcp $tempdir/cref $tempdir/cref_wholehead
    #    fslmaths $tempdir/cref -mas $tempdir/mask $tempdir/cref
    #fi

    # Copy or recalculate the sensitivity map (overwriting any created
    # by FAST) if we have a cref image or a sensitivity image has been supplied directly
    #if [ ! -z $isen ]; then
    #    # User-supplied sensitivity image
    #    Warn "User supplied sens image"
    #    fslmaths $isen -mas $tempdir/mask $outdir/native_space/sensitivity
    #elif [ ! -z $cref ]; then
    #    # Calculate sensitivty image using user-supplied cref image
    #    fslmaths $tempdir/calib -div $tempdir/cref -mas $tempdir/mask $outdir/native_space/sensitivity
    #fi

    # Senstivity correction cannot be done if this image hasn't been generated by this point
    #if [ -z $outdir/native_space/sensitivity ]; then
    #    senscorr=""
    #    Warn "sensitivity correction has been requested, but suitable map is not available, skipping that step"
    #fi

    #pre-processing for epochwise analysis
    # separate data into epochs here
    #if [ ! -z $epoch ]; then
    #if [ -z $eol ]; then
    #    eol=0
    #fi
    #    asl_file --data=$tempdir/diffdata --ntis=$ntis --ibf=tis --iaf=diff --epoch=$tempdir/epoch --elen=$elen --eol=$eol --eunit=tis
    #    eadv=`expr $elen - $eol`
    #fi

    # write options file for BASIL - these are the core options that are appropraite whether we are doing a single or epochwise analysis
    
#def basil(wsp, log=sys.stdout):
#    log.write("Setting up BASIL\n")

    # model_options = {}
    # # data acquired using CASL?
    # if casl:
    #     model_options["casl"] = True
    #     log.write("cASL model\n")
    # else
    #     log.write("pASL model\n")

    # #echo $tauslist >> $tempdir/basil_options.txt

    # # slice timing correction?
    # if [ ! -z $slicedt ]; then
    #     echo "--slicedt=$slicedt" >> $tempdir/basil_options.txt
    #     Log "Slice timing correction with delta: $slicedt"
    # fi

    # # Flip anlge for look-locker readout
    # if [ ! -z $fa ]; then
    #     echo "--FA=$fa" >> $tempdir/basil_options.txt
    #     Log "Flip angle (look-locker readout): $fa"
    # fi

    # # Multi-band setup (if not set then this is ignored)
    # if [ ! -z $sliceband ]; then
    #     echo "--sliceband=$sliceband" >> $tempdir/basil_options.txt
    #     Log "Multi-band setup with number of slices per band: $slicedband"
    # fi

    # # Infer arterial component?
    # if [ -z $artoff ]; then
    #     basil_options=$basil_options"--inferart "
    #     Log "Infer arterial component"
    # fi
    # # fix the bolus duration?
    # if [ -z $fixbolus ]; then
    #     basil_options=$basil_options"--infertau "
    #     Log "Varaiable bolus duration"
    # else
    #     Log "Fixed bolus duration"
    # fi

    # if [ -z $fixbat ]; then
    #     Log "Variable arterial arrival time"
    #     Log "Setting prior/initial (tissue/gray matter) bolus arrival time to $bat"

    #     # Tissue BAT SD
    #     #defaults
    #     if [ -z $singleti ]; then
    #     # multi TI/PLD data, set a more liberal prior for tissue ATT since we should be able to determine from data
    #     # NB this leave the arterial BAT alone.
    #     batsd=1;
    #     fi

    #     # if required add batsd option for basil
    #     if [ ! -z $batsd ]; then
    #     Log "Setting std dev of the (tissue) BAT prior to $batsd"
    #     echo "--batsd=$batsd" >> $tempdir/basil_options.txt
    #     fi
    # else
    #     basil_options=$basil_options"--fixbat " 
    #     Log "Fixed arterial arrival time"
    #     Log "Setting arterial arrival time to $bat"
    # fi

    # # Noise specification
    # if [ -z $snr ]; then
    #     snr=10; #default SNR
    # fi
    # if [ $tpoints -eq 1 ]; then
    #     # only a single time point in data, will use informative noise prior
    # noiseprior=1
    # Log "Single volume: informative noise prior will be used"
    # fi

    # if [ ! -z $noiseprior ]; then
    #     # use an informative nosie prior
    #     if [ -z $noisesd ]; then
    #     Log "Using SNR of $snr to set noise std dev"
    #     # estimate signal magntiude
    #     fslmaths $tempdir/diffdata_mean -Tmax $tempdir/datamax
    #     brain_mag=`fslstats $tempdir/datamax -k $tempdir/mask -M`
    #     # this will correspond to whole brain CBF (roughly) - about 0.5 of GM
    #     noisesd=`echo "scale=2;sqrt( $brain_mag * 2 / $snr )" | bc`
    #     fi

    #     Log "Using a prior noise sd of: $noisesd"
    #     echo "--prior-noise-stddev=$noisesd" >> $tempdir/basil_options.txt
    # fi


    # # Exteneded options for BASIL
    # if [ ! -z $spatial ]; then
    # # if we are using spatial smoothing on CBF then we will also do the analysis in a single step
    #     echo "Instructing BASIL to use automated spatial smoothing"
    #     basil_options=$basil_options"--spatial "
    #     Log "Employing spatial VB"

    # fi

    # if [ ! -z $infert1 ]; then
    #     echo "Instructing BASIL to infer variable T1 values"
    #     basil_options=$basil_options"--infert1 "
    #     Log "Including T1 uncertainty"
    # fi


    # if [ ! -z $exch ]; then
    #     # use a specific exchange model in BASIL
    #     Log "Using exchange model: $exch"
    #     basil_options=$basil_options"--exch=$exch "
    # fi

    # if [ ! -z $disp ]; then
    #     # use a specific dispersion model in BASIL
    #     Log "Using dispersion model: $disp"
    #     basil_options=$basil_options"--disp=$disp "
    # fi

    # if [ ! -z $devel ]; then
    #     basil_options=$basil_options" --devel "
    # fi

    # if [ ! -z $model_options ]; then
    #     echo "Appending additional BASIL options to $tempdir/basil_options.txt"
    #     cat $model_options >> $tempdir/basil_options.txt
    # fi

    # Log "BASIL options ($tempdir/basil_options.txt):"
    # Log "----"
    # Log "`cat $tempdir/basil_options.txt`"
    # Log "----"
    # # -- end of main basil options setting

# def basil_init():
#     ##### Analyse data using BASIL
#     ### First analysis on whole data, normal perfusion image
#     echo "Calling BASIL on data - conventional perusion image"
#     initbasil="" # Can be used to pass an intital MVN into the main run of BASIL from an intial run (below)
#     if [ $repeats -gt 1 ] || [ $repeats -eq 0 ]; then
#         # do an initial analysis using the data averaged at each TI
#         # NB repeats=0 is a special case of variable number of repeats at each TI
#         Log "Initial run of BASIL on data where we have avareged all repeats at each TI"
#         datafile=$tempdir/diffdata_mean
#         mkdir $tempdir/init
#         cat $tempdir/basil_options_core.txt > $tempdir/init/basil_options.txt
#         echo "--repeats=1" >> $tempdir/init/basil_options.txt
#         echo "$tislist" >> $tempdir/init/basil_options.txt
#         Dobasil $datafile $tempdir/init
#         initbasil=$tempdir/init/finalMVN
#     fi

# def basil_main()
#     # main analysis using full data
#     datafile=$tempdir/diffdata
#     if [ $repeats -gt 0 ]; then
#         echo "--repeats=$repeats" >> $tempdir/basil_options.txt
#     else
#         # variable number of repeats at each TI - tell basil
#         echo "$rptslist" >> $tempdir/basil_options.txt
#     fi
#     echo "$tislist" >> $tempdir/basil_options.txt

#     Log "Main run of BASIL on ASL data"
#     Dobasil $datafile $tempdir $initbasil
#     ### End of: First analysis on whole data

# ### Registration (2/2)
# # Revisit the registration now that we have a pefusion image (with good GM/WM contrast) using BBR
# # use existing registration for initial aligment
# if [ $register -eq 1 ]; then
#     if [ $finalreg -eq 1 ]; then
# 	echo "Performing final registration"
# 	Log "Final registration"
# 	cp $tempdir/asl2struct.mat $outdir/native_space/asl2struct_init.mat #preserve the intial registration for future reference
# 	Registration $tempdir/ftiss "-m $tempdir/mask --tissseg $tempdir/tissseg --imat $tempdir/asl2struct.mat --finalonly"
#     fi
# fi
# ### End of: Registration (2/2)

# ### Partial Volume Estimates
# # Note we do this here since we have the final registration now which we need to transform PV estimates into ASL space
# if [ ! -z $fasthasrun ] && [ -z $pvgm ]; then
#     # PVE in ASL space from strcutural segmentation results
#     # invert the transformation matrix
#     convert_xfm -omat $tempdir/struct2asl.mat -inverse $tempdir/asl2struct.mat
    
#     # Gray matter - assume this will be PVE 1
#     applywarp --ref=$tempdir/asldata --in=$tempdir/pvgm_struct --out=$tempdir/pvgm_inasl --premat=$tempdir/struct2asl.mat --super --interp=spline --superlevel=4
#     # white matter  - assume this will be PVE 2
#     applywarp --ref=$tempdir/asldata --in=$tempdir/pvwm_struct --out=$tempdir/pvwm_inasl --premat=$tempdir/struct2asl.mat --super --interp=spline --superlevel=4
#     # threshold (upper and lower) the PVE to avoid artefacts of spline interpolation and also ignore very low PVE that could cause numerical issues.
#     fslmaths $tempdir/pvgm_inasl -thr 0.1 -min 1 $tempdir/pvgm_inasl
#     fslmaths $tempdir/pvwm_inasl -thr 0.1 -min 1 $tempdir/pvwm_inasl
#     pvexist=1
# fi

# if [ ! -z $pvgm ]; then
#     #using supplied PV images
# 	Log "Loading supplied PV images"
# 	if [ -z $pvwm ]; then
# 	    echo "ERROR: no WM PV image has been supplied"
# 	fi
# 	Log "PV GM is: $pvgm"
# 	fslmaths $pvgm -thr 0.1 -min 1 $tempdir/pvgm_inasl
# 	Log "PV WM is: $pvwm"
# 	fslmaths $pvwm -thr 0.1 -min 1 $tempdir/pvwm_inasl
# 	pvexist=1
# fi

# if [ ! -z $pvexist ]; then
#     # make some masks 
#     # these are currently used for masking after model fitting
#     fslmaths $tempdir/pvgm_inasl -thr 0.1 -bin $tempdir/gmmask
#     fslmaths $tempdir/pvwm_inasl -thr 0.1 -bin $tempdir/wmmask
#     # these are for calculating mean perfusion within tissue types
#     fslmaths $tempdir/pvgm_inasl -thr 0.8 -bin $tempdir/gmmask_pure
#     fslmaths $tempdir/pvwm_inasl -thr 0.9 -bin $tempdir/wmmask_pure
# fi
# ### End of: Partial Volume Estimates

# ### Calibration
# # Do calibration here becuase we do not need it before this point & if we are generating a CSF mask we have a better registration at this point
# if [ -z $t1tset ]; then
#     t1tset=1.3;
# fi
# Log "T1t (for calibration): $t1tset"

# # TR (for calibration image)
# if [ -z $tr ]; then
#     tr=3.2
# fi


# # calibration image gain
# if [ -z $cgain ]; then
#     cgain=1;
# fi

# # Calibration if reqd
# if [ -z $alpha ]; then
#         # based on the ASL white paper
#     if [ -z $casl ]; then
# 	alpha=0.98;
#     else
# 	alpha=0.85;
#     fi
# fi

# if [ -z $cmethod ]; then
# # default calibration method is 'voxelwise' unless we have CSF PV estimates or CSF mask has been supplied
#     if [ ! -z $fasthasrun ] || [ ! -z $csf ]; then
# 	cmethod=single
#     else
# 	cmethod=voxel
#     fi
# fi

# if [ ! -z $calib ]; then

#     # Single M0 value for calibration
#     if [ $cmethod = 'single' ]; then
# 	Log "Calibration is using a single M0 value with a CSF reference"
# 	if [ -z $csf ] && [ -z $fasthasrun ]; then
# 	    echo "ERROR: Provide either a structural image or CSF mask for calibration when using --cmethod=single"
# 	    exit 1
# 	fi
#        # calcualte M0a from CSF
# 	Calibration
# 	Mo=`cat $outdir/calib/M0.txt`

#     # Voxelwise M0 values for calibration
#     elif [ $cmethod = 'voxel' ]; then
# 	Log "Calibration is voxelwise"
# 	mkdir $outdir/calib
#         # copy over the calibration image and apply the cgain setting - this increases the magntiude of M0 to match that of the ASL data (acquired with a higher gain - cgain>=1 normally)
# 	fslmaths $calib -mul $cgain $outdir/calib/M0
# 	Mo=$outdir/calib/M0 
# 	if [ 1 -eq `echo "$tr < 5" | bc`  ]; then
# 	 # correct the M0 image for short TR using the equation from the white paper
# 	    Log "Correcting the calibration (M0) image for short TR (using T1 of tissue $t1tset)"
# 	    ccorr=`echo "1 / (1 - e(- $tr / $t1tset) )" | bc -l`
# 	    fslmaths $Mo -mul $ccorr $Mo
# 	fi

# 	#inlcude partiition co-effcient in M0 image to convert from M0 tissue to M0 arterial
# 	fslmaths $Mo -div 0.9 $Mo

# 	if [ ! -z $edgecorr ]; then
# 	    # correct for (partial volume) edge effects
# 	    # median smoothing and erosion
# 	    fslmaths $Mo -fmedian -mas $tempdir/mask -ero $tempdir/calib_ero
# 	    # extrapolation to match mask
# 	    asl_file --data=$tempdir/calib_ero --ntis=1 --mask=$tempdir/mask --extrapolate --neighbour=5 --out=$Mo
# 	fi
#     else
# 	echo "Error unrecognised calibration method: $cmethod, (use single or voxel)"
#     fi 
# elif [ ! -z $M0 ]; then
#     # An M0 value has been supplied, use this
#     cmethod=single # we are in 'single' mode as a single value has been supplied
#     Mo=$M0
#     echo "M0: $Mo"
#     Log "Using supplied M0 value: $Mo"
# fi

# ### End of: Calibration


# ### Output main BASIL results
# # Note we do this here, as we have the registration done and masks created and calibration complete
# Dooutput

# # save the mask used to the (native space) output directory
# imcp $tempdir/mask $outdir/native_space/mask
# ### End of: Output main BASIL results

# ### Partial Volume Correction BASIL
# if [ ! -z $pvcorr ]; then
#     if [ -f $tempdir/struc_bet.* ]; then
# 	# redo the mask now that we have a better registration - as this might matter for PV correction
# 	# NB we dont use the PVE here since we dont (necessarily) want to exclude the ventricles from the mask as this has implications for the spatial priors
# 	fslmaths $tempdir/struc_bet -bin $tempdir/struc_bet_mask
# 	flirt -in $tempdir/struc_bet_mask -ref $regfrom -applyxfm -init $tempdir/struct2asl.mat -out $tempdir/mask -interp trilinear
# 	fslmaths $tempdir/mask -thr 0.25 -bin -fillh $tempdir/mask
# 	fslcpgeom $regfrom $tempdir/mask
# 	imcp $tempdir/mask $outdir/native_space/mask_pvcorr # copy new mask to output directory - indicate that it was used for PV correction analysis
#     fi
    
#     # intructions for BASIL
#     basil_options=$basil_options" --pgm $tempdir/pvgm_inasl --pwm $tempdir/pvwm_inasl "
#     mkdir $tempdir/pvcorr
#     cp $tempdir/basil_options.txt $tempdir/pvcorr/basil_options.txt #Dobasil expects the options file to be in the subdirectory
#     # Run BASIL
#     Dobasil $datafile $tempdir/pvcorr

#     imcp $tempdir/pvcorr/finalMVN $tempdir/finalMVN #just in case we are about to do a epochwise analysis

#     #output the results
#     Dooutput pvcorr
# fi
# ### End of: Partial Volume Correction

# ### Epoch BASIL
# if [ ! -z $epoch ]; then
#     # epochwise analysis
#     echo "Epochwise analysis"

#     #genereate a list of epochs
#     currdir=`pwd`
#     cd $tempdir
#     epochlist=`imglob epoch*`
#     cd $currdir

#     ecount=0
#     for e in $epochlist; do
# 	Log "Processing epoch: $e"
# 	etislist=""
#         # deal with the TIs
# 	for ((ei=0; ei<$elen; ei++)); do
# 	    ethis=`expr $ecount \* $eadv + $ei`
# 	    #echo $ethis
# 	    eidx=`expr $ei + 1`
# 	    #echo $ei
# 	    #echo ${alltis[$ethis]}
# 	    etislist=$etislist" --ti${eidx}=${alltis[$ethis]}"
# 	done
# 	Log "TIs for this epoch: "
# 	Log $etislist

# 	mkdir $tempdir/$e
# 	cp $tempdir/basil_options_core.txt $tempdir/$e/basil_options.txt # get the 'core' options and make a new basil_options file jsut for this TI
# 	echo "--repeats=1" >> $tempdir/$e/basil_options.txt #for epochs we specify all the TIs explicitly
# 	echo $etislist >>  $tempdir/$e/basil_options.txt #these are the basil options for this epoch

# 	fast=2 #we now switch BASIL to fast level '2' - this means it will only do analysis in a single step from here on in, but we will use our existing analysis for initialisation.
	
# 	Dobasil $tempdir/$e $tempdir/$e $tempdir/basil/finalMVN # init with results of first basil run

#         #output 
# 	Log "Saving results from epoch: $e"
# 	Dooutput $e

# 	ecount=`expr $ecount + 1`
#     done
# fi
# ### End of: Epoch BASIL




# #OUTPUTS
# # Setup option outputs - anything that would be common to all epochs
# # note that we now do directory creation right at the start
# #if [ ! -z $nativeout ]; then
# #fi
# if [ ! -z $structout ]; then
#     #cp $tempdir/asl2struct.mat $outdir/struct_space/asl2struct.mat
#     cp $tempdir/asl2struct.mat $outdir/native_space/asl2struct.mat #also provide the transformation matrix for reference
# fi
# #if [ ! -z $advout ]; then
# #fi

# #if [ -z $epoch ]; then
# # normal single analysis of data
# #Dooutput

# ##if [ ! -z $epoch ]; then
# # epochwise analysis
# #    for e in $epochlist; do
# #	Log "Saving results from epoch: $e"
# #	Dooutput $e
#  #   done
# #fi




# if [ ! -z $pvcorr ]; then
# # copy PVE in ASL space to output directory
# imcp $tempdir/pvgm_inasl $outdir/native_space/pvgm_inasl
# imcp $tempdir/pvwm_inasl $outdir/native_space/pvwm_inasl
# fi

# if [ ! -z $pvexist ]; then
#     # copy PV masks to output directory
#     imcp $tempdir/gmmask $outdir/native_space/gm_mask
#     imcp $tempdir/wmmask $outdir/native_space/wm_mask
#     imcp $tempdir/gmmask_pure $outdir/native_space/gm_roi
#     imcp $tempdir/wmmask_pure $outdir/native_space/wm_roi
# fi

# # clearup
# if [ ! -z $debug ]; then
#     mv $tempdir $outdir
# else
#     rm -r $tempdir
# fi

