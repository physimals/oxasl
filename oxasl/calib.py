#!/bin/env python
"""
ASL_CALIB: Calibration for ASL data

Michael Chappell & Brad MacIntosh, FMRIB Image Analysis & Physics Groups

Copyright (c) 2008-2013 Univerisity of Oxford
"""

import sys
import math
import traceback
from optparse import OptionParser

import numpy as np

from . import __version__, AslOptionGroup, fsl

def main():
    usage = """ASL_CALIB
    Calibration for ASL data

    asl_pycalib -i <perfusion image> -c <calibration image> --method <voxelwise|refregion> -o <output filename> [options]"""

    debug = False
    try:
        p = OptionParser(usage=usage, version=__version__)
        add_calib_options(p)
        options, _ = p.parse_args()
        
        # Convert to dictionary for easier handling
        options = vars(options)
        debug = options.pop("debug", False)

        perf_img = fsl.Image(options.pop("perf", None), "Perfusion")
        perf_img.summary()
        calib_img = fsl.Image(options.pop("calib", None), "Calibration")
        calib_img.summary()
        mask = options.pop("mask", None)
        if mask is not None:
            mask_img = fsl.Image(mask, "Mask")
        else:
            mask_img = fsl.Image("mask", role="Mask", data=np.ones(perf_img.shape))
            
        output_name = options.pop("output", None)
        if output_name is None:
            output_name = "%s_calib" % perf_img.iname
        print(output_name)

        calibrated_img = calib(perf_img, calib_img, output_name, **options)
        calibrated_img.summary()
        calibrated_img.save()
    except Exception as e:
        sys.stderr.write("ERROR: " + str(e) + "\n")
        if debug:
            traceback.print_exc()
        sys.exit(1)

def add_calib_options(parser, ignore=()):
    g = AslOptionGroup(parser, "Calibration", ignore=ignore)
    g.add_option("-c", dest="calib", help="Calibration image")
    g.add_option("-i", dest="perf", help="Perfusion image for calibration")
    g.add_option("-m", dest="mask", help="Mask in perfusion/calibration image space")
    g.add_option("-o", dest="output", help="Output filename for calibrated image (defaut=<input_filename>_calib")
    g.add_option("--method", dest="method", help="Calibration method: voxelwise or refregion")
    g.add_option("--debug", dest="debug", action="store_true", default=False, help="Debug mode")
    #g.add_option("-s", dest="struc", help="Structural image")
    #g.add_option("-t", dest="trans", help="ASL->Structural transformation matrix")
    #g.add_option("--mode", dest="mode", help="Calibration mode (longtr or satrevoc)", default="longtr")
    #g.add_option("--tissref", dest="tissref", help="Tissue reference type (csf, wm, gm or none)", default="csf")
    #g.add_option("--te", dest="te", help="Sequence TE", type=float, default=0.0)
    #g.add_option("--Mo", dest="mo", help="Save calculated M0 value to specified file")
    #g.add_option("--om", dest="om", help="Save CSF mask to specified file")

    parser.add_option_group(g)

    # echo ""
    # echo " Extended options (all optional):"
    # echo " -m         : Specify reference mask in calibration image space"
    # echo "              - strucutral image & transformation matrix are not required"
    # echo " --bmask    : Brain mask (in ASL data space) for sensitivity or tissue T1 estimation"
    # echo " --t2star   : Correct with T2* rather than T2"
    # echo "               (this alters the default values specified below to the T2* values)"
    # echo " --t1r      : T1 of reference tissue (defaults: csf 4.3, gm 1.3, wm 1.0 s) "
    # echo " --t2r      : T2(*) of reference tissue (defaults T2/T2*: csf 750/400, gm 100/60,  wm 50/50  ms)"
    # echo " --t2b      : T2(*) of blood (default T2/T2*: 150/50 ms)"
    # echo " --pc       : Partition co-efficient (defaults csf 1.15, gm 0.98,  wm 0.82)"
    # echo " --alpha    : Specify inversion efficiency - only applied to final CBF image calculation"
    # echo ""
    # echo " CSF masking options (only for --tissref csf)"
    # echo "  By default asl_calib extracts CSF from the structural image by segmentation and"
    # echo "  this is then masked using the ventricles in MNI152 space."
    # echo " --csfmaskingoff : turns off the ventricle masking, reference is based on segmentation only."
    # echo "  Registration between structural image and MNI152 is done automatically unless:"
    # echo "  --str2std  : Structural to MNI152 linear registration (.mat)"
    # echo "  --warp     : Structural to MNI152 non-linear registration (warp)"
    # echo ""
    # echo "MODES:"
    # echo "> longtr  Calibration image is a control image with a long TR."
    # echo "  {--tr}     : TR used in calibration sequence - {default: 3.2s}"
    # echo "  {--cgain}  : Relative gain between calibration and ASL data - {default: 1}"
    # echo ""
    # echo "> satrecov  Calibration image is a sequnce of control images at various TIs"
    # echo "            M0 is to be determined from a saturation recovery"
    # echo "            T1 of tissue (and FA correction) images are also calcualted"
    # echo " --tis       : comma separated list of inversion times, e.g. --tis 0.2,0.4,0.6"
    # echo " {--fa}      : Flip angle (in degrees) for Look-Locker readouts"
    # echo "   >> Look-Locker flip angle correction - to perform this provide:"
    # echo " {--lfa}     : Lower flip angle (in degrees) for dual FA calibration"
    # echo " {--nphases} : Number of phases (repetitions) of higher FA"
    # echo " {--fixa}    : Fix the saturation efficiency to 100% (useful if you have a low number of samples)"
    # echo ""
    # echo "Coil sensitivity correction:"
    # echo " Calculate and apply a voxel-wise correction for coil sensitivity"
    # echo " > using bias field from structural image (default)"
    # echo " {--osen}    : save sensitivity image to specified file."
    # echo " > using existing sensitivity image:"
    # echo "  --isen     : input coil sensitivity image"
    # echo " > using reference images (collected using same parameters):"
    # echo "  --cref     : Reference image from coil with minimal variation e.g. body."
    # echo "  {--cact}   : Image from coil used for actual ASL acquisition"
    # echo "               {default: calibration image - only in longtr mode}"
    # echo ""

def calib(perf_data, calib_data, output_name, method, multiplier=1.0, var=False, log=sys.stdout, **kwargs):
    """
    Do voxelwise calibration

    FIXME edge correction not applied

    :param data: fsl.Image containing data to calibrate
    :param m0: fsl.Image containing voxelwise m0 map
    :param alpha: Inversion efficiency
    :param var: If True, assume data represents variance rather than value
    """
    if method == "voxelwise":
        m0 = get_m0_voxelwise(calib_data, log=log, **kwargs)
    elif method == "refregion":
        m0 = ge_m0_refregion(calib_data, log=log, **kwargs)
    else:
        raise ValueError("Unknown calibration method: %s" % method)

    if var:
        m0 = m0**2
        multiplier = multiplier**2

    m0[m0 == 0] = 1

    calibrated = perf_data.data() / m0 
    calibrated *= multiplier
    return perf_data.derived(calibrated, name=output_name)

def get_m0_voxelwise(calib_data, gain=1.0, alpha=1.0, tr=None, t1=None, log=sys.stdout):
    """
    Calculate M0 value using voxelwise calibration

    FIXME edge correction not applied

    :param calib_data: fsl.Image containing voxelwise m0 map
    :param gain: Calibration gain
    :param alpha: Inversion efficiency
    :param tr: Sequence TR (s)
    """
    log.write("Doing voxelwise calibration\n")
    
    # Calculate M0 value
    m0 = calib_data.data() * alpha * gain

    if tr is not None and tr < 5:
        if t1 is not None:
    	    # correct the M0 image for short TR using the equation from the white paper
            log.write("Correcting the calibration (M0) image for short TR (using T1 of tissue %f)\n" % t1)
            ccorr = 1 / (1 - math.exp(-tr / t1))
            m0 *= ccorr
        else:
            log.write("WARNING: tr < 5 (%f) but reference tissue T1 not provided so cannot apply correction\n" % tr)

	# Include partiition co-effcient in M0 image to convert from M0 tissue to M0 arterial
    m0 /= 0.9

	#if edgecorr:
	    # Correct for (partial volume) edge effects
	    # median smoothing and erosion
	    #fslmaths $Mo -fmedian -mas $tempdir/mask -ero $tempdir/calib_ero
	    # extrapolation to match mask
	    #asl_file --data=$tempdir/calib_ero --ntis=1 --mask=$tempdir/mask --extrapolate --neighbour=5 --out=$Mo

    return m0

def ge_m0_refregion(calib_data, ref_mask, log=sys.stdout, **kwargs):
    """
    Do reference region calibration

    FIXME this is not yet complete

    :param calib_data: Calibration image
    :param ref_mask: Reference region mask image
    """
    log.write("Doing reference region calibration\n")
    
    mode = kwargs.get("mode", "longtr")
    tr = kwargs.get("tr", 3.2)
    te = kwargs.get("te", 0)
    taq = kwargs.get("taq", 0)
    log.write("Using TE=%f" % te)
        
    # Check if we have a sensitivity map
    sens = kwargs.get("sens", None)
    if sens is not None:
        log.write("Using sensitivity image: %s\n" % sens.iname)

    # Constants
    T2b = 150 # lu et a. 2012 MRM 67:42-49 have 154ms at 3T during normoxia

    # Parameters for reference tissue type: T1, T2, partition coeffs, FAST seg ID
    # Partition coeffs based on Herscovitch and Raichle 1985 with a blood water density of 0.87
    params = {
        "csf" : [4.5, 750, 1.15, 0],
        "wm" : [1.0, 50, 0.82, 2],
        "gm" : [1.3, 100, 0.98, 1],
    }
    if kwargs.get("t2star", False):
        # We need to correct for T2* not T2 so change the defaults
        # NB these will still be overridden by user-specified values 
        params["csf"][1] = 400
        params["gm"][1] = 60 # from Foucher 2011 JMRI 34:785-790
        params["wm"][1] = 50 # ditto
        T2b = 50 #from Petersen 2006 MRM 55(2):219-232 see discussion

    tissref = kwargs.get("tissref", "csf")
    log.write("Using tissue reference type: %s\n" % tissref)
         
    t1r_img, t2r_img = False, False
    if tissref in params:
        t1r, t2r, pc, fast_id = params[tissref]
    else:
        try:
            t1r, t2r, pc, fast_id = None, None, None, int(tissref)
        except:
            raise ValueError("Invalid tissue reference type: %s" % tissref)

    # Command line override of default T1 and T2
    if "t1r" in kwargs:
        t1r = kwargs.get("t1r", None)
        if isinstance(t1r, fsl.Image):
            log.write("Using T1 image for reference region: %s" % t1r.iname)
            t1r_img = True
        elif t1r is not None:
            log.write("Using user-specified T1r value: %f" % t1r)

    if "t2r" in kwargs:
        t2r = kwargs.get("t2r", None)
        if isinstance(t2r, fsl.Image):
            log.write("Using T2 image for reference region: %s" % t2r.iname)
            t2r_img = True
        elif t2r is not None:
            log.write("Using user-specified T2r value: %f" % t2r)

    if "t2b" in kwargs:
        t2b = kwargs["t2b"]
        log.write("Using user-specified T2b value: %f" % t2b)
        
    if "pc" in kwargs:
        pc = kwargs["pc"]
        log.write("Using user-specified partition coefficient: %f" % pc)

    if t1r is None:
        raise ValueError("T1 for reference tissue has not been set")
    if t2r is None:
        if te != 0:
            raise ValueError("T2 for reference tissue has not been set")
        else:
            t2r = 1.0
    if pc is None:
        raise ValueError("Partition coefficient  for reference tissue has not been set")

    log.write("T1r: %f; T2r: %f; T2b: %f; Part co-eff: %f\n" % (t1r, t2r, t2b, pc))
        
    # # sort out the M0 calib brain_mask
    # if [ -z $bmask ]; then
    #     echo "Creating brain mask from calibration image" >> $log
    #     #make a brain mask
    #     # take the mean
    #     fslmaths $calib -Tmean $temp_calib/calib_mean
    #     # bet
    #     bet $temp_calib/calib_mean $temp_calib/calib_mean -m #calib_mean_mask is the brain mask for the calib image
    #     bmask=$temp_calib/calib_mean_mask
    # fi

    # ### Sensitivity image calculation (if reqd)
    # if [ ! -z $crefim ]; then
    #     echo "Calculate sensitivity image" >> $log
    #     senson=1
    #     # take the mean (and mask with the mask from the main calib image)
    #     fslmaths $crefim -Tmean -mas $bmask $temp_calib/crefim
        
    #     # take the ratio to give the sensitivity image
    #     if [ -z $cactim ]; then
    #         # if the cact image has not been supplied then use the mean of the calib image
    #     if [ ! $mode = longtr ]; then
    #         echo "ERROR: You must supply an image from the actual coil used for ASL acquisition using --cact (unless you use longtr mode)"
    #         exit 1
    #     fi
    #     fslmaths $calib -Tmean $temp_calib/cactim
    #     fi 
    #     fslmaths $temp_calib/cactim -div $temp_calib/crefim -mas $bmask $temp_calib/sens
    # fi

    # if [ $tissref = "none" ]; then
    # # whole brain M0
    # # in this case use the brain mask
    #     imcp $bmask $temp_calib/refmask
    #     maskflag=1
    #     echo "Brain mask is being used as the reference tissue (beware!)" >> $log
    # fi

    # if [ -z $maskflag ]; then

    # # make brain mask from structural
    #     fslmaths $struc -bin $temp_calib/mask

    #     if [ -z $refpve ]; then
    #         # auto create tissue reference mask
    #     echo "FAST called to determine a reference tissue mask" >> $log
        
    #         # segment structural image
    #     fast -b -o $temp_calib/seg -p $struc
    #     fasthasrun=1;
    #     imcp  $temp_calib/seg_pve_$fastpve $temp_calib/refpve

    #     else
    #     # user supplied PV estimate for reference tissue
    #     echo "Using input reference PVE: $refpve" >> $log
    #     imcp $refpve $temp_calib/refpve
    #     fi

    #     if [ $tissref = "csf" ] & [ -z $csfmaskingoff ]; then
    #     echo "Ventricle selection" >> $log
    #     stdmaskfnirt=1  # by deafult now we do FNRIT transformation of ventricle mask

    # cut down brain mask so that it only covers middle of brain
    # sort out the roi
    #	  xsize=`fslinfo $struc | grep "^dim1" | sed 's:dim1[ ]*::'`
    #	  ysize=`fslinfo $struc | grep "^dim2" | sed 's:dim2[ ]*::'`
    #	  zsize=`fslinfo $struc | grep "^dim3" | sed 's:dim3[ ]*::'`
    #	  roisize="0.3";
    #     echo "$xsize $ysize $zsize $roisize"
    #	  delx=`echo "v = $roisize * $xsize; v /= 1; v" | bc`
    #	  xmin=`echo "v = 0.5 * $xsize - $delx / 2; v /= 1; v" | bc`
    #	  dely=`echo "v = $roisize * $ysize; v /= 1; v" | bc`
    #	  ymin=`echo "v = 0.5 * $ysize - $dely / 2; v /= 1; v" | bc`
    #	  delz=`echo "v = $roisize * $zsize; v /= 1; v" | bc`
    #	  zmin=`echo "v = 0.5 * $zsize - $delz / 2 + 0.1; v /= 1; v" | bc`
    #     echo "$xmin $delx $ymin $dely $zmin $delz" >> $log
    #	  fslmaths $temp_calib/mask -roi $xmin $delx $ymin $dely $zmin $delz 0 1 $temp_calib/mask

    #     # select ventricles based on standard space atlas
	# fslroi $FSLDIR/data/atlases/HarvardOxford/HarvardOxford-sub-prob-2mm $temp_calib/LVentricle 2 1
	# fslroi $FSLDIR/data/atlases/HarvardOxford/HarvardOxford-sub-prob-2mm $temp_calib/RVentricle 13 1
	# fslmaths $temp_calib/LVentricle -add $temp_calib/RVentricle -thr 0.1 -bin -ero $temp_calib/VentricleMask
	
    #     # register structural image to std space using FLIRT
	# if [ -z $str2std ] & [ -z $warp ]; then
	#     echo "Registering structural image to standard space using FLIRT" >> $log
	#     flirt -in $struc -ref $FSLDIR/data/standard/MNI152_T1_2mm_brain.nii.gz -omat $temp_calib/struc2std.mat
	# else
	#     if [ -f $str2std ]; then
	# 	cp $str2std $temp_calib/struc2std.mat
	#     else
	# 	echo "Error $str2std not found"
	# 	exit 1
	#     fi
	# fi
	
	# if [ -z $stdmaskfnirt ] & [ -z $warp ]; then
    #     # do std space masking with FLIRT registration
	    
	#     convert_xfm -omat $temp_calib/std2struc.mat -inverse $temp_calib/struc2std.mat
	#     flirt -in $temp_calib/VentricleMask -applyxfm -init $temp_calib/std2struc.mat -ref $temp_calib/refpve -out $temp_calib/mask
	# else
    #         # do std space masking using FNIRT registration
    #         # Register the structural image to standard space using FNIRT
	#     if [ -z $warp ]; then
	# 	echo "Registering structural image to standard space using FNIRT" >> $log
	# 	fnirt --in=$struc --aff=$temp_calib/struc2std.mat --config=T1_2_MNI152_2mm.cnf --cout=$temp_calib/coef_struc2MNI152
	#     else
	# 	if [ -f $warp ]; then
	# 	    cp $warp $temp_calib/coef_struc2MNI152
	# 	else
	# 	    echo "Error: $warp not found"
	# 	    exit 1
	# 	fi
	#     fi
    #         # Calculate the inverse warp using INVWARP and apply to standard space ventricle mask
	#     invwarp --ref=$struc --warp=$temp_calib/coef_struc2MNI152 --out=$temp_calib/coef_MNI2struc
	#     applywarp --ref=$temp_calib/refpve --in=$temp_calib/VentricleMask --warp=$temp_calib/coef_MNI2struc --out=$temp_calib/mask --interp=nn
	#     # using refpve as reference for applywarp as this is the image we apply the mask to, nit the structural
	# fi
	
	# echo "Masking FAST output with standard space derrived ventricle mask" >> $log
	# fslmaths $temp_calib/refpve -mas $temp_calib/mask $temp_calib/refpve
	# if [ ! -z $debug ]; then
	#     imcp $temp_calib/refmask $temp_calib/refmask_high_unthresh
	# fi
    # fi
    
    # echo "Transforming tissue reference mask into perfusion space" >> $log
    # #transform mask into perfusion space
    # convert_xfm -omat $temp_calib/high2low.mat -inverse $trans
    # #flirt -in $temp_calib/csf -applyxfm -init $temp_calib/high2low.mat -ref $calib -out $temp_calib/csf
    # # new conversion using applywarp, supersmapling and integration
    # applywarp --ref=$calib --in=$temp_calib/refpve --out=$temp_calib/refmask --premat=$temp_calib/high2low.mat --super --interp=spline --superlevel=4
    # if [ ! -z $debug ]; then
	# imcp $temp_calib/refmask $temp_calib/refmask_low_unthresh
    # fi
    
    # if [ ! -z $fasthasrun ] && [ -z $senson ]; then
    #     # also extract the bias field and convert to sensitivity image (as long as we have already been supplied by a sensivity iamge or reference)
	# applywarp --ref=$calib --in=$temp_calib/seg_bias --out=$temp_calib/biasfield --premat=$temp_calib/high2low.mat --super --interp=spline --superlevel=4
	# fslmaths $temp_calib/biasfield -recip $temp_calib/sens
	# senson=1
	# echo "Using bias field from structural image for sensitivity correction" >> $log
    # fi

    # threshold reference mask if it is not already binary
    #ref_mask = ref_mask > 0.9
    
    # Use supplied tissue reference mask
    log.write("Using supplied reference tissue mask: %s\n" % ref_mask.iname)
    
    mask = ref_mask.data()
    nonzero = np.count_nonzero(mask)
    log.write("Number of voxels in tissue reference mask: %i" % nonzero)
    if nonzero < 1:
        raise ValueError("Reference mask does not contain any unmasked voxels")

    log.write("MODE: %s\n" % mode)
    log.write("Calibration gain: %f" % gain)
    if mode == "longtr":
        # Calibration data is a long TR acquisition - all we need to do here is take the mean 
        if mask.ndim == 4:
            mask = np.mean(axis=-1)

        # if [ ! -z $senson ]; then
        #     echo "Apply sensitivity image" >> $log
        #     # apply sensitivity map to calibration image
        #     fslmaths $temp_calib/calib -div $temp_calib/sens $temp_calib/calib
        # fi

        #mask M0 map with tissue reference
        #calib_data[mask == 0] = 0
        
        # calcualte T1 of reference region (if a T1 image has been supplied)
        if t1r_img:
            t1r_data = t1r_img.data()
            t1r = np.mean(t1r_data)
            log.write("Calculated T1 of reference tissue: %f" % t1r)

        # calcualte T2 of reference region (if a T2 image has been supplied)
        if t2r_img:
            t2r_data = t2r_img.data()
            t2r = np.mean(t2r_data)
            log.write("Calculated T1 of reference tissue: %f" % t2r)

        # calculate M0_ref value
        m0 = np.mean(calib_data[mask > 0])
        m0 = m0 / (1 - math.exp(- (tr - taq) / t1r) )
        log.write("Mz of reference tissue: %f\n" % m0)
        
    elif mode == "satrecov":
        # Calibration image is control images and we want to do a saturation recovery fit
        # NB only do the fit in the CSF mask
        options = {
            "model" : "satrecov",
            "noise" : "white",
            "method" : "vb",
            "t1" : t1r,
        }
        
        #deal with TIs
        tis = kwargs.get("tis", [])
        log.write("TIs: %s" % str(tis))
        for idx, ti in enumerate(tis):
            options["ti%i" % (idx+1)] = ti

        # Extra options for Look Locker
        if "fa" in kwargs:
            options["FA"] = kwargs["fa"]
        if "nphases" in kwargs:
            options["phases"] = kwargs["nphases"]
        if "lfa" in kwargs:
            options["LFA"] = kwargs["lfa"]
        
        # Extra sat recovery options
        if kwargs.get("fixa", False):
            options["fixa"] = ""

        # do fabber within the tissue reference mask with a sensible T1 prior mean
        # if [ ! -z $senson ]; then
        # echo "Apply sensitivity image to data for reference tisse M0 estimation" >> $log
        #     # apply sensitivity map to calibration image - ONLY for the reference tissue calculations
        # fslmaths $calib -div $temp_calib/sens $temp_calib/calib_senscorr
        # else
        # # no sensitivity correction required, but copy image over ready for next command
        # imcp $calib $temp_calib/calib_senscorr
        # fi
        log("Running FABBER within reference tissue mask\n")
        wsp.fabber(calib_img, ref_mask, options)

        # calculate M0 value
        # Moval=`fslstats $temp_calib/satrecov/mean_M0t -M` # this is M0 of CSF at the TE of the sequence
        # echo "M0 of reference tissue: $Moval" >> $log

        # if [ ! -z $outdir ]; then
        #     # save useful results to specified output directory
        # imcp $temp_calib/satrecov/mean_T1t $outdir/T1_ref
        # imcp $temp_calib/satrecov/mean_M0t $outdir/M0_ref
        # fi

        # # do fabber again within whole brain to get estimated T1 of tissue and FA correction (if LL)
        # # (note that we do not apply sensitivity correction to the data here - thius is 'built-into' the M0t map)
        # echo "FABBER (again) within whole brain mask" >> $log
        # if [ ! -z $outdir ]; then #NB we only bother with this if we have an output directory to put the results in
        # $fabber --data=$calib --mask=$bmask --output=$temp_calib/satrecovT --data-order=singlefile --model=satrecov --noise=white --method=vb $tislist $llopts $sropts 


        #     # save useful results to specified output directory
        #     imcp $temp_calib/satrecovT/mean_T1t $outdir/T1t
        #     imcp $temp_calib/satrecovT/mean_M0t $outdir/M0t
        #     if [ ! -z $lfa ]; then
        #     imcp $temp_calib/satrecovT/mean_g $outdir/facorr
        #     fi
        # fi

    # use equation to get the M0 value that is needed
    m0 = m0 / math.exp(- te / t2r) #  T2 correction
    m0 = m0 * gain / pc
    m0 = m0 * math.exp(- te / t2b)
    log.write("M0: %f\n" % m0)

    # Apply calibration to input image
    # if [ ! -z $senson ]; then
    # # apply sensitivity image
    # fslmaths $infile -div $temp_calib/sens $temp_calib/infile
    # else
    # imcp $infile $temp_calib/infile
    # fi

    if alpha != 1.0:
        log.write("Applying inversion efficiency of: %f\n")
        # apply the inversion efficiency supplied to M0 prior to final calculation
        m0 = m0 * alpha
    
    return m0

if __name__ == "__main__":
    main()
