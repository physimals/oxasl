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
import scipy.ndimage

from . import __version__, AslOptionGroup, fslwrap as fsl

def main():
    usage = """ASL_CALIB
    Calibration for ASL data

    asl_pycalib -i <perfusion image> -c <calibration image> --calib-method <voxelwise|refregion> -o <output filename> [options]
    """

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
        brain_mask = options.pop("brain_mask", None)
        if brain_mask is not None:
            brain_mask_img = fsl.Image(brain_mask, "Mask")
        else:
            brain_mask_img = fsl.Image("brain_mask", role="Mask", data=np.ones(perf_img.shape))
            
        output_name = options.pop("output", None)
        if output_name is None:
            output_name = "%s_calib" % perf_img.iname

        calibrated_img = calib(perf_img, calib_img, output_name, brain_mask=brain_mask_img, **options)
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
    g.add_option("-i", dest="perf", help="Perfusion image for calibration, in same image space as calibration image")
    g.add_option("-m", "--brain-mask", dest="brain_mask", help="brain mask in perfusion/calibration image space")
    g.add_option("-o", dest="output", help="Output filename for calibrated image (defaut=<input_filename>_calib")
    g.add_option("--calib-method", dest="calib_method", help="Calibration method: voxelwise or refregion")
    g.add_option("--alpha", help="Inversion efficiency", type=float, default=1.0)
    g.add_option("--cgain", dest="gain", help="Relative gain between calibration and ASL data", type=float, default=1.0)
    g.add_option("--tr", help="TR used in calibration sequence (s)", type=float, default=3.2)
    g.add_option("--debug", dest="debug", action="store_true", default=False, help="Debug mode")
    parser.add_option_group(g)

    g = AslOptionGroup(parser, "Voxelwise calibration", ignore=ignore)
    g.add_option("--pct", help="Tissue/arterial partition coefficiant", type=float, default=0.9)
    g.add_option("--t1t", help="T1 of tissue (s)", type=float, default=1.3)
    parser.add_option_group(g)

    g = AslOptionGroup(parser, "Reference region calibration", ignore=ignore)
    g.add_option("--mode", help="Calibration mode (longtr or satrevoc)", default="longtr")
    g.add_option("--tissref", help="Tissue reference type (csf, wm, gm or none)", default="csf")
    g.add_option("--te", help="Sequence TE", type=float, default=0.0)
    g.add_option("--t1r", help="T1 of reference tissue (defaults: csf 4.3, gm 1.3, wm 1.0 s)", type=float, default=None)
    g.add_option("--t2r", help="T2 of reference tissue (defaults T2/T2*: csf 750/400, gm 100/60,  wm 50/50  ms)", type=float, default=None)
    g.add_option("--t2b", help="T2(*) of blood (default T2/T2*: 150/50 ms)", type=float, default=None)
    g.add_option("--refmask", dest="ref_mask", help="Reference tissue mask in perfusion/calibration image space")
    g.add_option("--t2star", action="store_true", default=False, help="Correct with T2* rather than T2 (alters the default T2 values)")
    g.add_option("--pcr", help="Reference tissue partition coefficiant (defaults csf 1.15, gm 0.98,  wm 0.82)", type=float, default=None)
    #g.add_option("--Mo", dest="mo", help="Save calculated M0 value to specified file")
    #g.add_option("--om", dest="om", help="Save CSF mask to specified file")
    parser.add_option_group(g)
    
    g = AslOptionGroup(parser, "longtr mode (calibration image is a control image with a long TR)", ignore=ignore)
    parser.add_option_group(g)

    g = AslOptionGroup(parser, "satrecov mode (calibration image is a sequnce of control images at various TIs)", ignore=ignore)
    g.add_option("--tis", help="Comma separated list of inversion times, e.g. --tis 0.2,0.4,0.6")
    g.add_option("--fa", help="Flip angle (in degrees) for Look-Locker readouts", type=float)
    g.add_option("--lfa", help="Lower flip angle (in degrees) for dual FA calibration", type=float)
    g.add_option("--nphases", help="Number of phases (repetitions) of higher FA", type=int)
    g.add_option("--fixa", action="store_true", default=False, help="Fix the saturation efficiency to 100% (useful if you have a low number of samples)")
    parser.add_option_group(g)

    g = AslOptionGroup(parser, "Coil sensitivity correction, either using existing sensitivity image or reference images collected using same parameters", ignore=ignore)
    g.add_option("--isen", help="Input coil sensitivity image")
    g.add_option("--cref", help="Reference image from coil with minimal variation e.g. body.")
    g.add_option("--cact", help="Image from coil used for actual ASL acquisition (default: calibration image - only in longtr mode)")
    parser.add_option_group(g)

    # echo " CSF masking options (only for --tissref csf)"
    # echo "  By default asl_calib extracts CSF from the structural image by segmentation and"
    # echo "  this is then masked using the ventricles in MNI152 space."
    #g.add_option("-s", dest="struc", help="Structural image")
    #g.add_option("-t", dest="trans", help="ASL->Structural transformation matrix")
    # echo " --csfmaskingoff : turns off the ventricle masking, reference is based on segmentation only."
    # echo "  Registration between structural image and MNI152 is done automatically unless:"
    # echo "  --str2std  : Structural to MNI152 linear registration (.mat)"
    # echo "  --warp     : Structural to MNI152 non-linear registration (warp)"
    # echo ""
    # echo "            M0 is to be determined from a saturation recovery"
    # echo "            T1 of tissue (and FA correction) images are also calcualted"
    # echo "   >> Look-Locker flip angle correction - to perform this provide:"

def calib(perf_img, calib_img, calib_method, output_name=None, multiplier=1.0, var=False, log=sys.stdout, **kwargs):
    """
    Do calibration

    :param data: fsl.Image containing data to calibrate
    :param calib_img: fsl.Image containing voxelwise m0 map
    :param calib_method: ``voxelwise`` or ``refregion``
    :param multiplier: Multiplication factor to turn result into desired units
    :param var: If True, assume data represents variance rather than value

    Additional parameters are required for each method.
    """
    if not perf_img:
        raise ValueError("Perfusion data cannot be None")
    if not calib_img:
        raise ValueError("Calibration data cannot be None")

    if calib_method == "voxelwise":
        m0 = get_m0_voxelwise(calib_img, log=log, **kwargs)
    elif calib_method == "refregion":
        m0 = get_m0_refregion(calib_img, log=log, **kwargs)
    else:
        raise ValueError("Unknown calibration method: %s" % calib_method)

    if var:
        log.write("Treating data as variance - squaring M0 correction and multiplier\n")
        m0 = m0**2
        multiplier = multiplier**2

    if isinstance(m0, np.ndarray):
        # If M0 is zero, make calibrated data zero
        calibrated = np.zeros(perf_img.shape)
        calibrated[m0 > 0] = perf_img.data()[m0 > 0] / m0[m0 > 0]
    else:
        calibrated = perf_img.data() / m0

    log.write("Using multiplier for physical units: %f\n" % multiplier)
    calibrated *= multiplier

    if output_name is None:
        output_name = perf_img.iname + "_calib"
    return perf_img.derived(calibrated, name=output_name)

def get_m0_voxelwise(calib_img, gain=1.0, alpha=1.0, tr=None, t1t=None, pct=0.9, brain_mask=None, edgecorr=False, log=sys.stdout):
    """
    Calculate M0 value using voxelwise calibration

    :param calib_img: fsl.Image containing voxelwise m0 map
    :param gain: Calibration gain
    :param alpha: Inversion efficiency
    :param tr: Sequence TR (s) (optional for short TR correction)
    :param t1t: Tissue T1 (optional for short TR correction)
    """
    log.write("Doing voxelwise calibration\n")
    
    # Calculate M0 value
    m0 = calib_img.data() * alpha * gain

    if tr is not None and tr < 5:
        if t1t is not None:
    	    # Correct the M0 image for short TR using the equation from the white paper
            log.write("Correcting the calibration (M0) image for short TR (using T1 of tissue %f)\n" % t1t)
            m0 /= (1 - math.exp(-tr / t1t))
        else:
            log.write("WARNING: tr < 5 (%f) but reference tissue T1 not provided so cannot apply short TR correction\n" % tr)

	# Include partiition co-effcient in M0 image to convert from M0 tissue to M0 arterial
    log.write("Using partition coefficient: %f\n" % pct)
    m0 /= pct
    log.write("Mean M0: %f\n" % np.mean(m0))

    if edgecorr and brain_mask is not None:
        log.write("Doing edge correction\n")
        m0 = edge_correct(m0, brain_mask)
        
    return m0

def edge_correct(m0, brain_mask):
    """
    Correct for (partial volume) edge effects
    """
    brain_mask = brain_mask.data()

    # Median smoothing
    #fslmaths $Mo -fmedian -mas $tempdir/mask -ero $tempdir/calib_ero
    m0 = scipy.ndimage.median_filter(m0, size=3)

    # Erode mask using 3x3x3 structuring element
    mask_ero = scipy.ndimage.morphology.binary_erosion(brain_mask, structure=np.ones([3, 3, 3]), border_value=1)
    m0[mask_ero == 0] = 0

    # Extrapolate remaining data to fit original mask
    # ASL_FILE works slicewise using a mean 5x5 filter on nonzero values, so we will do the same
    for z in range(m0.shape[2]):
        zslice = m0[..., z]
        zslice_extrap = scipy.ndimage.filters.generic_filter(zslice, _masked_mean, footprint=np.ones([5, 5]))
        m0[..., z] = zslice_extrap
    m0[brain_mask == 0] = 0
    
    return m0

def _masked_mean(vals):
    """
    Called by scipy.ndimage.filters.generic_filter

    :param vals: Values in the kernel (a 5x5 square patch centered on the voxel in this case)

    For nonzero voxels, returns the voxel value (i.e. the middle element of vals).
    For zero voxels, returns the mean of non zero voxels in the kernel
    """
    voxel_val = vals[int((len(vals)-1) / 2)]
    if voxel_val == 0:
        nonzero = vals[vals != 0]
        if len(nonzero) > 0:
            return np.mean(nonzero)
        else:
            return 0
    else:
        return voxel_val

def get_tissue_defaults(tiss_type=None):
    """
    Get default T1, T2, T2* and PC for different tissue types
    
    :return: If tiss_type given, tuple of T1, T2, T2* and PC for tissue type.
             Otherwise return dictionary of tissue type name (csf, wm, gm) to tuple
             of T1, T2, T2* and PC for all known types
    """
    tiss_defaults = {
        "csf" : [4.3, 750, 400, 1.15],
        "gm" : [1.3, 100, 50, 0.98],
        "wm" : [1.0, 50, 60, 0.82],
    }
    if tiss_type is None:
        return tiss_defaults
    elif tiss_type.lower() in tiss_defaults:
        return tiss_defaults[tiss_type.lower()]
    else:
        raise ValueError("Invalid tissue type: %s" % tiss_type)

def get_m0_refregion(calib_img, ref_mask=None, brain_mask=None, mode="longtr", gain=1.0, alpha=1.0, log=sys.stdout, **kwargs):
    """
    Do reference region calibration

    FIXME this is not yet complete

    :param calib_img: Calibration image as fsl.Image
    :param ref_mask: Reference region mask image
    :param mode: Calibration mode, ``longtr`` or ``satrecov``
    :param gain: Calibration gain
    :param alpha: Inversion efficiency
    """
    log.write("Doing reference region calibration\n")

    tr = kwargs.get("tr", 3.2)
    te = kwargs.get("te", 0)
    taq = kwargs.get("taq", 0)
    log.write("Using TE=%f\n" % te)

    t2star = kwargs.get("t2star", False)
    if t2star:
        # From Petersen 2006 MRM 55(2):219-232 see discussion
        t2b = 50
    else:
        # From Lu et al. 2012 MRM 67:42-49 have 154ms at 3T during normoxia
        t2b = 150

    # Parameters for reference tissue type: T1, T2, T2*, partition coeffs, FAST seg ID
    # Partition coeffs based on Herscovitch and Raichle 1985 with a blood water density of 0.87
    # GM/WM T2* from Foucher 2011 JMRI 34:785-790
    tissref = kwargs.get("tissref", "csf")
    log.write("Using tissue reference type: %s\n" % tissref)
         
    t1r, t2r, t2rstar, pcr = get_tissue_defaults(tissref)
    if t2star:
        t2r = t2rstar

    # Command line override of default T1, T2, PC
    t1r_img, t2r_img = False, False

    if "t1r" in kwargs:
        t1r = kwargs.get("t1r", None)
        if isinstance(t1r, fsl.Image):
            log.write("Using T1 image for reference region: %s\n" % t1r.iname)
            t1r_img = True
        elif t1r is not None:
            log.write("Using user-specified T1r value: %f\n" % t1r)

    if "t2r" in kwargs:
        t2r = kwargs.get("t2r", None)
        if isinstance(t2r, fsl.Image):
            log.write("Using T2 image for reference region: %s\n" % t2r.iname)
            t2r_img = True
        elif t2r is not None:
            log.write("Using user-specified T2r value: %f\n" % t2r)

    if "t2b" in kwargs:
        t2b = kwargs["t2b"]
        log.write("Using user-specified T2b value: %f\n" % t2b)
        
    if "pcr" in kwargs:
        pcr = kwargs["pcr"]
        log.write("Using user-specified partition coefficient: %f\n" % pcr)

    if t1r is None:
        raise ValueError("T1 for reference tissue has not been set")
    if t2r is None:
        if te != 0:
            raise ValueError("T2 for reference tissue has not been set")
        else:
            t2r = 1.0
    if pcr is None:
        raise ValueError("Partition coefficient for reference tissue has not been set")

    log.write("T1r: %f; T2r: %f; T2b: %f; Part co-eff: %f\n" % (t1r, t2r, t2b, pcr))

    # Check the data and masks
    calib_data = calib_img.data()
    if calib_data.ndim == 4:
        log.write("Taking mean across calibration images\n")
        calib_data = np.mean(calib_data, -1)

    if brain_mask is not None:
        brain_mask = brain_mask.data()
    else:
        brain_mask = np.ones(calib_img.shape[:3])

    if ref_mask is not None:
        log.write("Using supplied reference tissue mask: %s\n" % ref_mask.iname)
        ref_mask = ref_mask.data()
    else:
        # In this case use the brain mask
        log.write("Brain mask is being used as the reference tissue (beware!)\n")
        ref_mask = brain_mask
        
    nonzero = np.count_nonzero(ref_mask)
    if nonzero < 1:
        raise ValueError("Reference mask does not contain any unmasked voxels")
    else:
        log.write("Number of voxels in tissue reference mask: %i\n" % nonzero)

    ### Sensitivity image calculation (if we have a sensitivity image)
    sens_img = kwargs.get("sens_img", None)
    cref_img = kwargs.get("cref_img", None)
    cact_img = kwargs.get("cact_img", None)

    sens_corr = False
    if sens_img:
        log.write("Using sensitivity image: %s\n" % sens_img.iname)
        sens_corr = True
        sens_data = sens_img.data()
    elif cref_img:
        log.write("Calculate sensitivity image from reference image\n")
        sens_corr = True

        # Take the mean (and mask with the mask from the main calib image)
        cref_data = cref_img.data()
        if cref_data.ndim == 4:
            cref_data = np.mean(cref_data, -1)
        cref_data[brain_mask == 0] = 0
        
        if cact_img:
            cact_data = cact_img.data()
        elif mode == "longtr":
            # If the cact image has not been supplied then use the mean of the calib image in longtr mode
            cact_data = calib_data
        else:
            raise ValueError("You must supply an image from the actual coil used for ASL acquisition using --cact (unless you use longtr mode)")
            
        # Take the ratio to give the sensitivity image
        sens_data = cact_data / cref_data
        sens_data[brain_mask == 0] = 0

    log.write("MODE: %s\n" % mode)
    log.write("Calibration gain: %f\n" % gain)

    if mode == "longtr":
        if sens_corr:
            log.write("Applying sensitivity image\n")
            calib_data /= sens_data
        
        # Mask M0 map with tissue reference
        calib_data[ref_mask == 0] = 0
        
        # calcualte T1 of reference region (if a T1 image has been supplied)
        if t1r_img:
            t1r_data = t1r_img.data()
            t1r = np.mean(t1r_data[ref_mask != 0])
            log.write("Calculated T1 of reference tissue: %f\n" % t1r)

        # calcualte T2 of reference region (if a T2 image has been supplied)
        if t2r_img:
            t2r_data = t2r_img.data()
            t2r = np.mean(t2r_data[ref_mask != 0])
            log.write("Calculated T1 of reference tissue: %f\n" % t2r)

        # calculate M0_ref value
        m0 = np.mean(calib_data[ref_mask != 0])
        log.write("mean of reference tissue: %f\n" % m0)
        m0 = m0 / (1 - math.exp(- (tr - taq) / t1r) )
        log.write("M0 of reference tissue: %f\n" % m0)
        
    elif mode == "satrecov":
        # Calibration image is control images and we want to do a saturation recovery fit
        # NB only do the fit in the CSF mask
        options = {
            "method" : "vb",
            "noise" : "white",
            "model" : "satrecov",
            "t1" : t1r,
        }
        
        #deal with TIs
        tis = kwargs.get("tis", [])
        log.write("TIs: %s\n" % str(tis))
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

        # Do fabber within the tissue reference mask with a sensible T1 prior mean

        if sens_corr:
            log.write("Apply sensitivity image to data for reference tisse M0 estimation\n")
            # apply sensitivity map to calibration image - ONLY for the reference tissue calculations
            # fslmaths $calib -div $temp_calib/sens $temp_calib/calib_senscorr
        else:
            pass
            # no sensitivity correction required, but copy image over ready for next command
            # imcp $calib $temp_calib/calib_senscorr

        log.write("Running FABBER within reference tissue mask\n")
        wsp = fsl.Workspace(workdir=kwargs.get("workdir", None))
        wsp.fabber(calib_img, ref_mask, options)
        mean_m0 = fsl.Image("%s/mean_M0t" % wsp.workdir)

        # Calculate M0 value - this is mean M0 of CSF at the TE of the sequence
        m0_value = np.mean(mean_m0.data()[ref_mask != 0])

        log.write("M0 of reference tissue: %f\n" % m0_value)

        if kwargs.get("save_results"):
            # Save useful results
            t1_ref = fsl.Image("%s/mean_T1t" % wsp.workdir)
            m0_ref = fsl.Image("%s/mean_M0t" % wsp.workdir)

            # Do fabber again within whole brain to get estimated T1 of tissue and FA correction (if LL)
            # (note that we do not apply sensitivity correction to the data here - thius is 'built-into' the M0t map)
            log.write("FABBER (again) within whole brain mask\n")

            wsp.fabber(calib_img, brain_mask, options)
            # $fabber --data=$calib --mask=$bmask --output=$temp_calib/satrecovT --data-order=singlefile --model=satrecov --noise=white --method=vb $tislist $llopts $sropts 

            # save useful results to specified output directory
            #imcp $temp_calib/satrecovT/mean_T1t $outdir/T1t
            #imcp $temp_calib/satrecovT/mean_M0t $outdir/M0t
            #if [ ! -z $lfa ]; then
            #imcp $temp_calib/satrecovT/mean_g $outdir/facorr
    else:
        raise ValueError("Unknown reference region mode: %s (Should be satrecov or longtr)" % mode)

    # Use equation to get the M0 value that is needed
    m0 = m0 / math.exp(- te / t2r) #  T2 correction
    m0 = m0 * gain / pcr
    m0 = m0 * math.exp(- te / t2b)
    log.write("M0: %f\n" % m0)

    # Apply calibration to input image
    if sens_corr:
        # Apply sensitivity image
        # fslmaths $infile -div $temp_calib/sens $temp_calib/infile
        pass
    else:
        # imcp $infile $temp_calib/infile
        pass

    if alpha != 1.0:
        log.write("Applying inversion efficiency of: %f\n")
        m0 = m0 * alpha
    
    return m0

def get_csf_mask():
    pass
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
    
if __name__ == "__main__":
    main()
