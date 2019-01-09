OXASL command reference
=======================

The main ``OXASL`` command has many options, however most have sensible defaults so in practice 
only a few typically need to be used.

Example usages
--------------

Basic run without calibration or structural data
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In this case we only need specify the structure and acquisition parameters for our ASL data. The data 
we are using is a multi-PLD pCASL acquisition::

    $oxasl -i mpld_asltc --casl --plds=0.25,0.5,0.75,1.0,1.25,1.5 --iaf=tc --ibf=tis --slicedt=0.0452 \
           -o oxasl_out

Here the ``--iaf=tc`` option indicates that the data contains tag-control pairs. ``--ibf=tis`` indicates 
that the TC pairs are stored in blocks corresponding to each TI/PLD, so the first two volumes are 
tag/control for PLD=0.25, while the next two volumes are repeats at the same PLD. If the data 
instead consisted of label-control pairs for *all* the PLDs followed by another block of repeated
images at *all* the PLDs, this would be indicated using ``--ibf=rpt`` (blocks of repeats).

``--slicedt=0.0452`` indicates that this was a 2D readout and each slice in the Z direction adds
0.0452s to the PLD.

``-o oxasl_out`` sets the output directory.

Adding structural information
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

By providing structural information we get the following benefits:

 - Better brain extraction 
 - Output in both native (ASL data) space and also in structural space for overlaying onto structural image
 - Possibility of automatic reference-region calibration
 - Possibility of partial volume correction (see below)

Structural data may be supplied as a T1 weighted image or, better, an output folder from the FSL_ANAT
command run on a structural image. This is preferred because the structural image is already segmented
and bias-corrected so these steps do not need to be done by OXASL. If a structural image is supplied
directly the FSL FAST segmentation tool will be used to do a structural segmentation which can be
a slow process. Here we use FSL_ANAT output::

    $oxasl -i mpld_asltc --casl --plds=0.25,0.5,0.75,1.0,1.25,1.5 --iaf=tc --ibf=tis --slicedt=0.0452 \
           --fslanat T1.anat --senscorr \
           -o oxasl_out --overwrite 
           
``--senscorr`` indicates that the bias-correction field from the FSL_ANAT should be used. We have also
added the ``--overwrite`` option - otherwise OXASL will refuse to run since the output 
directory already exists from our previous run.

Adding calibration
~~~~~~~~~~~~~~~~~~

Calibration enables the output of perfusion maps in physical units, enabling cross-subject and 
cross-session comparisons::

    $oxasl -i mpld_asltc --casl --plds=0.25,0.5,0.75,1.0,1.25,1.5 --iaf=tc --ibf=tis --slicedt=0.0452 \
           --fslanat T1.anat --senscorr \
           -c aslcalib --tr=4.8 --cmethod=single  \
           -o oxasl_out --overwrite 

The calibration image is specified using ``-c aslcalib``. ``--cmethod=single`` indicates that we
want to use a single M0 value for calibration, derived from a reference region. By default OXASL
uses CSF from the ventricles, identified by registering the structural image to a standard 
brain image and using this to mask the ventricles from the CSF segmentation output from either
FAST or (in this case) FSL_ANAT. ``--tr=4.8`` allows a correction to be made for differing T1
value in the tissue and reference. TE can also be similarly provided to correct for differing T2
values.

Distortion correction
~~~~~~~~~~~~~~~~~~~~~

In this case a phase-encoding reversed calibration image (`Blipped` image) was obtained which 
can be used to apply distortion correction using the FSL TOPUP tool::

    $oxasl -i mpld_asltc --casl --plds=0.25,0.5,0.75,1.0,1.25,1.5 --iaf=tc --ibf=tis --slicedt=0.0452 \
           --fslanat T1.anat --senscorr \
           -c aslcalib --tr=4.8 --cmethod=single  \
           --cblip=aslcalib_PA --echospacing=0.06 --pedir=y \
           -o oxasl_out --overwrite 

The echo spacing is given in seconds and the phase encoding direction must also be given ``--pedir=y``. 
Normally this corresponds to scanner co-ordinates, however it is important to view the results of
distortion correction to make sure it is as expected.
 
Full option list
---------------- 

Input ASL image:

-i ASLDATA, --asldata=ASLDATA
                    ASL data file
--iaf=IAF           input ASl format: diff=differenced,tc=tag-control,ct
                    =control-tag,mp=multiphase,ve=vessel-encoded
--order=ORDER       Data order as sequence of 2 or 3 characters:
                    t=TIs/PLDs, r=repeats, l=labelling (tag/control/phases
                    etc). First character is fastest varying
--tis=TIS           TIs (s) as comma-separated list
--plds=PLDS         PLDs (s) as comma-separated list - alternative to
                    --tis
--ntis=NTIS         Number of TIs (for use when processing does not
                    require actual values)
--nplds=NPLDS       Equivalent to --ntis
--rpts=RPTS         Variable repeats as comma-separated list, one per
                    TI/PLD
--nphases=NPHASES   For --iaf=mp, number of phases (assumed to be evenly
                    spaced)
--nenc=NENC         For --iaf=ve, number of encoding cycles
--casl              Acquisition was pseudo cASL (pcASL) rather than pASL
--tau=TAU, --taus=TAU, --bolus=TAU
                    Bolus duration (s). Can be single value or comma
                    separated list, one per TI/PLD
--slicedt=SLICEDT   Timing difference between slices (s) for 2D readout
--sliceband=SLICEBAND
                    Number of slices per pand in multi-band setup
--artsupp           Arterial suppression (vascular crushing) was used
--ibf=IBF           input block format - alternative to --order for
                    compatibility. rpt=Blocks of repeats (i.e. repeats are
                    slowest varying), tis=Blocsk of TIs/PLDs

Structural image:

-s STRUC, --struc=STRUC
                    Structural image
--struc-brain=STRUC_BRAIN, --sbet=STRUC_BRAIN, --struc-bet=STRUC_BRAIN
                    Structural image (brain extracted)
--struc2asl=STRUC2ASL
                    Structural->ASL transformation matrix
--asl2struc=ASL2STRUC
                    ASL->Structural transformation matrix
--wm-seg=WM_SEG     White matter segmentation of structural image
--gm-seg=GM_SEG     Grey matter segmentation of structural image
--csf-seg=CSF_SEG   CSF segmentation of structural image
--fslanat=FSLANAT   FSL_ANAT output directory for structural information
--fastsrc=FASTSRC   Images from a FAST segmentation - if not set FAST will
                    be run on structural image
--struc2std=STRUC2STD
                    Structural to MNI152 linear registration (.mat)
--struc2std-warp=STRUC2STD_WARP
                    Structural to MNI152 non-linear registration (warp)

Main Options:

--wp                Analysis which conforms to the 'white papers' (Alsop
                    et al 2014)
--mc                Motion correct data
--fixbat            Fix bolus arrival time
--fixbolus          Fix bolus duration
--artoff            Do not infer arterial component
--spatial-off       Do not include adaptive spatial smoothing on CBF

Acquisition/Data specific:

--bat=BAT           Estimated bolus arrival time (s) - default=0.7 (pASL),
                    1.3 (cASL)
--batsd=BATSD       Bolus arrival time standard deviation (s)
--t1=T1             Tissue T1 (s)
--t1b=T1B           Blood T1 (s)

Output options:

--save-corrected    Save corrected input data
--save-reg          Save registration information (transforms etc)
--save-basil        Save Basil modelling output
--save-calib        Save calibration output
--save-all          Save all output (enabled when --debug specified)
--output-stddev, --output-std
                    Output standard deviation of estimated variables
--output-var, --vars
                    Output variance of estimated variables
--no-report         Don't try to generate an HTML report

Calibration:

-c CALIB, --calib=CALIB
                    Calibration image
--calib-method=CALIB_METHOD, --cmethod=CALIB_METHOD
                    Calibration method: voxelwise or refregion
--calib-alpha=CALIB_ALPHA, --alpha=CALIB_ALPHA
                    Inversion efficiency
--calib-gain=CALIB_GAIN, --cgain=CALIB_GAIN
                    Relative gain between calibration and ASL data
--tr=TR             TR used in calibration sequence (s)

Voxelwise calibration:

--pct=PCT           Tissue/arterial partition coefficiant
--t1t=T1T           T1 of tissue (s)

Reference region calibration:

--mode=MODE         Calibration mode (longtr or satrevoc)
--tissref=TISSREF   Tissue reference type (csf, wm, gm or none)
--te=TE             Sequence TE (ms)
--t1r=T1R           T1 of reference tissue (s) - defaults: csf 4.3, gm
                    1.3, wm 1.0
--t2r=T2R           T2/T2* of reference tissue (ms) - defaults T2/T2*: csf
                    750/400, gm 100/60,  wm 50/50
--t2b=T2B           T2/T2* of blood (ms) - default T2/T2*: 150/50)
--refmask=REFMASK   Reference tissue mask in perfusion/calibration image
                    space
--t2star            Correct with T2* rather than T2 (alters the default T2
                    values)
--pcr=PCR           Reference tissue partition coefficiant (defaults csf
                    1.15, gm 0.98,  wm 0.82)

longtr mode (calibration image is a control image with a long TR):

satrecov mode (calibration image is a sequnce of control images at various TIs):

--fa=FA             Flip angle (in degrees) for Look-Locker readouts
--lfa=LFA           Lower flip angle (in degrees) for dual FA calibration
--calib-nphases=CALIB_NPHASES
                    Number of phases (repetitions) of higher FA
--fixa              Fix the saturation efficiency to 100% (useful if you
                    have a low number of samples)

Registration:

--regfrom=REGFROM   Registration image (e.g. perfusion weighted image)

Distortion correction using fieldmap:

--fmap=FMAP         fieldmap image (in rad/s)
--fmapmag=FMAPMAG   fieldmap magnitude image - wholehead extracted
--fmapmagbrain=FMAPMAGBRAIN
                    fieldmap magnitude image - brain extracted
--nofmapreg         Do not perform registration of fmap to T1 (use if fmap
                    already in T1-space)

Distortion correction using phase-encode-reversed calibration image (TOPUP):

--cblip=CBLIP       phase-encode-reversed (blipped) calibration image

General distortion correction options:

--echospacing=ECHOSPACING
                    Effective EPI echo spacing (sometimes called dwell
                    time) - in seconds
--pedir=PEDIR       Phase encoding direction, dir = x/y/z/-x/-y/-z
--gdcwarp=GDCWARP   Additional warp image for gradient distortion
                    correction - will be combined with fieldmap or TOPUP
                    distortion correction

Sensitivity correction:

--cref=CREF         Reference image for sensitivity correction
--cact=CACT         Image from coil used for actual ASL acquisition
                    (default: calibration image - only in longtr mode)
--isen=ISEN         User-supplied sensitivity correction in ASL space
--senscorr-auto, --senscorr
                    Apply automatic sensitivity correction using bias
                    field from FAST
--senscorr-off      Do not apply any sensitivity correction

Partial volume correction:
--pvcorr            Apply partial volume correction

Generic:

-o OUTPUT, --output=OUTPUT
                    Output directory
--overwrite         Overwrite output directory if it already exists
-m MASK, --mask=MASK
                    Brain mask image in native ASL space
--optfile=OPTFILE   File containing additional options
--debug             Debug mode
--version           show program's version number and exit
-h, --help          show help message and exit
