OXFORD_ASL analysis
===================

Functions
---------

    Output -> output(basil_name, output_name, subdir)

Outputs data from Basil, renaming it and regridding it to structural/standard space

    OutputMasked -> output_masked(output_name, subdir, mask)

As above, but with masking. Used in PV correction to output parameters by tissue type

    Calibrate -> calibrate(basil_name, output_name, m0val, multiplier, subdir)

Calibrate Basil parameter and output in all spaces

    Report -> report(output_name, subdir, masktype)

Report mean value of parameter within specified named mask - used to report PV values

    Normalise -> normalise(output_name, subdir, masktype)

Output parameter value normalised by dividing by mean within a mask- used in PVC

    Normalise_var -> normalise_var(output_name, subdir, masktype)

Output parameter variance normalised by dividing by mean squared within a mask - used in PVC

    Registration -> registration(regbase, transopt, distout)

Register an image to the structural image

    Calibration -> calibration()

Calculate M0 value for calibration in reference region method

    Dobasil -> do_basil()

Run BASIL

    Dooutput -> do_output(subdir)

Do all the output to a specified output space subdir

    Log
    Warn 

Info or warnings.

Control flow
------------

0. Define usage text and functions as above -> *Line 609/2104*
1. Parse options, display usage if required -> *Line 873/2104*
2. Create directories, decide on output spaces, print welcome, set option defaults including WP mode -> *Line 994/2104*
3. Preprocessing: -> *Line 1275/2104*
  - Handle FSL_ANAT or structural image
  - Remove vol 1 of [calib, cref, cblip)] image, MC and take mean
  - MC asl data using calib as reference
  - Tag control subtraction and reordering using ASL_FILE
  - Get perfusion weighted image (mean over TIs)
4. Registration of ASL data to structural image -> *Line 1320/2104*
5. Segmentation of structural image and pull out sensitivity, PV maps -> *Line 1385/2104*
6. Generate mask -> *Line 1429/2104*
7. Distortion correction and repeat TC subtraction after correction -> *Line 1557/2104
8. Set up BASIL options -> *Line 1795/2104*
9. Run BASIL -> *Line 1827/2104*
10. Re-run registration on perfusion map -> *Line 1840/2104*
11. Do partial volume estimates -> *Line 1881/2104*
12. Calibration -> *Line 1968/2104*
13. Do output
14. Re-run BASIL for PVC -> *Line 2003/2104*
15. Do epoch analysis -> *Line 2048/2104*
16. Final misc outputs and cleanup -> *Line 2104/2104*

Proposed modules
----------------

The main `oxasl` program should end up as a series of calls to separate modules, each
of which has it's own command line tool as well as the direct API

---------------------------
|`image`   | Basic image class for ASL data                         | DONE
|`calib`   | Calibration (voxelwise and ref region)                 | DONE?
|`mask`    | Mask generation                                        | DONE?
|`reg`     | Registration to structural data                        | DONE?
|`preproc` | Preprocessing (differencing, smoothing, moco?)         | DONE?
|`basil`   | Model fitting                                          | DONE?
|`struc`   | Handling structural data                               | DONE?
|`discorr` | Distortion correction                                  | TODO
|`pvc`     | Partial volume corrections                             | TODO
|`epoch`   | Epoch analysis                                         | TODO
|`report`  | Human-readable report document                         | PARTIAL
--------------------------------

Status
------

Steps 1-6 are at least partially implemented and functional
Step 7 (distcorr) is not yet implemented
Steps 8-9 (basil) are implemented
Steps 10-11 (re-registration and PVC) are not implemented
Step 12 (calib) is implemented
Step 13 (output) is incomplete
Step 14 (re-fit for PVC) is not implemented
Step 15 (epoch) is not implemented
Step 16 is not implemented

Data structures
---------------

How to pass data between modules? In particular data generated once should be re-used elsewhere, e.g.

 - Brain-extracted versions of images
 - Segmentations of structural image

Handle this using the `Workspace` class. Modules which need processed versions of data (e.g.
brain extracted structural) call the relevant module function to do this (e.g. `preproc_struc`).
Modules are responsible for:

 1. Not re-running themselves unnecessarily
 2. Not overwriting data provided explicitly (e.g. user-supplied brain extracted structural image)
