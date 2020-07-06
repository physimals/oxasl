"""
Modelling support for the QUASAR ASL sequence
---------------------------------------------

Like pCASL, the QUASAR sequence consists of a set of label-control image
pairs acquired at single or multiple delay times, and with the potential
for repeated measurements at each PLD. However each label-control pair is
additionally acquired 7 times, each with a different combination of
arterial flow suppression ('Crusher') gradients. The total number of volumes
in the data set is therefore 2 x 7 x Number of PLDs x Number of repeats.

The standard QUASAR sequence uses 13 PLDs and 6 repeats giving a total of
1092 volumes. TC pairs alternate for the 7 cycles, and one TI is processed
(with all its repeats) at a time (--ibf=tis)

The 'F' form of QUASAR data is averaged over repeats and with ordered so that
the TC pairs for each cycle (at all PLDs) are together (effectively --ibf=cycles)

For processing we require to convert to the 'F' form which would amount to::

  asldata.mean_across_repeats().reorder("ltc")

(c is a ordering character which refers to 'cycles' - this is not yet implemented
in AslImage)

Data in the 'F' form is then differenced to yield a dataset with ordering "tc" and the
final low flip angle cycle is discarded. The resulting dataset will have 13x6=78 
volumes. A mask may be generated from the mean of the control images.

In the calibration step, Fabber is run on the 91 control images using the satrecov
model (including the low-flip angle phase) to obtain the M0 image for calibration.
An edge correction step may be performed using a median filter on the edge voxels.
(alternatively ASL_CALIB can be called instead on the control images using the
satrecov method - note that this is not yet implemented in OXASL)

If PVC required, there is now a registration and segmentation of a structural 
image (or FSL_ANAT output). OXASL should have all this done already. The endpoint
shold be GM and WM PV estimates in native space. Interestingly the T1t output
from the satrecov calibration is used as the native space reference image.

Model based analysis runs the subtracted data on the quasar model with considerable
complexity in prior specification!

Model-free analysis uses asl_mfree.
"""












"""