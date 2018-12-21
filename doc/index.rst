OXASL - ASL-MRI analysis pipeline
=================================

OXASL is a package for performing Bayesian analysis of Arterial Spin Labelling
MRI data. Features of the toolbox include:

 - Support for single or multi delay (inversion time) data
 - pCASL or PASL acquisitions
 - Support for acquisitions with variable repeats and/or labelling durations
 - Handles label-control or subtracted input data in various ordering conventions
 - Calibration using the reference region or voxelwise methods
 - Structural registration and output in structural space
 - GM/WM Partial volume correction
 - HTML summary report

Plugins are also available for:

 - Vessel-encoded ASL data

OXASL works within an FSL environment which must be installed.

.. toctree::
   :maxdepth: 2
   :caption: Contents:
  
   asl
   oxasl
   api
