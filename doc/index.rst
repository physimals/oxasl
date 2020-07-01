
.. image:: images/oxasl.png
   :scale: 30 %
   :alt: OXASL logo
   :align: right

OXASL - ASL-MRI analysis pipeline
=================================

OXASL is a package for performing Bayesian analysis of Arterial Spin Labelling
MRI data. Features of the toolbox include:

 - Support for single or multi delay (inversion time) pCASL or PASL data including acquisitions with variable repeats and/or labelling durations
 - Handles label-control or subtracted input data in various ordering conventions
 - Calibration using the reference region or voxelwise methods
 - Structural registration and output in structural space
 - GM/WM Partial volume correction
 - HTML summary report

Plugins are also available for:

 - Vessel-encoded ASL data
 - Multiphase ASL data

OXASL works within an FSL environment which must be installed.

.. toctree::
   :maxdepth: 1
   :caption: Contents:
  
   asl
   download
   walkthrough_cl
   walkthrough_gui
   example_ve
   example_mp
   oxasl
   api
