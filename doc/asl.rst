Bayesian Inference for Arterial Spin Labelling MRI
==================================================

.. image:: images/basil_perfusion.jpg
   :scale: 100 %
   :alt: BASIL perfusion image
   :align: right

Arterial Spin Labeling (ASL) MRI is a non-invasive method for the quantification 
of perfusion. Analysis of ASL data typically requires the inversion of a kinetic 
model of labeled blood-water inflow along with a separate calculation of the equilibrium 
magnetization of arterial blood. 

The OXASL toolbox derives from the BASIL package which performs the kinetic modelling
using Bayesian inference principles. The toolbox was orginally developed for 
multi delay (inversion time) data where it can be used to greatest effect, but 
is also sufficiently fleixble to deal with the widely used single delay form 
of acquisition.

.. note::
   If you want to perform analysis of a functional experiment with ASL data, i.e. one where 
   you want to use a GLM, then you should consult the perfusion section of 
   `FEAT <https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FEAT/UserGuide>`_, 
   or if you have dual-echo (combined BOLD and ASL) data then consult 
   `FABBER <https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FABBER>`_.

For single delay ASL data kinetic model inversion is relatively trivial and 
solutions to the standard model have been described in the literature. However,
there are various advantages to aquiring ASL data at multiple times 
post-inversion and fitting the resultant data to a kinetic model. This 
permits problems in perfusion estimation associated with variable bolus arrival 
time to be avoided, since this becomes a parameter of the model whose value is 
determined from the data. Commonly the model fitting will be performed with a 
least squares technique providing parameter estimates, e.g. perfusion and bolus 
arrival time. In contrast to this BASIL uses a (fast) Bayesian inference method 
for the model inversion, this provides a number of advantages:

 - Voxel-wise estimation of perfusion and bolus arrival time along with parameter 
   variance (allowing confidence intervals to be calculated).

 - Incorporation of natural varaibility of other model parameters, e.g. values of T1,
   T1b and labeling/bolus duration.

 - Spatial regularization of the estimated perfusion image.

 - Correction for partial volume effects (where the appropriate segmentation 
   information is available).

While the first two apply specfically to the case of mulitple delay data, the latter 
are also applicable to single delay ASL and are only available using the Bayesian 
technique employed by BASIL.