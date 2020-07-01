Example using multiphase ASL data
=================================

This example shows how to process multiphase data within OXASL. Multiphase data
is collected at a range of inversion phases, so rather than pairs of tag/control
images, one has a set of images covering angles between 0 and 180.

The multiphase plugin performs a multiphase fitting step after preprocessing but
before running the kinetic model inversion. This step essentially reduces the 
multiphse data to a form equivalent to subtracted tag-control data.

Running OXASL on multiphase ASL data
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The key distinction between running multiphase data and regular ASL data is the ``--iaf``
option which should be set to ``--iaf=mp``.

In addition **one** of the following must be specified

 - ``--nphases=<Number of phases>``
 - ``--phases=<Comma separated list of phases in degrees>``

If ``--nphases`` is used, the phases are spread equally between 0 and 180 degrees.

A minimal example command line for multiphase ASL data would be::

    oxasl -i data.nii.gz --casl --tau=1.4 --plds=0.4 \
          --iaf=mp --nphases=8 \
          --output=oxasl_mp

Other standard OXASL options can also be used to enable calibration, provide structural data,
or apply preprocessing corrections to the data.

Additional options for multiphase data
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``--mp-spatial``
----------------

This option uses spatial regularization during the multiphase fitting step. It is
independent of the ``--spatial`` option which controls spatial regularization during
the kinetic model fitting step.

``--mp-spatial-phase``
----------------------

When using ``--mp-spatial`` this option places the spatial prior on the phase offset
parameter which can give better spatial regularization than the default which performs
spatial regularization on the magnitude.

``--mp-biascorr``
~~~~~~~~~~~~~~~~~

This option applies a bias correction to the multiphase fitting step. In the presence
of significant noise, the multiphase fitting is biased towards higher magnitudes with
larger amounts of noise. This option performs a multi-step process which aims to 
reduce the bias by initially fitting the phase offset and then fixing this while fitting
the magnitude. The full method is described in:

 - Msayib, Y., et al. *Robust Estimation of Quantitative Perfusion from Multi-Phase 
   Pseudo-Continuous Arterial Spin Labelling.* Magnetic Resonance in Medicine, 2019.
 
``--mp-biascorr-sv``
~~~~~~~~~~~~~~~~~~~~

Number of supervoxels in the bias-correction step

``--mp-biascorr-comp``
~~~~~~~~~~~~~~~~~~~~~~

Supervoxel compactness in the bias correction step

``--mp-biascorr-sigma``
~~~~~~~~~~~~~~~~~~~~~~~

Pre-supervoxel smoothing parameter in the bias correction step

``--mp-options=<options file>``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This allows additional options to be passed to the multiphase fitting step
(for example ``max-iterations=20`` to increase the number of iterations).
The options should be placed in a text file.

The log output
~~~~~~~~~~~~~~

The command line output is similar to the label-control case, however there is an additional
multiphase fitting step between the preprocessing and the kinetic model inversion, which will
look something like this (note depending on the multiphase options used there may be fewer
steps performed)::

    Performing multiphase decoding:
    - Using supervoxel-based bias correction
    - Number of supervoxels: 20  - Compactness: 0.050000  - Pre-smoothing width: 0  - Step 1: Running initial biased fit 100% - DONE
    - Step 2: Fitting mean signal in supervoxels 100% - DONE
    - Step 3: Running final fit with fixed phase 100% - DONE

    DONE multiphase decoding

Output images
~~~~~~~~~~~~~

The output images are as usual found in the ``oxasl_out/output`` directory. Output files from the multiphase
fitting step are found in ``oxasl_out/mp``.

The usual OXASL output images are produced, for example:

 - ``perfusion.nii.gz`` - This is the relative perfusion image
 - ``arrival.nii.gz`` - This is the inferred bolus arrival time image
 - ``aCBV.nii.gz`` - This is the inferred macrovascular signal image containing arterial volume fraction as a percentage
 - ``mask.nii.gz`` - This is the binary brain mask used in the analysis

Calibrated outputs are also produced if calibration data is supplied, and structural space outputs are
also produced where structural data is available.

Summary report
~~~~~~~~~~~~~~

Currently the multiphase fitting step does not generate any information in the summary report. This will be 
improved in the future!
