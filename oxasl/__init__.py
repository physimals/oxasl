"""
Python package for ASL-MRI analysis
===================================

The modules in this package are mostly Python replacements for existing shell script
and C++ code from the oxford_asl FSL module.

For many tools, FSL is required. The ``fslpy`` package is used to wrap required tools
such as BET and FAST.

Design
------

AslImage
~~~~~~~~
The ``oxasl.image.AslImage`` class represents a captured ASL data file. It 
contains information about the acquisition (number/values of TIs/PLDs, repeats,
ordering of label/control image etc). It also has methods which act directly
on the data, for example performing tag-control subtraction, or generation of
a perfusion-weighted image.

    img = AslImage("mpld.nii.gz", plds=[0.25, 0.5, 0.75, 1.0], iaf="tc", order='lrt')
    diffdata = img.diff()
    pwi = img.perf_weighted()

Workspaces
~~~~~~~~~~
The ``workspace`` module contains the ``Workspace`` class which
can be used to build a higher-level interface for complex workflows. A 
workspace is simply an object whose attributes are images, data, etc being
used as part of the workflow. Unlike a normal object, requesting an attribute
that does not exist returns None rather than raising an exception.

    wsp = Workspace()
    print(wsp.asldata) # prints None

A workspace can optionally be associated with a physical directory. If it is,
then setting attributes causes files to be saved in this directory for supported
data types, such as images or 2D matrices. 

    wsp = Workspace(savedir="myoutput")
    wsp.asldata = AslImage("mpld.nii.gz", ntis=1) # Saves myoutput/asldata.nii.gz

A workspace is also associated with 
a log stream (sys.stdout by default) and a prepared logging dictionary ``fsllog``
for passing to FSL Python wrapper commands.

    wsp = Workspace()
    wsp.log.write("Hello World\n")
    wsp.mask = fslmaths(img).bin().run(log=wsp.fsllog)

Module functions
~~~~~~~~~~~~~~~~
Other modules typically contains one or more functions which operate on a 
workspace, in some cases with additional parameters (but not always).

For example the ``calib`` module contains the ``calibrate`` function which
calibrates a perfusion image to physical units using either voxelwise or reference 
region methods. It reads parameters required for this from the workspace, including 
the calibration method to use.

Most of these functions write their output back into the workspace under a standard
name, however in some cases the function might be called on different input images
and might therefore return an image directly, which can be added to the workspace
by the caller under whatever name they prefer 

Command line tools
~~~~~~~~~~~~~~~~~~
Most modules include a ``main()`` function which implements a command line
tool to wrap the module functionality. For example the ``preproc`` module
implements the ``oxasl_preproc`` command line tool which can be used to
do simple preprocessing of ASL data, such as the following to perform
label-control subtraction:

    oxasl_preproc -i asldata.nii.gz --nplds=5 --diff -o asldata_diff.nii.gz

"""

from ._version import __version__, __timestamp__

from .image import AslImage, AslImageOptions
from .preproc import AslPreprocOptions, preprocess
from .workspace import Workspace

__all__ = [
    "__version__",
    "AslImage",
    "AslImageOptions", 
    "AslPreprocOptions",
    "Workspace",
]
