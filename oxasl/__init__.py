"""
Python package for ASL-MRI analysis
===================================

The modules in this package are mostly Python versions of existing shell script
and C++ code from the oxford_asl FSL module.

For many tools, FSL is required.

Design
------

AslImage
~~~~~~~~
The ``oxasl.image.AslImage`` class represents a captured ASL data file. It 
contains information about the acquisition (number/values of TIs/PLDs, repeats,
ordering of label/control image etc). It also has methods which act directly
on the data, for example performing tag-control subtraction, or generation of
a perfusion-weighted image.

    img = AslImage("mpld.nii.gz", plds=[0.25, 0.5, 0.75, 1.0], order='prt')
    diffdata = img.diff()
    pwi = img.perf_weighted()

Module functions
~~~~~~~~~~~~~~~~
Other modules typically contains one or more functions which take data images
and process them in some way. For example the ``calib`` module calibrates a
perfusion image using either voxelwise or reference region methods.

These functions take parameters they need and return data and/or images as
appropriate.

    perf_calib_img = calib(perf_img, calib_img, method="voxelwise")

Command line tools
~~~~~~~~~~~~~~~~~~
Most modules include a ``main()`` function which implements a command line
tool to wrap the module functionality. For example the ``preproc`` module
implements the ``oxasl_preproc`` command line tool which can be used to
do simple preprocessing of ASL data, such as the following to perform
label-control subtraction:

    oxasl_preproc -i asldata.nii.gz --nplds=5 --diff -o asldata_diff.nii.gz

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

The base ``Workspace`` class simply does the above and can be subclassed to
include additional methods which act on data already in the workspace to produce
new data. The ``AslWorkspace`` class has methods used by the ``oxasl`` tool,
for example the ``generate_mask`` method which examines the data present in the
workspace (whether there is a structural image, etc) and generates a brain
mask in the space of the ASL data and sets the ``mask`` attribute to this
image.
"""

from ._version import __version__, __timestamp__

from .image import AslImage, AslImageOptions
from .calib import calib
from .preproc import AslPreprocOptions, preprocess
from . import basil

__all__ = [
    "__version__",
    "AslImage",
    "AslImageOptions", 
    "AslPreprocOptions",
    "preprocess",
    "basil", 
    "calib", 
]
