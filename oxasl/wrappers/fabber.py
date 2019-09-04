"""
FSL-compatible wrappers for Fabber tools
"""
from __future__ import absolute_import

import sys
import os

import six
import numpy as np
import nibabel as nib

from fsl.data.image import Image
from fsl.wrappers import LOAD, wrapperutils  as wutils
import fsl.utils.assertions as asrt
from fabber import Fabber, FabberException, percent_progress

from oxasl.utils import Tee

def _matching_image(base_img, img):
    if isinstance(base_img, nib.Nifti1Image):
        return img.nibImage
    elif isinstance(base_img, np.ndarray):
        return img.nibImage.get_data()
    else:
        return img

class _Results(dict):
    """
    Nicked from fsl.wrapperutils
    """
    def __init__(self, output):
        dict.__init__(self)
        self.__output = output

    @property
    def output(self):
        """Access the return value of the decorated function. """
        return self.__output

def fabber(options, output=LOAD, ref_nii=None, progress_log=None, **kwargs):
    """
    Wrapper for Fabber tool

    This is not a 'conventional' FSL command line tool wrapper. Rather it is
    using the Fabber Python API which interfaces to the C++ code using either
    the pure-C api and Python ctypes or its own CLI wrapper.

    The main reason for not using a conventional wrapper is that Fabber
    can take an arbitrary number of inputs including image inputs so the
    @fileOrImage type decorators don't really work well. Also you can't
    tell what image inputs you might have until you query the individual
    model for its options. All of this complexity is therefore hidden in
    the generic Fabber python API.

    Nevertheless we aim to replicate the interface of FSL tools as much as
    possible, e.g. accepting fsl.data.image.Image instances for input
    data, and respecting the LOAD special parameter to indicate which
    data items should be returned as Image instances.

    :param options: Fabber run options
    :param output: Name of output directory to put results in. The special value
                   LOAD is supported and will cause output to be returned as
                   a dictionary instead.
    :param ref_nii: Optional reference Nibabel image to use when writing output
                    files. Not required if main data is FSL or Nibabel image.
    :param progress_log: File-like stream to logging progress percentage to
    :return: Dictionary of output data items name:image. The image matches the
             type of the main input data unless this was a file in which case
             an fsl.data.image.Image is returned.
    """
    extra_search_dirs = kwargs.pop("fabber_dirs", ())
    fab = Fabber(*extra_search_dirs)

    options = dict(options)
    main_data = options.get("data", None)
    if main_data is None:
        raise ValueError("Main data not specified")

    if output != LOAD and not os.path.exists(output):
        os.makedirs(output)

    # Get a reference Nibabel image to use when generating output
    if not ref_nii:
        if isinstance(main_data, Image):
            ref_nii = main_data.nibImage
        else:
            ref_nii = Image(main_data).nibImage

    if ref_nii:
        header = ref_nii.header
        affine = ref_nii.header.get_best_affine()
    else:
        header = None
        affine = np.identity(4)

    # Replace fsl.Image objects with the underlying nibabel object. The Fabber
    # Python API can already handle Numpy arrays, nibabel images and filenames
    for key in list(options.keys()):
        value = options[key]
        if isinstance(value, Image):
            options[key] = value.nibImage

    # Streams to capture stdout and stderr and maybe send them elsewhere too
    stdout = Tee()
    stderr = Tee()

    # Deal with standard keyword arguments
    ret_exitcode = kwargs.pop("exitcode", False)
    ret_stdout = kwargs.pop("stdout", True)
    ret_stderr = kwargs.pop("stderr", False)
    log = kwargs.pop("log", {})

    stdout.add(log.get("stdout", None))
    stderr.add(log.get("stderr", None))
    if log.get("tee", False):
        stdout.add(sys.stdout)
        stderr.add(sys.stderr)

    if kwargs.pop("submit", False):
        raise ValueError("submit not supported for Fabber")

    exception = None
    cmd_output = []
    ret = _Results(cmd_output)
    try:
        ret["paramnames"] = fab.get_model_params(options)
        if log.get("cmd", None):
            log["cmd"].write("Using fabber:\n  core lib=%s\n  core_exe=%s\n  model libs=%s\n  model exes=%s\n" % (fab.core_lib, fab.core_exe, fab.model_libs, fab.model_exes))
            log["cmd"].write("fabber ")
            for key, value in options.items():
                if not isinstance(value, six.string_types) and not isinstance(value, (int, float)):
                    value = str(type(value))
                log["cmd"].write("--%s=%s " % (key.replace("_", "-"), value))
            log["cmd"].write("\n")
        progress_cb = None
        if progress_log:
            progress_cb = percent_progress(progress_log)
        run = fab.run(options, progress_cb)
        ret["logfile"] = run.log

        # Write output data or save it as required
        for data_name, data in run.data.items():
            img = Image(nib.Nifti1Image(data, header=header, affine=affine))
            if output == LOAD:
                # Return in-memory data items as the same type as image as the main data
                ret[data_name] = _matching_image(main_data, img)
            else:
                fname = os.path.join(output, data_name)
                img.save(fname)

    except FabberException as exc:
        # Error while actually running Fabber - may raise later
        # or replace with exit code
        exception = exc
        stderr.write(str(exc) + "\n")

    if ret_stdout:
        cmd_output.append(str(stdout))
    if ret_stderr:
        cmd_output.append(str(stderr))
    if ret_exitcode:
        cmd_output.append(int(exception is not None))

    if exception is not None and not ret_exitcode:
        raise exception

    return ret

@wutils.fileOrImage('mvn', 'output', 'valim', 'varim', 'mask')
@wutils.fslwrapper
def mvntool(mvn, param, **kwargs):
    """
    Wrapper for MVNTOOL.

    :param mvn: Input MVN data
    :param param: Parameter index (or name, but only if param-list is specified)

    For other arguments, see command line tool for now
    """
    asrt.assertIsNifti(mvn)

    valmap = {
        'write' : wutils.SHOW_IF_TRUE,
        'new' : wutils.SHOW_IF_TRUE,
    }

    cmd = ['mvntool', '--input={}'.format(mvn), '--param={}'.format(param)]
    cmd += wutils.applyArgStyle('--=', valmap=valmap, **kwargs)

    return cmd
