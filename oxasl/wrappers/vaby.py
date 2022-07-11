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
from fsl.wrappers import LOAD

from vaby import run

from oxasl.utils import Tee

FABBER_OPTION_MAPPING = {
    "model" : "model_name",
    "continue_from_mvn" : "initial_posterior",
    "save_mvn" : "save_posterior",
}

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

def vaby(options, output=LOAD, ref_nii=None, progress_log=None, log=None):
    """
    Wrapper for VABY model fitting

    :param options: Fabber run options
    :param output: Name of output directory to put results in. The special value
                   LOAD is supported and will cause output to be returned as
                   a dictionary instead.
    :param ref_nii: Optional reference Nibabel image to use when writing output
                    files. Not required if main data is FSL or Nibabel image.
    :param progress_log: File-like stream to log progress percentage to
    :param log: Dictionary of log information as per FSL wrappers
    :return: Dictionary of output data items name:image. The image matches the
             type of the main input data unless this was a file in which case
             an fsl.data.image.Image is returned.
    """
    options = dict(options)
    main_data = options.get("data", None)
    if main_data is None:
        raise ValueError("Main data not specified")

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
        value = options.pop(key)
        newkey = key.replace("-", "_")
        if isinstance(value, Image):
            options[newkey] = value.nibImage
        else:
            options[newkey] = value

    # Streams to capture stdout and stderr and maybe send them elsewhere too
    stdout = Tee()
    stderr = Tee()
    stdout.add(log.get("stdout", None))
    stderr.add(log.get("stderr", None))
    if log.get("tee", False):
        stdout.add(sys.stdout)
        stderr.add(sys.stderr)

    # Deal with differences between Fabber options and Vaby options
    options["method"] = "avb"
    for fabber_key, vaby_key in FABBER_OPTION_MAPPING.items():
        if fabber_key in options:
            options[vaby_key] = options.pop(fabber_key)

    # This is rather silly because fabber can't cope with TIs/rpts as a list but vaby can
    tis, n = [], 1
    while "ti%i" % n in options:
        tis.append(options.pop("ti%i" % n))
        n += 1
    options["tis"] = tis
    rpts, n = [], 1
    while "rpt%i" % n in options:
        rpts.append(options.pop("rpt%i" % n))
        n += 1
    options["rpts"] = rpts

    cmd_output = []
    outdict = {}
    _write_cmd_to_log(log, options)
    runtime, state = run(outdict=outdict, log_stream=sys.stdout, **options)
    print(state["model_mean"])
    print(state.keys())
    print(outdict.keys())

    ret = _Results(cmd_output)
    for k, v in outdict.items():
        k = k.rstrip("_native")
        if k == "posterior":
            ret["finalMVN"] = v
        else:
            ret[k] = v
    print(ret)
    print("Final MVN shape: ", ret["finalMVN"].shape)

    #ret["paramnames"] = fab.get_model_params(options)
    #ret["logfile"] = run.log

    # Write output data or save it as required
    # for data_name, data in run.data.items():
    #     nii = nib.Nifti1Image(data, header=header, affine=affine)
    #     nii.update_header()
    #     nii.header.set_data_dtype(data.dtype)
    #     img = Image(nii)
    #     if output == LOAD:
    #         # Return in-memory data items as the same type as image as the main data
    #         ret[data_name] = _matching_image(main_data, img)
    #     else:
    #         fname = os.path.join(output, data_name)
    #         img.save(fname)

    return ret

def _write_cmd_to_log(log, options):
    if log.get("cmd", None):
        log["cmd"].write("vaby ")
        for key, value in options.items():
            if not isinstance(value, six.string_types) and not isinstance(value, (int, float)):
                value = str(type(value))
            log["cmd"].write("--%s=%s " % (key.replace("_", "-"), value))
        log["cmd"].write("\n")
        
def _matching_image(base_img, img):
    if isinstance(base_img, nib.Nifti1Image):
        return img.nibImage
    elif isinstance(base_img, np.ndarray):
        return img.nibImage.get_data()
    else:
        return img
