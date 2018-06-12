"""
Python API for the the FSL Fabber tool using the C API via the Python ctypes library
"""

import os
import sys
import warnings
import datetime
import time
import glob
from ctypes import CDLL, c_int, c_char_p, c_void_p, c_uint, CFUNCTYPE, create_string_buffer

import six
import numpy as np
import numpy.ctypeslib as npct
import nibabel as nib

if sys.platform.startswith("win"):
    _LIB_FORMAT = "bin\\%s.dll"
    _BIN_FORMAT = "bin\\%s.exe"
elif sys.platform.startswith("darwin"):
    _LIB_FORMAT = "lib/lib%s.dylib"
    _BIN_FORMAT = "bin/%s"
else:
    _LIB_FORMAT = "lib/lib%s.so"
    _BIN_FORMAT = "bin/%s"

def percent_progress(voxel, nvoxels):
    """
    Convenience progress callback which updates a percentage
    """
    complete = 100*voxel/nvoxels
    sys.stdout.write("\b\b\b\b%3i%%" % complete)
    sys.stdout.flush()

def _find_file(current_value, envdir, search_for):
    if current_value is not None:
        return current_value
    elif envdir in os.environ:
        newfpath = os.path.join(os.environ[envdir], search_for)
        if os.path.isfile(newfpath):
            return newfpath
        else:
            return current_value
    else:
        return None

def find_fabber():
    """
    Find the Fabber executable, core library and model libraries, or return None if not found

    :return: A tuple of executable, core library, sequence of model libraries
    """
    ex, lib, models = None, None, []
    for envdir in ("FABBERDIR", "FSLDEVDIR", "FSLDIR"):
        ex = _find_file(ex, envdir, _BIN_FORMAT % "fabber")
        lib = _find_file(lib, envdir, _LIB_FORMAT % "fabbercore_shared")
        models += glob.glob(os.path.join(os.environ.get(envdir, ""), _LIB_FORMAT % "fabber_models_*"))

    return ex, lib, models

def load_options_files(fname):
    """ 
    Load options for a Fabber run from an .fab options file

    :param fname: File name of options file
    """
    options = {}
    with open(fname, "r") as fabfile:
        for line in fabfile.readlines():
            line = line.strip()
            if line and line[0] != "#":
                keyval = line.split("=", 1)
                key = keyval[0].strip()
                if len(keyval) > 1:
                    value = keyval[1].strip()
                else:
                    value = True
                options[key] = value

    return options

def save_options_file(options, fname):
    """
    Save options as a .fab file.
    """
    with open(fname, "w") as fabfile:
        dump_options_file(options, fabfile)
        
def dump_options_file(options, stream):
    """
    Dump to an output stream

    :param stream: Output stream (e.g. stdout or fileobj)
    """
    for key in sorted(options.keys()):
        value = options[key]
        if value == "" or (isinstance(value, bool) and value):
            stream.write("%s" % key)
        elif not isinstance(value, bool):
            stream.write("%s=%s" % (key, value))
        stream.write("\n")

class FabberException(RuntimeError):
    """
    Thrown if there is an error using the Fabber executable or library
    """
    def __init__(self, msg, errcode=None, log=None):
        self.errcode = errcode
        self.log = log
        if errcode is not None:
            RuntimeError.__init__(self, "%i: %s" % (errcode, msg))
        else:
            RuntimeError.__init__(self, msg)

class FabberRun(object):
    """
    The result of a Fabber run
    """
    def __init__(self, data, log):
        self.data = data
        self.log = log
        self.timestamp, self.timestamp_str = self._get_log_timestamp(self.log)

    def write_to_dir(self, dirname, ref_nii=None, extension=".nii.gz"):
        """
        Write the run output to a directory

        This aspires to write the output in a form as close to the command line tool
        as possible, however exact agreement is not guaranteed
        
        :param dirname: Name of directory to write to, will be created if it does not exist
        """
        if not os.path.exists(dirname):
            os.makedirs(dirname)
        
        if not os.path.isdir(dirname):
            raise IOError("Specified directory '%s' exists but is not a directory" % dirname)

        if ref_nii:
            header = ref_nii.header
            affine = ref_nii.header.get_best_affine()
        else:
            header = None
            affine = np.identity(4)

        for data_name, arr in self.data.items():
            nii = nib.Nifti1Image(arr, header=header, affine=affine)
            nii.to_filename(os.path.join(dirname, "%s%s" % (data_name, extension)))
        
        with open(os.path.join(dirname, "logfile"), "w") as logfile:
            logfile.write(self.log)

    def _get_log_timestamp(self, log):
        prefixes = ["start time:", "fabberrundata::start time:"]
        timestamp_str = ""
        for line in log.splitlines():
            line = line.strip()
            for prefix in prefixes:
                if line.lower().startswith(prefix):
                    timestamp_str = line[len(prefix):].strip()
                    try:
                        timestamp = time.strptime(timestamp_str)
                        return timestamp, timestamp_str
                    except ValueError:
                        warnings.warn("Failed to parse timestamp: '%s'" % timestamp_str)
        if log != "":
            warnings.warn("Could not find timestamp in log")
        return datetime.datetime.now(), timestamp_str

class Fabber(object):
    """
    Interface to Fabber in library mode using simplified C API
    """
    def __init__(self, core_lib=None, model_libs=None):
        self.ex, self.core_lib, self.model_libs = find_fabber()
            
        if core_lib:
            self.core_lib = core_lib

        if self.core_lib is None or not os.path.isfile(self.core_lib):
            raise FabberException("Invalid core library - file not found: %s" % self.core_lib)

        if model_libs:
            self.model_libs = set(model_libs)

        for lib in self.model_libs:
            if not os.path.isfile(lib):
                raise FabberException("Invalid models library - file not found: %s" % lib)
           
        self._errbuf = create_string_buffer(255)
        self._outbuf = create_string_buffer(1000000)
        self._progress_cb_type = CFUNCTYPE(None, c_int, c_int)
        self._clib = self._init_clib()
        self._handle = None
        self._init_handle()

    def get_methods(self):
        """ 
        Get known inference methods
        
        :return: Sequence of known inference method names
        """
        self._trycall(self._clib.fabber_get_methods, self._handle, len(self._outbuf), self._outbuf, self._errbuf)
        return self._outbuf.value.splitlines()

    def get_models(self):
        """ 
        Get known models
        
        :return: Sequence of known model names
        """
        self._trycall(self._clib.fabber_get_models, self._handle, len(self._outbuf), self._outbuf, self._errbuf)
        return self._outbuf.value.splitlines()

    def get_options(self, method=None, model=None):
        """
        Get known Fabber options

        :param method: If specified, return options for this method
        :param model: If specified, return options for this model

        Only one of method and model should be specified. If neither are specified, generic
        Fabber options are returned.

        :return: Tuple of options, description. Options is a list of options, each in the form of a dictionary.
        Description is a simple text description of the method or model
        """
        if method and model:
            raise ValueError("get_options: Can't get options for method and model at same time")
        elif method:
            key = "method"
            value = method
        elif model:
            key = "model"
            value = model
        else:
            key = None
            value = None

        self._trycall(self._clib.fabber_get_options, self._handle, key, value, len(self._outbuf), self._outbuf, self._errbuf)
        opt_keys = ["name", "description", "type", "optional", "default"]
        opts = []
        lines = self._outbuf.value.split("\n")
        for opt in lines[1:]:
            if opt:
                opt = dict(zip(opt_keys, opt.split("\t")))
                opt["optional"] = opt["optional"] == "1"
                opts.append(opt)
        return opts, lines[0]

    def get_model_params(self, options):
        """ 
        Get model parameters
        
        :param options: Options dictionary
        :return: Sequence of model parameter names
        """
        self._init_handle()
        for key, value in options.items():
            self._trycall(self._clib.fabber_set_opt, self._handle, str(key), str(value), self._errbuf)

        self._trycall(self._clib.fabber_get_model_params, self._handle, len(self._outbuf), self._outbuf, self._errbuf)

        return self._outbuf.value.splitlines()

    def get_model_outputs(self, options=None):
        """ 
        Get additional model timeseries outputs
        
        :param options: Fabber options
        :return: Sequence of names of additional model timeseries outputs
        """
        self._init_handle()
        if options:
            for key, value in options.items():
                self._trycall(self._clib.fabber_set_opt, self._handle, str(key), str(value), self._errbuf)

        self._trycall(self._clib.fabber_get_model_outputs, self._handle, len(self._outbuf), self._outbuf, self._errbuf)
        return self._outbuf.value.splitlines()

    def model_evaluate(self, options, param_values, nvols, indata=None):
        """
        Evaluate the model with a specified set of parameters

        :param options: Fabber options as key/value dictionary
        :param param_values: Parameter values as a dictionary of param name : param value
        :param nvols: Length of output array - equivalent to number of volumes in input data set
        """
        self._init_handle()
        for key, value in options.items():
            self._trycall(self._clib.fabber_set_opt, self._handle, str(key), str(value), self._errbuf)

        # Get model parameter names and form a sequence of the values provided for them
        self._trycall(self._clib.fabber_get_model_params, self._handle, len(self._outbuf), self._outbuf, self._errbuf)
        model_params = self._outbuf.value.splitlines()
        plist = []
        for param in model_params:
            if param not in param_values:
                raise FabberException("Model parameter %s not specified" % param)
            else:
                plist.append(param_values[param])

        if len(param_values) != len(model_params):
            raise FabberException("Incorrect number of parameters specified: expected %i (%s)" % (len(model_params), ",".join(model_params)))

        ret = np.zeros([nvols,], dtype=np.float32)
        if indata is None: 
            indata = np.zeros([nvols,], dtype=np.float32)

        # Call the evaluate function in the C API
        self._trycall(self._clib.fabber_model_evaluate, self._handle, len(plist), np.array(plist, dtype=np.float32), nvols, indata, ret, self._errbuf)

        return ret

    def _is_data_option(self, key, model_options):
        if key in ("data", "mask", "suppdata", "continue-from-mvn"):
            return True
        elif key.startswith("PSP_byname") and key.endswith("_image"):
            return True
        else:
            model_data_options = [option["name"] for option in model_options if option["type"] in ("IMAGE", "TIMESERIES")]
            return key in model_data_options

    def run(self, options, progress_cb=None):
        """
        Run fabber

        :param options: Fabber options as key/value dictionary. Data may be passed as Numpy arrays, Nifti 
                        images or strings (which are interpreted as filenames)
        :param progress_cb: Callable which will be called periodically during processing

        :return: On success, a FabberRun instance
        """
        if not options.has_key("data"):
            raise ValueError("Main voxel data not provided")

        # Allow input data to be given as Numpy array, Nifti image or filename
        input_data = {}
        model_options, _ = self.get_options(model=options.get("model", "poly"))
        for key in options.keys():
            if self._is_data_option(key, model_options):
                value = options.pop(key)
                if value is None:
                    pass
                elif isinstance(value, nib.Nifti1Image):
                    input_data[key] = nib.get_data()
                elif self._is_data_option(key, model_options) and isinstance(value, six.string_types):
                    input_data[key] = nib.load(value).get_data()
                elif isinstance(value, np.ndarray):
                    input_data[key] = value
                else:
                    raise ValueError("Unsupported type for input data: %s = %s" % (key, type(value)))

        shape = input_data["data"].shape
        nvoxels = shape[0] * shape[1] * shape[2]

        # Make mask suitable for passing to int* c function
        mask = input_data.pop("mask", np.ones(nvoxels))
        mask = np.ascontiguousarray(mask.flatten(order='F'), dtype=np.int32)

        self._init_handle()
        for key, value in options.items():
            # Fabber interprets boolean values as 'option given=True, not given=False. Option value must be blank
            if isinstance(value, bool):
                if value:
                    value = ""
                else:
                    continue
            self._trycall(self._clib.fabber_set_opt, self._handle, str(key), str(value), self._errbuf)
        self._trycall(self._clib.fabber_get_model_params, self._handle, len(self._outbuf), self._outbuf, self._errbuf)
        params = self._outbuf.value.splitlines()

        output_items = []
        if "save-mean" in options:
            output_items += ["mean_" + p for p in params]
        if "save-std" in options:
            output_items += ["std_" + p for p in params]
        if "save-zstat" in options:
            output_items += ["zstat_" + p for p in params]
        if "save-noise-mean" in options:
            output_items.append("noise_means")
        if "save-noise-std" in options:
            output_items.append("noise_stdevs")
        if "save-free-energy" in options:
            output_items.append("freeEnergy")
        if "save-model-fit" in options:
            output_items.append("modelfit")
        if "save-residuals" in options:
            output_items.append("residuals")
        if "save-mvn" in options:
            output_items.append("finalMVN")
        if "save-model-extras" in options:
            output_items += self.get_model_outputs()

        retdata, log = {}, ""
        self._trycall(self._clib.fabber_set_extent, self._handle, shape[0], shape[1], shape[2], mask, self._errbuf)
        for key, item in input_data.items():
            if len(item.shape) == 3:
                size = 1
            else:
                size = item.shape[3]
            item = np.ascontiguousarray(item.flatten(order='F'), dtype=np.float32)
            self._trycall(self._clib.fabber_set_data, self._handle, key, size, item, self._errbuf)

        progress_cb_func = self._progress_cb_type(0)
        if progress_cb is not None:
            progress_cb_func = self._progress_cb_type(progress_cb)

        self._trycall(self._clib.fabber_dorun, self._handle, len(self._outbuf), self._outbuf, self._errbuf, progress_cb_func)
        log = self._outbuf.value
        for key in output_items:
            size = self._trycall(self._clib.fabber_get_data_size, self._handle, key, self._errbuf)

            arr = np.ascontiguousarray(np.empty(nvoxels * size, dtype=np.float32))
            self._trycall(self._clib.fabber_get_data, self._handle, key, arr, self._errbuf)
            if size > 1:
                arr = arr.reshape([shape[0], shape[1], shape[2], size], order='F')
            else:
                arr = arr.reshape([shape[0], shape[1], shape[2]], order='F')
            retdata[key] = arr

        return FabberRun(retdata, log)

    def __del__(self):
        self._destroy_handle()

    def _init_clib(self):
        try:
            clib = CDLL(str(self.core_lib))

            # Signatures of the C functions
            c_int_arr = npct.ndpointer(dtype=np.int32, ndim=1, flags='CONTIGUOUS')
            c_float_arr = npct.ndpointer(dtype=np.float32, ndim=1, flags='CONTIGUOUS')

            clib.fabber_new.argtypes = [c_char_p]
            clib.fabber_new.restype = c_void_p
            clib.fabber_load_models.argtypes = [c_void_p, c_char_p, c_char_p]
            clib.fabber_set_extent.argtypes = [c_void_p, c_uint, c_uint, c_uint, c_int_arr, c_char_p]
            clib.fabber_set_opt.argtypes = [c_void_p, c_char_p, c_char_p, c_char_p]
            clib.fabber_set_data.argtypes = [c_void_p, c_char_p, c_uint, c_float_arr, c_char_p]
            clib.fabber_get_data_size.argtypes = [c_void_p, c_char_p, c_char_p]
            clib.fabber_get_data.argtypes = [c_void_p, c_char_p, c_float_arr, c_char_p]
            clib.fabber_dorun.argtypes = [c_void_p, c_uint, c_char_p, c_char_p, self._progress_cb_type]
            clib.fabber_destroy.argtypes = [c_void_p]

            clib.fabber_get_options.argtypes = [c_void_p, c_char_p, c_char_p, c_uint, c_char_p, c_char_p]
            clib.fabber_get_models.argtypes = [c_void_p, c_uint, c_char_p, c_char_p]
            clib.fabber_get_methods.argtypes = [c_void_p, c_uint, c_char_p, c_char_p]

            clib.fabber_get_model_params.argtypes = [c_void_p, c_uint, c_char_p, c_char_p]
            clib.fabber_get_model_outputs.argtypes = [c_void_p, c_uint, c_char_p, c_char_p]
            clib.fabber_model_evaluate.argtypes = [c_void_p, c_uint, c_float_arr, c_uint, c_float_arr, c_float_arr, c_char_p]
            return clib
        except Exception as exc:
            raise RuntimeError("Error initializing Fabber library: %s" % str(exc))

    def _init_handle(self):
        # This is required because currently there is no CAPI function to clear the options.
        # So we destroy the old Fabber handle and create a new one
        self._destroy_handle()    
        self._handle = self._clib.fabber_new(self._errbuf)
        if self._handle is None:
            raise RuntimeError("Error creating fabber context (%s)" % self._errbuf.value)

        for lib in self.model_libs:
            self._trycall(self._clib.fabber_load_models, self._handle, lib, self._errbuf)

    def _destroy_handle(self):
        if hasattr(self, "_handle"):
            if self._handle:
                self._clib.fabber_destroy(self._handle)
                self._handle = None

    def _trycall(self, call, *args):
        ret = call(*args)
        if ret < 0:
            raise FabberException(self._errbuf.value, ret, self._outbuf.value)
        else:
            return ret
