"""
Workspace for executing FSL commands
"""

import os
import sys
import errno
import six

import numpy as np
import yaml

from fsl.data.image import Image

from oxasl import AslImage
from oxasl.reporting import Report

class Workspace(object):
    """
    A workspace for data processing

    The contents of a workspace are modelled as attributes on the workspace
    object. An initial set of contents can be specified using keyword arguments.
    For command line tools, typically these are provided directly from the
    OptionParser results, e.g.

        options, args = parser.parse_args(sys.argv)
        wsp = Workspace(**vars(options))

    A workspace may optionally be associated with a physical directory. In this
    case certain types of objects will automatically be saved to the workspace.
    Supported types are currently:

         - ``fsl.data.image.Image`` - Saved as Nifti
         - 2D Numpy array - Saved as ASCII matrix

    To avoid saving a particular item, use the ``add`` method rather than
    directly setting an attribute, as it supports a ``save`` option.
    """

    def __init__(self, savedir=None, input_wsp="input", parent=None, defaults=("corrected", "input"), auto_asldata=False, **kwargs):
        """
        Create workspace

        :param savedir: If specified, use this path to save data. Will be created
                        if it does not not already exist
        :param separate_input: If True a sub workspace named 'input' will be created
                               and initial data will be stored there. This sub workspace
                               will also be checked by default when attributes are 
                               requested from the main workspace.
        :param auto_asldata: If True, automatically create an AslImage attribute 
                             from the input keyword arguments
        :param log:     File stream to write log output to (default: sys.stdout)
        """
        self.savedir = None # Have to set this first otherwise setattr fails!
        self._parent = parent
        self._defaults = list(defaults)
        self._stuff = {}
        if savedir is not None:
            self.savedir = os.path.abspath(savedir)
            mkdir(savedir, log=self.ifnone("log", kwargs.get("log", sys.stdout)))

        # Defaults - these can be overridden by kwargs but might be
        # already defined in parent workspace
        if self.log is None:
            self.log = kwargs.pop("log", sys.stdout)
        if self.debug is None:
            self.debug = kwargs.pop("debug", False)
        if self.report is None:
            self.report = kwargs.pop("report", Report())
        
        # Default log configuration for FSL wrapper commands
        if self.fsllog is None:
            self.fsllog = kwargs.pop("fsllog", None)
            if not self.fsllog:
                self.fsllog = {"stderr" : self.log}
                if self.debug:
                    self.fsllog.update({"stdout" : self.log, "cmd" : self.log})

        # Set kwargs as attributes in input workspace (if configured)
        if input_wsp:
            input_wsp = self.sub(input_wsp)
        else:
            input_wsp = self

        for key, value in kwargs.items():
            setattr(input_wsp, key, value)

        # Auto-generate ASLImage object
        if auto_asldata:
            if kwargs.get("asldata", None) is None:
                raise ValueError("Input ASL file not specified\n")
            input_wsp.asldata = AslImage(self.asldata, **kwargs)

    def __getattr__(self, name):
        if name in self._defaults:
            return None
        for wsp in self._defaults:
            default_wsp = getattr(self, wsp)
            if isinstance(default_wsp, Workspace):
                val = getattr(default_wsp, name)
                if val is not None:
                    return val
        if self._parent is not None:
            return getattr(self._parent, name)
        return None

    def __setattr__(self, name, value):
        self.set_item(name, value)

    def ifnone(self, attr, alternative):
        """
        Return the value of an attribute, if set and not None, or 
        otherwise the supplied alternative
        """
        ret = getattr(self, attr, None)
        if ret is None:
            ret = alternative
        return ret

    def set_item(self, name, value, save=True, save_name=None, save_fn=None):
        """
        Add an item to the workspace

        Normally this is achieved by assigning to an attribute
        directly, however this function exists to allow greater control
        where required.

        The item will be set as an attribute on the workspace
        If a save directory is configured and the value is a supported
        type it will be saved there. This can be disabled by setting save=False

        :param name: Name, must be a valid Python identifier
        :param value: Value to set
        :param save: If False, do not save item even if savedir defined
        :param save_name: If specified, alternative name to use for saving this item
        :param save_fn: If specified, Callable which generates string representation of
                        value suitable for saving the item to a file
        """
        if value is not None and save and self.savedir:
            if not save_name:
                save_name = name

            if save_fn is not None:
                with open(os.path.join(self.savedir, save_name), "w") as tfile:
                    tfile.write(save_fn(value))
            elif isinstance(value, Image):
                # Save as Nifti file
                fname = os.path.join(self.savedir, save_name)
                value.save(fname)
                value.name = save_name
                if type(value) == Image:
                    # Create fresh image to avoid data caching
                    # NB do not want subclasses
                    value = Image(fname, name=save_name)
                else:
                    # Release cached data otherwise
                    value.nibImage.uncache()
            elif isinstance(value, np.ndarray) and value.ndim == 2:
                # Save as ASCII matrix
                with open(os.path.join(self.savedir, save_name), "w") as tfile:
                    tfile.write(matrix_to_text(value))
            elif not name.startswith("_") and isinstance(value, (int, float, six.string_types)):
                # Save other attributes in JSON file
                self._stuff[name] = value
                self._save_stuff()

        super(Workspace, self).__setattr__(name, value)

    def uncache(self):
        """
        Release in-memory data for all images in this workspace
        """
        for attr in dir(self):
            val = getattr(self, attr)
            if isinstance(val, Image):
                val.nibImage.uncache()
            elif isinstance(val, Workspace) and val != self and not attr.startswith("_"):
                val.uncache()

    def sub(self, name, parent_default=True, **kwargs):
        """
        Create a sub-workspace, (i.e. a subdir of this workspace)

        This inherits the log configuration from the parent workspace. If
        a savedir is defined, it will be created with a savedir which is
        a subdirectory of the original workspace. Additional data may be 
        set using keyword arguments. The child-workspace will be available 
        as an attribute on the parent workspace.
        
        :param name: Name of sub workspace
        :param parent_default: If True, attribute values default to the parent workspace
                               if not set on the sub-workspace 
        """
        if self.savedir:
            savedir = os.path.join(self.savedir, name)
        else:
            savedir = None
        
        if parent_default and name not in self._defaults:
            parent = self
        else:
            parent = None
            kwargs["log"] = self.log
            kwargs["debug"] = self.debug

        sub_wsp = Workspace(savedir=savedir, parent=parent, input_wsp=None, **kwargs)
        setattr(self, name, sub_wsp)
        return sub_wsp

    def _save_stuff(self):
        with open(os.path.join(self.savedir, "_oxasl.yml"), "w") as tfile:
            yaml.dump(self._stuff, tfile, default_flow_style=False)

def matrix_to_text(mat):
    """
    Convert matrix array to text using spaces/newlines as col/row delimiters
    """
    rows = []
    for row in mat:
        rows.append(" ".join([str(v) for v in row]))
    return "\n".join(rows)

def text_to_matrix(text):
    """
    Convert space or comma separated file to matrix
    """
    fvals = []
    ncols = -1
    lines = text.splitlines()
    for line in lines:
        # Discard comments
        line = line.split("#", 1)[0].strip()
        # Split by commas or spaces
        vals = line.replace(",", " ").split()
        # Ignore empty lines
        if not vals: continue
        # Check correct number of columns
        if ncols < 0: ncols = len(vals)
        elif len(vals) != ncols:
            raise ValueError("File must contain a matrix of numbers with fixed size (rows/columns)")
        # Check all data is numeric
        for val in vals:
            try:
                float(val)
            except:
                raise ValueError("Non-numeric value '%s' found in matrix text" % val)
        fvals.append([float(v) for v in vals])     
    return np.array(fvals)

def mkdir(dirname, fail_if_exists=False, warn_if_exists=True, log=sys.stdout):
    """
    Create a directory, including necessary subdirs
    """
    try:
        os.makedirs(dirname)
    except OSError as exc:
        if exc.errno == errno.EEXIST:
            if fail_if_exists: raise
            elif warn_if_exists: log.write("WARNING: mkdir - Directory %s already exists\n" % dirname)
    return os.path.abspath(dirname)
