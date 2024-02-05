"""
OXASL - Workspace to store images and data associated with a processing pipeline

This class is conceptually simple - you can store pretty much any data by setting
an attribute on the workspace and retrieve the resulting data as the same attribute.
The workspace is backed by a directory and image data is save to files rather than
store in memory.

This hides considerable complexity. Here are a few of the issues:

 - To ensure that image data really is kicked out of memory we create an ImageProxy
   object for each image. This simply stores the filename and metadata associated
   with an image. ``__getattribute__`` is overridden to return ImageProxy attributes
   as the underlying Image.

 - There is a special ImageProxy for an AslImage. This might go away if we can
   represent the full state of an AslImage using metadata alone.

Copyright (c) 2008-2020 Univerisity of Oxford
"""

import os
import sys
import errno
import glob
import shutil
import tempfile

import six
import numpy as np
import pandas as pd
import yaml

from fsl.data.image import Image

from oxasl import AslImage
from oxasl.reporting import Report
from oxasl.utils import Tee

def current_umask():
    """
    Bonkers workaround to get the current umask which is impossible without setting it too
    """
    tmp = os.umask(0o022)
    os.umask(tmp)
    return tmp

class ImageProxy(object):
    """
    Reference to a saved Image and it's metadata
    """
    def __init__(self, fname, md=None):
        self._fname = fname
        self._md = md

    def img(self):
        """
        Return an Image object by loading from the file
        """
        img = Image(self._fname, loadData=False)
        if self._md:
            for key, value in self._md:
                img.setMeta(key, value)
        return img

class AslImageProxy(ImageProxy):
    """
    Reference to a saved AslImage and it's metadata
    """

    def img(self):
        """
        Return an AslImage object from the file and stored metadata
        """
        return AslImage(self._fname, loadData=False, **self._md)

class Workspace(object):
    """
    A workspace for data processing

    The contents of a workspace are modelled as attributes on the workspace
    object. An initial set of contents can be specified using keyword arguments.
    For command line tools, typically these are provided directly from the
    OptionParser results, e.g.

        options, args = parser.parse_args(sys.argv)
        wsp = Workspace(**vars(options))

    A workspace is always associated with a physical directory. Certain types of
    objects are automatically be saved to the workspace. If no save directory
    is specified a temporary directory is created

    Supported types are currently:

         - ``fsl.data.image.Image`` - Saved as Nifti
         - 2D Numpy array - Saved as ASCII matrix

    All other attributes are serialized to YAML and stored in a special
    ``_oxasl.yml`` file.

    To avoid saving a particular item, use the ``add`` method rather than
    directly setting an attribute, as it supports a ``save`` option.
    """

    def __init__(self, savedir=None, input_wsp="input", parent=None, search_childs=("prequantify", "filter", "corrected", "preproc", "input"), auto_asldata=False, **kwargs):
        """
        Create workspace

        :param savedir: If specified, use this path to save data. Will be created
                        if it does not not already exist. If not specified a temporary
                        directory will be created so data can be flushed to disk rather
                        than held in memory.
        :param input_wsp: Name of child workspace to be created and initial data stored within
        :param parent: Parent Workspace which will be searched if attribute not found within
                       this workspace
        :param search_childs: Sequence of names of child workspaces to be searched for attributes
                              (in order - first to appear is first to be checked)
        :param auto_asldata: If True, automatically create an AslImage attribute
                             from the input keyword arguments
        :param log: Optional file stream to write log output to (default: parent workspace log or sys.stdout)
        """
        # Have to set this first otherwise setattr fails!
        if savedir is not None:
            savedir = os.path.abspath(savedir)
            create_savedir = True
        else:
            savedir = tempfile.mkdtemp(prefix="oxasl_wsp")
            create_savedir = False
        self.set_item("savedir", savedir, save=False)

        self._parent = parent
        self._search_childs = list(search_childs)
        self._stuff = {}
        if create_savedir:
            if parent is not None:
                warn_overwrite = not parent._overwrite_warned
            else:
                warn_overwrite = True
            mkdir(savedir, log=self.ifnone("log", kwargs.get("log", sys.stdout)), 
                  warn_if_exists=warn_overwrite)
            self._overwrite_warned = True
    
        # Defaults - these can be overridden by kwargs but might be
        # already defined in parent workspace
        if "log" in kwargs or self.log is None:
            self.log = kwargs.pop("log", None)
            if self.log is None:
                logfile = open(os.path.join(self.savedir, "logfile"), "w")
                self.log = Tee(sys.stdout, logfile)
        if "debug" in kwargs or self.debug is None:
            self.debug = kwargs.pop("debug", False)
        if "log_cmds" in kwargs or self.log_cmds is None:
            self.log_cmds = kwargs.pop("log_cmds", False)
        if "log_cmdout" in kwargs or self.log_cmdout is None:
            self.log_cmdout = kwargs.pop("log_cmdout", False)
        if "report" in kwargs or self.report is None:
            self.report = kwargs.pop("report", Report())

        # Default log configuration for FSL wrapper commands
        if self.fsllog is None:
            self.fsllog = kwargs.pop("fsllog", None)
            if not self.fsllog:
                self.fsllog = {"stderr" : self.log}
                if self.debug or self.log_cmds:
                    self.fsllog.update({"cmd" : self.log})
                if self.debug or self.log_cmdout:
                    self.fsllog.update({"stdout" : self.log})

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

    def __getattribute__(self, name):
        ret = super(Workspace, self).__getattribute__(name)
        if isinstance(ret, ImageProxy):
            return ret.img()
        else:
            return ret

    def __getattr__(self, name):
        ret = None
        if name not in self._search_childs:
            for wsp in self._search_childs:
                default_wsp = getattr(self, wsp)
                if isinstance(default_wsp, Workspace) and default_wsp != self:
                    val = getattr(default_wsp, name)
                    if val is not None:
                        ret = val
                        break

        if ret is None and self._parent is not None:
            ret = getattr(self._parent, name)

        return ret

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
        :param save: If False, do not save item
        :param save_name: If specified, alternative name to use for saving this item
        :param save_fn: If specified, Callable which generates string representation of
                        value suitable for saving the item to a file
        """
        if save:
            if not save_name:
                save_name = name

            # Remove any existing file first - it could be left behind if
            # the extension is different or the new value is None
            existing_files = glob.glob(os.path.join(self.savedir, "%s.*" % save_name))
            if not isinstance(value, Workspace):
                existing_files += glob.glob(os.path.join(self.savedir, save_name))
            for existing_file in existing_files:
                if os.path.isdir(existing_file):
                    shutil.rmtree(existing_file)
                else:
                    os.remove(existing_file)

            if value is not None:
                if save_fn is not None:
                    with open(os.path.join(self.savedir, save_name), "w") as tfile:
                        tfile.write(save_fn(value))
                elif isinstance(value, Image):
                    # Save as Nifti file. Note need to work around fslpy lack of respect for umask
                    fname = os.path.join(self.savedir, save_name)
                    value.save(fname)
                    os.chmod(value.dataSource, 0o666 ^ current_umask())
                    value.name = save_name
                    # Replace images with ImageProxy objects to avoid excess in-memory storage
                    if isinstance(value, AslImage):
                        value = AslImageProxy(fname, md=dict(value.metaItems()))
                    elif isinstance(value, Image):
                        value = ImageProxy(fname, md=dict(value.metaItems()))

                elif isinstance(value, np.ndarray) and value.ndim in (1, 2):
                    # Save as ASCII matrix
                    with open(os.path.join(self.savedir, save_name + ".mat"), "w") as tfile:
                        tfile.write(matrix_to_text(value))
                elif not name.startswith("_") and isinstance(value, pd.DataFrame):
                    # Save data frame in CSV file
                    value.to_csv(os.path.join(self.savedir, save_name + ".csv"), 
                                 index=False, header=True, float_format='%.4g')
                elif not name.startswith("_") and isinstance(value, (int, float, six.string_types)):
                    # Save other attributes in JSON file
                    self._stuff[name] = value
                    self._save_stuff()

        super(Workspace, self).__setattr__(name, value)

    def sub(self, name, parent_default=True, **kwargs):
        """
        Create a sub-workspace, (i.e. a subdir of this workspace)

        This inherits the log configuration from the parent workspace. The savedir
        will be a subdirectory of the original workspace. Additional data may be
        set using keyword arguments. The child-workspace will be available
        as an attribute on the parent workspace.

        :param name: Name of sub workspace
        :param parent_default: If True, attribute values default to the parent workspace
                               if not set on the sub-workspace
        """
        savedir = os.path.join(self.savedir, name)
        if parent_default and name not in self._search_childs:
            parent = self
        else:
            parent = None
            kwargs["log"] = self.log
            kwargs["debug"] = self.debug
            kwargs["report"] = self.report

        # We don't want sub-workspaces to have default search children unless
        # explicitly asked for FIXME probably better if always needs to be
        # explicit
        if "search_childs" not in kwargs:
            kwargs["search_childs"] = []
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
    # Pad 1D array to 2D matrix
    if mat.ndim == 1:
        mat = np.array([mat,])
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
        if not vals:
            continue
        # Check correct number of columns
        if ncols < 0:
            ncols = len(vals)
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
            if fail_if_exists:
                raise
            elif warn_if_exists:
                log.write("WARNING: mkdir - Directory %s already exists\n" % dirname)
    return os.path.abspath(dirname)
