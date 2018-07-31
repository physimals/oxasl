"""
Workspace for executing FSL commands
"""

import os
import sys
import errno
import tempfile

import numpy as np

from fsl.data.image import Image

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

    def __init__(self, savedir=None, **kwargs):
        """
        Create workspace

        :param savedir: If specified, use this path to save data. Will be created
                        if it does not not already exist
        :param log:     File stream to write log output to (default: sys.stdout)
        """
        if savedir is not None:
            self._savedir = os.path.abspath(savedir)
            mkdir(savedir, log=kwargs.get("log", sys.stdout))
        else:
            self._savedir = None

        # Defaults might be overridden by kwargs
        self.log = sys.stdout
        self.debug = False
        self.fsllog = None
        for key, value in kwargs.items():
            setattr(self, key, value)
        
        # Default log configuration for FSL wrapper commands
        if not self.fsllog:
            self.fsllog = {"stderr" : self.log}
            if self.debug:
                self.fsllog.update({"stdout" : self.log, "cmd" : self.log})

    def __getattr__(self, name):
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

    def set_item(self, name, value, save=True, save_name=None):
        """
        Add an item to the workspace

        The item will be set as an attribute on the workspace
        If a save directory is configured and the value is a supported
        type it will be saved there. This can be disabled by setting save=False

        :param name: Name, must be a valid Python identifier
        :param value: Value to set
        :param save: If False, do not save item even if savedir defined
        :param save: If specified, alternative name to use for saving this item
        """
        if value is not None and save and self._savedir:
            if not save_name:
                save_name = name
            if isinstance(value, Image):
                # Save as Nifti file
                value.save(os.path.join(self._savedir, save_name))
                value.name = save_name
            elif isinstance(value, np.ndarray) and value.ndim == 2:
                # Save as ASCII matrix
                with open(os.path.join(self._savedir, save_name), "w") as tfile:
                    tfile.write(matrix_to_text(value))
        super(Workspace, self).__setattr__(name, value)

    def sub(self, name, **kwargs):
        """
        Create a sub-workspace, (i.e. a subdir of this workspace)

        This inherits the log configuration from the parent workspace. If
        a savedir is defined, it will be created with a savedir which is
        a subdirectory of the original workspace. Additional data may be 
        set using keyword arguments.
        
        :param name: Name of subdir
        """
        if self._savedir:
            savedir = os.path.join(self._savedir, name)
        else:
            savedir = None
        setattr(self, name, Workspace(savedir=savedir, log=self.log, debug=self.debug, **kwargs))
        return getattr(self, name)

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

def tempdir(suffix, debug=False, log=sys.stdout):
    """
    Create a temporary directory

    :param debug: If True, creates directory in current working directory
    """
    if debug:
        tmpdir = os.path.join(os.getcwd(), "tmp_%s" % suffix)
        mkdir(tmpdir, log=log)
    else:
        tmpdir = tempfile.mkdtemp("_%s" % suffix)
    return tmpdir
