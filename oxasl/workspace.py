"""
Workspace for executing FSL commands
"""

import os
import sys
import shutil
import errno
import tempfile

import numpy as np

from fsl.data.image import Image

from .fabber import Fabber, percent_progress

class Workspace(object):
    """
    A workspace for FSL processing

    Based on a working directory, which may be temporary or persistent.
    """

    def __init__(self, workdir=None, log=sys.stdout, imgs=(), path=None, debug=False, echo=False, use_local_dir=True):
        """
        Create workspace

        :param workdir: If specified, use this path for the working directory. Will be created
                        if it does not not already exist
        :param log:     File stream to write log output to (default: sys.stdout)
        :param imgs:    List of Image object which will be save to working directory
        :param path:    Optional list of directories to search for binaries. If not specified, will
                        look in $FSLDEVDIR/bin, $FSLDIR/bin
        :param debug:   If True, enable debugging messages
        :param echo:    If True, echo commands and their output to the log
        :param use_local_dir: If True, use the directory of the calling program (sys.argv[0])
        as the first entry in the binary search path. This is enabled by default because
        it gives a useful way for scripts to replace standard FSL programs with their own copies
        by putting them in the same directory as the script.
        """
        if workdir is not None:
            self.workdir = os.path.abspath(workdir)
            mkdir(workdir)
            self.is_temp = False
        else:
            self.workdir = tempdir("fsl")
            self.is_temp = True
        
        for img in imgs:
            self.add_img(img)

        if path is None:
            path = []
            if use_local_dir:
                path.append(os.path.dirname(os.path.abspath(sys.argv[0])))
            for env_var in ("FSLDEVDIR", "FSLDIR"):
                if env_var in os.environ:
                    path.append(os.path.join(os.environ[env_var], "bin"))
        self._path = path
        self.log = log
        self.debug_enabled = debug
        self.echo_enabled = echo
        
        if "FSLOUTPUTTYPE" not in os.environ:
            os.environ["FSLOUTPUTTYPE"] = "NIFTI_GZ"

    def __del__(self):
        if self.is_temp:
            shutil.rmtree(self.workdir)

    def img(self, name):
        """
        Get an fsl.data.image.Image from the workspace
        """
        return Image(os.path.join(self.workdir, name))

    def add_img(self, img, name=None):
        """
        Save an image to the workspace

        :param img: Image object. It will be saved in the workspace and its filepath modified
        """
        if name is None:
            name = img.name
        img.save(os.path.join(self.workdir, name))
        return img

    def add_file(self, fpath, name=None):
        """
        Add a random named file to the workspace
        """
        if name is None:
            name = os.path.basename(fpath)
        shutil.copyfile(fpath, os.path.join(self.workdir, name))

    def add_text(self, text, name):
        """
        Write a text file to the workspace

        :param text: Text to put in the file
        :param name: File name
        """
        with open(os.path.join(self.workdir, name), "w") as tfile:
            tfile.write(text)

    def del_img(self, img):
        """
        Delete an image from the workspace

        FIXME not working
        """
        try:
            image, _ = self._input_img(img)
            os.remove(image.fpath)
        except IOError:
            self.log.write("WARNING: failed to delete %s\n" % img)

    def sub(self, name, imgs=()):
        """
        Create a sub-workspace, (i.e. a subdir of this workspace)

        This inherits the log output from the parent workspace

        :param name: Name of subdir
        :imgs imgs: Images to copy to the sub workspace
        """
        return Workspace(os.path.join(self.workdir, name), imgs=imgs, log=self.log)

    def _wsp_fn(method):
        def _func(self, *args, **kwargs):
            cwd = os.getcwd()
            try:
                os.chdir(self.workdir)
                return method(self, *args, **kwargs)
            finally:
                os.chdir(cwd)
        return _func

    @_wsp_fn
    def fabber(self, options, **kwargs):
        """
        Fabber model fitting

        :param options: Fabber options
        """
        options = dict(options)
        output_name = options.pop("output", "fabber_output")

        # Replace fsl.Image objects with the underlying data
        for key in options.keys():
            value = options[key]
            if isinstance(value, Image):
                options[key] = value.nibImage.get_data()

        fab = Fabber()
        run = fab.run(options, percent_progress)
        run.write_to_dir(os.path.join(self.workdir, output_name))
        return run

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

def mkdir(dirname, fail_if_exists=False, warn_if_exists=True):
    """
    Create a directory, including necessary subdirs
    """
    try:
        os.makedirs(dirname)
    except OSError as exc:
        if exc.errno == errno.EEXIST:
            if fail_if_exists: raise
            elif warn_if_exists: print("WARNING: mkdir - Directory %s already exists" % dirname)
    return os.path.abspath(dirname)

def tempdir(suffix, debug=False):
    """
    Create a temporary directory

    :param debug: If True, creates directory in current working directory
    """
    if debug:
        tmpdir = os.path.join(os.getcwd(), "tmp_%s" % suffix)
        mkdir(tmpdir)
    else:
        tmpdir = tempfile.mkdtemp("_%s" % suffix)
    return tmpdir
