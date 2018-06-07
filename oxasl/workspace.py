#!/usr/bin/env python
"""
FSL wrapper functions - named fslwrap so as not to clash with any potential
'official' FSL python modules which would probably be named 'fsl'
"""

import os
import sys
import shutil
import shlex
import subprocess
import errno
import tempfile
import collections
import re

import nibabel as nib
import numpy as np

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
        :param imgs: List of Image object which will be save to working directory
        :param path: Optional list of directories to search for binaries. If not specified, will
        look in $FSLDEVDIR/bin, $FSLDIR/bin
        :param debug: If True, enable debugging messages
        :param echo: If True, echo commands and their output to the log
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

    def fabber(self, model_group, img, mask, options, overwrite=True, **kwargs):
        """
        Fabber model fitting

        :param img: Input image
        :param mask: Mask image
        :return: MVN output image
        """
        options = dict(options)

        img, itype = self._input_img(img)
        output_name, extra_args = self._get_std(img, "_fabber", kwargs)
        options["data"] = img.ipath

        if mask:
            mask, _ = self._input_img(mask)
            options["mask"] = mask.ipath
            
        options["output"] = output_name
        options["save-mvn"] = ""
        if overwrite:
            options["overwrite"] = ""
        option_args = ""
        for k in sorted(options.keys()):
            v = options[k]
            if v == "" or (v and type(v) == bool):
                option_args += " --%s" % k
            elif v:
                option_args += " --%s=%s" % (k, str(v))

        args = "%s %s" % (option_args, extra_args)
        imgs, _ = self.run("fabber_%s" % model_group.lower(), args=args, expected=[output_name + "/finalMVN"])
        return self._output_img(imgs[0], itype)
        
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
    except OSError as e:
        if e.errno == errno.EEXIST:
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