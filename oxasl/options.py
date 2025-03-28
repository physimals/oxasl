"""
Classes for parsing command line options for ASL tools

Copyright (c) 2008 University of Nottingham
"""
import sys
from optparse import OptionGroup, OptionParser, Option, OptionValueError
from collections import defaultdict
from copy import copy
import tempfile

import numpy as np

from fsl.data.image import Image

from oxasl import __version__

class AslOptionParser(OptionParser):
    """
    OptionParser which is extended to include the concept of option categories

    A category contains one or more OptionGroup instances. An option
    dictionary can be filtered by category, i.e. to extract the options
    belonging to a particular category (e.g. 'calibration options'). The
    existing OptionGroup class is not sufficient for this because groups
    cannot be nested
    """
    def __init__(self, usage="", version=__version__, **kwargs):
        OptionParser.__init__(self, usage=usage, version=version, option_class=_ImageOption, **kwargs)
        self._categories = defaultdict(list)

    def parse_args(self, argv=None, values=None):
        if argv is None:
            argv = sys.argv[1:]
        options, args = OptionParser.parse_args(self, argv, values)
        if options.optfile:
            # When an option file is specifeid, extract the options, build
            # a new argv vector and re-parse it. This is the only way to ensure
            # that options in the file work identically to CLI options.
            new_argv = self._add_from_file(argv, options.optfile)
            options, args = OptionParser.parse_args(self, new_argv, values)

        # Deal with case where asldata is given as separate files
        if args and options.asldata is None:
            merged_data = None
            for idx, fname in enumerate(args):
                img = Image(fname)
                shape = list(img.shape)
                if img.ndim == 3:
                    shape += [1,]
                if merged_data is None:
                    merged_data = np.zeros(shape[:3] + [shape[3] * len(args)])
                merged_data[..., idx*shape[3]:(idx+1)*shape[3]] = img.data
            merged_img = Image(merged_data, header=img.header)
            temp_asldata = tempfile.NamedTemporaryFile(prefix="oxasl", delete=True)
            options.asldata = temp_asldata.name
            merged_img.save(options.asldata)

        return options, args

    def _add_from_file(self, argv, optfile):
        new_argv = list(argv)
        with open(optfile, "r") as f:
            optlines = f.readlines()

        for line in optlines:
            line = line[:line.find("#")].strip()
            if line:
                line = line.lstrip("-").replace(":", " ").replace("=", " ")
                kv = line.split(None, 1)
                key = "-" + kv[0]
                if len(kv[0]) > 1:
                    key = "-" + key
                new_argv.append(key)
                if len(kv) == 2:
                    new_argv.append(kv[1])
        return new_argv

    def add_category(self, category):
        """
        Add an OptionCategory to the parser
        """
        for group in category.groups(self):
            self.add_option_group(group)
            for option in group.option_list:
                self._categories[category.name].append(option)

class OptionCategory(object):
    """
    A named category of options.
    """
    def __init__(self, name):
        self.name = name

    def groups(self, parser):
        """
        :param parser: OptionParser instance
        :return: Sequence of OptionGroup instances for this category of options
        """
        return []

class GenericOptions(OptionCategory):
    """
    OptionCategory which contains generic options common to many command line tools
    """

    def __init__(self, title="Generic", output_type="directory", default_output="oxasl", **kwargs):
        OptionCategory.__init__(self, "generic", **kwargs)
        self.title = title
        self.output_type = output_type
        self.default_output = default_output

    def groups(self, parser):
        group = OptionGroup(parser, self.title)
        group.add_option("--output", "-o", help="Output %s" % self.output_type, default=self.default_output)
        group.add_option("--overwrite", help="Overwrite output %s if it already exists" % self.output_type, action="store_true", default=False)
        group.add_option("--mask", "-m", help="Brain mask image in ASL space", default=None, type="image")
        group.add_option("--optfile", help="File containing additional options")
        group.add_option("--log-cmds", help="Log all external commands run", action="store_true", default=False)
        group.add_option("--log-cmdout", help="Log the standard output of all external commands run", action="store_true", default=False)
        group.add_option("--debug", help="Debug mode - log all command output and keep all output files", action="store_true", default=False)
        return [group, ]

def load_options_file(fname):
    options = {}
    if fname:
        with open(fname) as options_file:
            for line in options_file:
                keyval = line.strip().rstrip("\n").lstrip("--").split("=", 1)
                key = keyval[0].strip()
                if key != "":
                    if len(keyval) == 1:
                        options[key] = True
                    else:
                        options[key] = keyval[1].strip()
    return options

def load_matrix(fname):
    """
    Load a matrix from a text file

    Matrix should be a series of ASCII floating point numbers separated by whitespace
    on each line.

    :param fname: Filename
    """
    matrix = []
    with open(fname, "r") as f:
        for line in f.readlines():
            if line.strip():
                matrix.append([float(v) for v in line.strip().split()])
    return np.array(matrix, dtype=np.float32)

def _check_image(option, opt, value):
    try:
        return Image(value, loadData=False)
    except ValueError:
        raise OptionValueError("option %s: invalid Image value: %r" % (opt, value))

def _check_matrix(option, opt, value):
    try:
        return load_matrix(value)
    except ValueError:
        raise OptionValueError("option %s: invalid matrix value: %r" % (opt, value))

def _check_optfile(option, opt, value):
    try:
        return load_options_file(value)
    except ValueError:
        raise OptionValueError("option %s: invalid options file: %r" % (opt, value))

class _ImageOption(Option):
    TYPES = Option.TYPES + ("image", "matrix", "optfile",)
    TYPE_CHECKER = copy(Option.TYPE_CHECKER)
    TYPE_CHECKER["image"] = _check_image
    TYPE_CHECKER["matrix"] = _check_matrix
    TYPE_CHECKER["optfile"] = _check_optfile
