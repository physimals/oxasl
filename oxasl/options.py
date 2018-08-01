"""
Classes for parsing command line options for ASL tools

Copyright (c) 2008-2018 University of Oxford
"""
from optparse import OptionGroup, OptionParser, Option, OptionValueError
from collections import defaultdict
from copy import copy

from fsl.data.image import Image

from ._version import __version__

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

    def add_category(self, category):
        """ 
        Add an OptionCategory to the parser
        """
        for group in category.groups(self):
            self.add_option_group(group)
            for option in group.option_list:
                self._categories[category.name].append(option)

    def filter(self, options, category, consume=True):
        """
        Filter options, returning only those in a specified category

        :param options: Options dictionary or namespace
        :param category: Category name
        :param consume: If True, remove filtered options from the passed dictionary
        
        :return: Dictionary of options that were found in the specified category
        """
        if not isinstance(options, dict):
            options = vars(options)
            
        filtered = {}
        for option in self._categories[category]:
            if option.dest in options:
                filtered[option.dest] = options[option.dest]
                if consume:
                    options.pop(option.dest)
        return filtered

class OptionCategory(object):
    """
    A named category of options. 
    """
    def __init__(self, name, ignore=()):
        self.name = name
        self.ignore = ignore
    
    def groups(self, parser):
        """
        :param parser: OptionParser instance
        :return: Sequence of OptionGroup instances for this category of options
        """
        return []

class IgnorableOptionGroup(OptionGroup):
    """
    OptionGroup with support for ignoring certain options
    """
    def __init__(self, *args, **kwargs):
        """
        Create AslOptionGroup

        :param ignore: Sequence of option names/destinations to ignore
        """
        self._ignore = kwargs.pop("ignore", [])
        OptionGroup.__init__(self, *args, **kwargs)

    def add_option(self, *args, **kwargs):
        """
        Add option - overridden to skip options we have been asked to ignore
        """
        name = args[0]
        if name not in self._ignore and name.lstrip("-") not in self._ignore and ("dest" not in kwargs or kwargs["dest"] not in self._ignore):
            OptionGroup.add_option(self, *args, **kwargs)

class GenericOptions(OptionCategory):
    """
    OptionCategory which contains generic options common to many command line tools
    """

    def __init__(self, title="Generic", output_type="directory", **kwargs):
        OptionCategory.__init__(self, "generic", **kwargs)
        self.title = title
        self.output_type = output_type

    def groups(self, parser):
        group = IgnorableOptionGroup(parser, self.title, ignore=self.ignore)
        group.add_option("-o", "--output", dest="output", help="Output %s" % self.output_type, default=None)
        group.add_option("-m", "--mask", dest="mask", help="Brain mask image in native ASL space", default=None, type="image")
        group.add_option("--debug", help="Debug mode", action="store_true", default=False)
        return [group, ]

def _check_image(option, opt, value):
    try:
        return Image(value)
    except ValueError:
        raise OptionValueError("option %s: invalid Image value: %r" % (opt, value))
            
class _ImageOption(Option):
    TYPES = Option.TYPES + ("image",)
    TYPE_CHECKER = copy(Option.TYPE_CHECKER)
    TYPE_CHECKER["image"] = _check_image
