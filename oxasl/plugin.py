"""
Plugin system for OXASL
=======================

Design
~~~~~~

A plugin is registered as a Python entry point. It must be an instance of
OxaslPlugin and provide the following information:

 - Name
 - Brief description
 - Version
 - Optional OptionCategory defining command line options supported
 - Optional list of IAF values which the plugin supports modelling for
 - Optional PVC estimation method?
"""

from oxasl.options import OptionCategory

class OxaslPlugin:
    def __init__(self, name, description, version, **kwargs):
        self.name = name
        self.description = description
        self.version = version
        self.options = kwargs.get("options", OptionCategory())
