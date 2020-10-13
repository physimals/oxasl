#!/bin/env python
"""
OXASL - Prefiltering module

This is run before any model fitting and is intended for volume selection, denoising etc.

Copyright (c) 2008-2020 Univerisity of Oxford
"""

try:
    import oxasl_enable
except ImportError:
    oxasl_enable = None
    
def run(wsp):
    if oxasl_enable and wsp.use_enable:
        oxasl_enable.init()
        #wsp.sub("enable")
        oxasl_enable.run(wsp)
        #wsp.corrected.asldata = wsp.enable.asldata_enable
