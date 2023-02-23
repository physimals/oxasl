#!/bin/env python
"""
OXASL - Prefiltering module

This is run before any model fitting and is intended for volume selection, denoising etc.
Pre-filtering processes do not depend on the type of ASL data (e.g. multiphase, vessel-encoded)

Copyright (c) 2008-2020 Univerisity of Oxford
"""

try:
    import oxasl_enable
except ImportError:
    oxasl_enable = None

try:
    import oxasl_deblur
except ImportError:
    oxasl_deblur = None

def run(wsp):
    wsp.sub("filter")

    if oxasl_enable and wsp.use_enable:
        oxasl_enable.run(wsp, output_wsp=wsp.filter)

    if oxasl_deblur and wsp.deblur:
        oxasl_deblur.run(wsp, output_wsp=wsp.filter)
