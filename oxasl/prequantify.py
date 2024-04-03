#!/bin/env python
"""
OXASL - Pre-quantification module

This is run before model fitting and should generate effectively differenced ASL data
e.g. by multiphase decoding

Copyright (c) 2008 University of Nottingham
"""

try:
    import oxasl_mp
except ImportError:
    oxasl_mp = None
    
def run(wsp):
    if wsp.asldata.iaf == "mp":
        if oxasl_mp is None:
            raise ValueError("Multiphase data supplied but oxasl_mp is not installed")
        wsp.sub("prequantify")
        oxasl_mp.run(wsp, output_wsp=wsp.prequantify)
