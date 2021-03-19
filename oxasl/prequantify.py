#!/bin/env python
"""
OXASL - Pre-quantification module

This is run before model fitting and should generate effectively differenced ASL data
e.g. by multiphase decoding

Copyright (c) 2008-2020 Univerisity of Oxford
"""

try:
    import oxasl_ve
except ImportError:
    oxasl_ve = None
    
try:
    import oxasl_mp
except ImportError:
    oxasl_mp = None
    
def run(wsp):
    if wsp.asldata.iaf in ("ve", "vediff"):
        if oxasl_ve is None:
            raise ValueError("Vessel encoded data supplied but oxasl_ve is not installed")
        oxasl_ve.run(wsp)
    elif wsp.asldata.iaf == "mp":
        if oxasl_mp is None:
            raise ValueError("Multiphase data supplied but oxasl_mp is not installed")
        oxasl_mp.run(wsp)
