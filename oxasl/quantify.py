#!/bin/env python
"""
OXASL -Quantification module

Fits a kinetic model to differenced ASL data

Copyright (c) 2008-2020 Univerisity of Oxford
"""

from oxasl import basil, reg, output

try:
    import oxasl_multite
except ImportError:
    oxasl_multite = None

def run(wsp):
    wsp.basildirs = []

    # Quantification in native space
    quantify = _get_quantify_method(wsp)

    quantify.run(wsp.sub("basil"))
    wsp.basildirs.append("")

    # Re-do registration using PWI as reference
    reg.run(wsp, redo=True, struc_bbr=True, struc_flirt=False)

    # Quantification in alternate spaces
    for quantify_space in ("struc", "std", "custom"):
        if wsp.ifnone("quantify_%s" % quantify_space, False):
            basil_wsp = wsp.sub("basil_%s" % quantify_space)
            wsp.basildirs.append(quantify_space)
            basil_wsp.image_space = quantify_space
            quantify.run(basil_wsp) 

def _get_quantify_method(wsp):
    if wsp.asldata.iaf in ("tc", "ct", "diff"):
        if wsp.asldata.ntes == 1:
            return basil
        elif oxasl_multite is None:
            raise ValueError("Multi-TE data supplied but oxasl_multite is not installed")
        else:
            return oxasl_multite
