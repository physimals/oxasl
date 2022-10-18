#!/bin/env python
"""
OXASL -Quantification module

Fits a kinetic model to differenced ASL data

Copyright (c) 2008-2020 Univerisity of Oxford
"""

from oxasl import basil, reg

try:
    import oxasl_multite
except ImportError:
    oxasl_multite = None

try:
    import oxasl_ve
except ImportError:
    oxasl_ve = None

def run(wsp):
    wsp.basildirs = []

    # Quantification in ASL space
    quantify = _get_quantify_method(wsp)
    quantify(wsp)

    # Re-do registration using PWI as reference
    reg.run(wsp, redo=True, struc_bbr=True, struc_flirt=False, use_basil_wsp=wsp.basil)

    # Quantification in alternate spaces
    for basildir in wsp.basildirs:
        for quantify_space in ("struc", "std", "custom"):
            if wsp.ifnone("quantify_%s" % quantify_space, False):
                if basildir:
                    quantify_name = "%s_%s" % (basildir, quantify_space)
                else:
                    quantify_name = quantify_space
                basil_wsp = wsp.sub("basil_ " + quantify_name)
                wsp.basildirs.append(quantify_name)
                basil_wsp.image_space = quantify_space
                quantify.run(basil_wsp) 

def _default_quantify(wsp):
    basil.run(wsp.sub("basil"))
    wsp.basildirs.append("")

def _already_quantified(wsp):
    wsp.log.write(" - Skipping quantification, data is already quantified\n")
    basildir = wsp.sub("basil")
    # Two possible locations for compatibility
    if wsp.rois is not None and wsp.rois.mask is not None:
        basildir.analysis_mask = wsp.rois.mask
    else:
        basildir.analysis_mask = wsp.mask
    finaldir = basildir.sub("finalstep")
    finaldir.mean_ftiss = wsp.asldata
    wsp.basildirs.append("")

def _get_quantify_method(wsp):
    if wsp.asldata.iaf == "quant":
        return _already_quantified
    if wsp.asldata.iaf in ("tc", "ct", "diff"):
        if wsp.asldata.ntes == 1:
            return _default_quantify
        elif oxasl_multite is None:
            raise ValueError("Multi-TE data supplied but oxasl_multite is not installed")
        else:
            return oxasl_multite
    elif wsp.asldata.iaf == "ve":
        if oxasl_ve is None:
            raise ValueError("VE data supplied but oxasl_ve is not installed")
        else:
            return oxasl_ve.run

