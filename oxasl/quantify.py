#!/bin/env python
"""
OXASL - Quantification module

Fits a kinetic model to differenced ASL data

Copyright (c) 2008-2020 Univerisity of Oxford
"""

from oxasl import basil

try:
    import oxasl_multite
except ImportError:
    oxasl_multite = None

try:
    import oxasl_ve
except ImportError:
    oxasl_ve = None

def run(wsp):
    wsp.quantify_wsps = []
    quantify = _get_quantify_method(wsp)
    wsp.log.write("\nQuantifying perfusion from ASL data\n")
    quantify(wsp)

def _already_quantified(wsp):
    # Make fake basil-like directory for already quantified data
    wsp.log.write(" - Skipping quantification, data is already quantified\n")
    quantify_wsp = wsp.sub("prequant")
    # Two possible locations for compatibility
    if wsp.rois is not None and wsp.rois.mask is not None:
        quantify_wsp.analysis_mask = wsp.rois.mask
    else:
        quantify_wsp.analysis_mask = wsp.mask
    finaldir = quantify_wsp.sub("finalstep")
    finaldir.mean_ftiss = wsp.asldata
    wsp.quantify_wsps.append("prequant")
    return quantify_wsp

def _get_quantify_method(wsp):
    if wsp.asldata.iaf == "quant":
        return _already_quantified
    if wsp.asldata.iaf in ("tc", "ct", "diff"):
        if wsp.asldata.ntes == 1:
            return basil.run
        elif oxasl_multite is None:
            raise ValueError("Multi-TE data supplied but oxasl_multite is not installed")
        else:
            return oxasl_multite.run
    elif wsp.asldata.iaf == "ve":
        if oxasl_ve is None:
            raise ValueError("VE data supplied but oxasl_ve is not installed")
        else:
            return oxasl_ve.run
