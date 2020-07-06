#!/bin/env python
"""
Brain extraction

Copyright (c) 2008-2013 Univerisity of Oxford
"""

import fsl.wrappers as fsl

def brain(wsp, img, thresh=0.5):
    """
    Extract brain from wholehead image
    """
    bet_result = fsl.bet(img, seg=True, mask=False, fracintensity=thresh, output=fsl.LOAD, log=wsp.fsllog)
    return bet_result["output"]

def mask(wsp, img, thresh):
    """
    Extract brain mask from wholehead image
    """
    bet_result = fsl.bet(img, seg=False, mask=True, fracintensity=thresh, output=fsl.LOAD, log=wsp.fsllog)
    return bet_result["output_mask"]
