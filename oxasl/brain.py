#!/bin/env python
"""
Brain extraction

Copyright (c) 2008-2013 Univerisity of Oxford
"""

import fsl.wrappers as fsl

def brain(wsp, img, thresh=0.5):
    bet_result = fsl.bet(img, seg=True, mask=False, fracintensity=thresh, output=fsl.LOAD, log=wsp.fsllog)
    return bet_result["output"]

def mask(wsp, img, thresh):
    bet_result = fsl.bet(img, seg=False, mask=True, fracintensity=thresh, output=fsl.LOAD, log=wsp.fsllog)
    return bet_result["output_mask"]
