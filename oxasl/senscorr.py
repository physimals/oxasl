#!/bin/env python
"""
OXASL - Sensitivity correction module

Copyright (c) 2008-2020 Univerisity of Oxford
"""
from __future__ import unicode_literals

import numpy as np

from fsl.data.image import Image

from oxasl import reg
from oxasl.options import OptionCategory, OptionGroup
from oxasl.reporting import LightboxImage

class Options(OptionCategory):
    """
    Options for corrections of the input data
    """

    def __init__(self):
        OptionCategory.__init__(self, "corrections")

    def groups(self, parser):
        ret = []

        g = OptionGroup(parser, "Sensitivity correction")
        g.add_option("--cref", help="Reference image for sensitivity correction", type="image")
        g.add_option("--cact", help="Image from coil used for actual ASL acquisition (default: calibration image - only in longtr mode)", type="image")
        g.add_option("--isen", help="User-supplied sensitivity correction in ASL space")
        g.add_option("--senscorr-auto", "--senscorr", help="Apply automatic sensitivity correction using bias field from FAST", action="store_true", default=False)
        g.add_option("--senscorr-off", help="Do not apply any sensitivity correction", action="store_true", default=False)
        ret.append(g)

        return ret

def run(wsp):
    """
    Get sensitivity correction image

    Required workspace attributes
    -----------------------------

     - ``asldata`` : ASL data

    Optional workspace attributes
    -----------------------------

     - ``isen`` : User supplied sensitivity image
     - ``cact`` : Calibration image. Used in conjunction with ``cref`` to calculate sensitivity map
     - ``calib`` : Calibration image. Used as alternative to cact provided ``mode`` is ``longtr``
     - ``cref`` : Calibration reference image
     - ``senscorr_auto`` : If True, automatically calculate sensitivity correction using FAST
     - ``senscorr_off`` If True, do not apply sensitivity correction

    Updated workspace attributes
    ----------------------------

     - ``sensitivity``    : Sensitivity correction image in ASL space
    """
    wsp.log.write("\nCalculating Sensitivity correction\n")
    sensitivity = None
    bias = None
    if wsp.senscorr_off:
        wsp.log.write(" - Sensitivity correction disabled\n")
    elif wsp.isen is not None:
        wsp.log.write(" - Sensitivity image supplied by user\n")
        sensitivity = wsp.isen
    elif wsp.cact is not None and wsp.cref is not None:
        wsp.log.write(" - Sensitivity image calculated from calibration actual and reference images\n")
        cref_data = np.copy(wsp.cref.data)
        cref_data[cref_data == 0] = 1
        sensitivity = Image(wsp.cact.data.astype(np.float32) / cref_data, header=wsp.calib.header)
    elif wsp.calib is not None and wsp.cref is not None:
        if wsp.ifnone("mode", "longtr") != "longtr":
            raise ValueError("Calibration reference image specified but calibration image was not in longtr mode - need to provided additional calibration image using the ASL coil")
        wsp.log.write(" - Sensitivity image calculated from calibration and reference images\n")
        cref_data = np.copy(wsp.cref.data)
        cref_data[cref_data == 0] = 1
        sensitivity = Image(wsp.calib.data.astype(np.float32) / cref_data, header=wsp.calib.header)
    elif wsp.senscorr_auto and wsp.structural.bias is not None:
        wsp.log.write(" - Sensitivity image calculated from bias field\n")
        bias = reg.change_space(wsp, wsp.structural.bias, "asl")
        sensitivity = Image(np.reciprocal(bias.data), header=bias.header)
    else:
        wsp.log.write(" - No source of sensitivity correction was found\n")

    if sensitivity is not None:
        sdata = sensitivity.data
        sdata[sdata < 1e-12] = 1
        sdata[np.isnan(sdata)] = 1
        sdata[np.isinf(sdata)] = 1
        wsp.sub("senscorr")
        wsp.senscorr.sensitivity = Image(sdata, header=sensitivity.header)

        page = wsp.report.page("sensitivity")
        page.heading("Sensitivity correction", level=0)
        page.heading("Sensitivity map", level=1)
        page.image("sensitivity", LightboxImage(wsp.senscorr.sensitivity))

    if bias is not None:
        wsp.senscorr.bias = bias
