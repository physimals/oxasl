#!/bin/env python
"""
OXASL - Module to calibrate a perfusion output using previously calculated M0 value or image

Copyright (c) 2008 University of Nottingham
"""
import numpy as np

from fsl.data.image import Image

from oxasl.options import OptionCategory

class Options(OptionCategory):
    """
    Options for calibration
    """

    def __init__(self):
        OptionCategory.__init__(self, "calibrate")

    def groups(self, parser):
        groups = []
        return groups

def run(wsp, perf_img, multiplier=1.0, var=False):
    """
    Do calibration of a perfusion image from a calibration (M0) image or value

    :param wsp: Workspace object
    :param perf_img: Image containing perfusion data to calibrate
    :param multiplier: Scalar multiple to convert output to physical units
    :param alpha: Inversion efficiency
    :param var: If True, assume data represents variance rather than value

    :return: Image containing calibrated data

    Required workspace attributes
    -----------------------------

     - ``m0`` : M0 single value or voxelwise Image
    """
    if not perf_img:
        raise ValueError("Perfusion data cannot be None")
    if not wsp.m0:
        raise ValueError("No calibration data supplied")

    wsp.log.write("\nCalibrating perfusion data: %s\n" % perf_img.name)
    alpha = wsp.ifnone("calib_alpha", 1.0 if wsp.asldata.iaf in ("ve", "vediff") else 0.85 if wsp.asldata.casl else 0.98)
    m0 = wsp.m0
    if isinstance(m0, Image):
        m0 = m0.data

    if var:
        wsp.log.write(" - Treating data as variance - squaring M0 correction, multiplier and inversion efficiency\n")
        m0 = m0**2
        multiplier = multiplier**2
        alpha = alpha**2

    if isinstance(m0, np.ndarray):
        # If M0 is <= zero, make calibrated data zero
        calibrated = np.zeros(perf_img.shape)
        calibrated[m0 > 0] = perf_img.data[m0 > 0] / m0[m0 > 0]
    else:
        calibrated = perf_img.data / m0

    if alpha != 1.0:
        wsp.log.write(" - Using inversion efficiency correction: %f\n" % alpha)
        calibrated /= alpha

    if multiplier != 1.0:
        wsp.log.write(" - Using multiplier for physical units: %f\n" % multiplier)
        calibrated *= multiplier

    perf_calib = Image(calibrated, name=perf_img.name, header=perf_img.header)
    return perf_calib
