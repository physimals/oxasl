"""
Tests for CALIB module

FIXME need satrecov tests
FIXME need sensitivity correction tests
FIXME need test with edge correction
"""
import math
import StringIO

import pytest
import numpy as np

from fsl.data.image import Image

from oxasl import Workspace, reg, AslImage

def get_wsp():
    wsp = Workspace()
    wsp.struc = Image(np.random.rand(10, 10, 10))
    wsp.asl.mean_brain = AslImage(np.random.rand(5, 5, 5, 4), tis=[1, 2], iaf="tc", ibf="rpt")
    wsp.calibration.brain = Image(np.random.rand(5, 5, 5))
    wsp.done("preproc_asl")
    wsp.done("calib.init")
    return wsp

def test_get_regfrom_supplied():
    """
    Test a pre-supplied regfrom is preferred to anything else
    """
    wsp = get_wsp()
    user_regfrom = Image(np.random.rand(5, 5, 5))
    wsp.reg.regfrom = user_regfrom
    reg.get_regfrom(wsp)
    assert(np.all(user_regfrom.data == wsp.reg.regfrom.data))

def test_get_regfrom_asldata_mean():
    """
    Test brain extracted ASL mean is used by default
    """
    wsp = get_wsp()
    reg.get_regfrom(wsp)
    assert(np.all(wsp.asl.mean_brain.data == wsp.reg.regfrom.data))

def test_get_regfrom_calib():
    """
    Test brain extracted calibration image is used if mean ASL data not available
    """
    wsp = get_wsp()
    wsp.asl.mean_brain = None
    reg.get_regfrom(wsp)
    assert(np.all(wsp.calibration.brain.data == wsp.reg.regfrom.data))
