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
    wsp.asldata_mean_brain = AslImage(np.random.rand(5, 5, 5, 4), tis=[1, 2], iaf="tc", ibf="rpt")
    wsp.calib_brain = Image(np.random.rand(5, 5, 5))
    wsp.done("preproc_asl")
    wsp.done("preproc_calib")
    return wsp

def test_get_regfrom_supplied():
    """
    Test a pre-supplied regfrom is preferred to anything else
    """
    wsp = get_wsp()
    user_regfrom = Image(np.random.rand(5, 5, 5))
    wsp.regfrom = user_regfrom
    reg.get_regfrom(wsp)
    assert(np.all(user_regfrom.data == wsp.regfrom.data))

def test_get_regfrom_asldata_mean():
    """
    Test brain extracted ASL mean is used by default
    """
    wsp = get_wsp()
    reg.get_regfrom(wsp)
    assert(np.all(wsp.asldata_mean_brain.data == wsp.regfrom.data))

def test_get_regfrom_calib():
    """
    Test brain extracted calibration image is used if mean ASL data not available
    """
    wsp = get_wsp()
    wsp.asldata_mean_brain = None
    reg.get_regfrom(wsp)
    assert(np.all(wsp.calib_brain.data == wsp.regfrom.data))
