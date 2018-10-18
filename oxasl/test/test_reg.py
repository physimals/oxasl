"""
Tests for CALIB module

FIXME need satrecov tests
FIXME need sensitivity correction tests
FIXME need test with edge correction
"""
import math

import pytest
import numpy as np

from fsl.data.image import Image

from oxasl import Workspace, reg, AslImage, brain

def get_wsp():
    wsp = Workspace(debug=True)
    wsp.struc = Image(np.random.rand(10, 10, 10))
    wsp.calib = Image(np.random.rand(5, 5, 5))
    return wsp

def test_get_regfrom_supplied():
    """
    Test a pre-supplied regfrom is preferred to anything else
    """
    wsp = get_wsp()
    user_regfrom = Image(np.random.rand(5, 5, 5))
    wsp.regfrom = user_regfrom
    reg.get_regfrom(wsp)
    assert(np.allclose(user_regfrom.data, wsp.reg.regfrom.data))

def test_get_regfrom_asldata_mean_tc():
    """
    Test brain extracted ASL mean is used for TC data
    """
    wsp = get_wsp()
    wsp.asldata = AslImage(np.random.rand(5, 5, 5, 4), tis=[1, 2], iaf="tc", ibf="rpt")
    reg.get_regfrom(wsp)
    meanasl_brain = brain.brain(wsp, wsp.asldata.mean(), thresh=0.2)
    assert(np.allclose(meanasl_brain.data, wsp.reg.regfrom.data))

def test_get_regfrom_asldata_mean_ct():
    """
    Test brain extracted ASL mean is used for CT data
    """
    wsp = get_wsp()
    wsp.asldata = AslImage(np.random.rand(5, 5, 5, 4), tis=[1, 2], iaf="ct", ibf="rpt")
    reg.get_regfrom(wsp)
    meanasl_brain = brain.brain(wsp, wsp.asldata.mean(), thresh=0.2)
    assert(np.allclose(meanasl_brain.data, wsp.reg.regfrom.data))

def test_get_regfrom_calib():
    """
    Test brain extracted calibration image is used for differenced data
    """
    wsp = get_wsp()
    wsp.asldata = AslImage(np.random.rand(5, 5, 5, 4), tis=[1, 2], iaf="diff", ibf="rpt")
    reg.get_regfrom(wsp)
    calib_brain = brain.brain(wsp, wsp.calib, thresh=0.2)
    assert(np.allclose(calib_brain.data, wsp.reg.regfrom.data))
