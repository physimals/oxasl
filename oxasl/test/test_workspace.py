"""
Tests for workspace module
"""
#import StringIO

#import pytest
#import numpy as np

from oxasl import Workspace

def test_default_attr():
    """ Check attributes are None by default """
    wsp = Workspace()
    assert(wsp.wibble is None)

def test_set_attr():
    """ Check attributes are None by default """
    wsp = Workspace()
    assert(wsp.wibble is None)
    wsp.wibble = 7
    assert(wsp.wibble == 7)
