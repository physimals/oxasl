"""
Tests for workspace module
"""
import os
import shutil
import tempfile
import StringIO

import numpy as np
import pytest

from fsl.data.image import Image

from oxasl import Workspace
from oxasl.workspace import text_to_matrix

def test_default_attr():
    """ Check attributes are None by default """
    wsp = Workspace()
    assert(wsp.wibble is None)

def test_set_attr():
    """ Check attributes can bet set """
    wsp = Workspace()
    assert(wsp.wibble is None)
    wsp.wibble = 7
    assert(wsp.wibble == 7)

def test_ctor_attributes():
    """ Check attributes specified in constructor """
    wsp = Workspace(wobble="hi")
    assert(wsp.wobble == "hi")

def test_log():
    """ Check that the log is picked up """
    log = StringIO.StringIO()
    wsp = Workspace(log=log)
    wsp.log.write("hello")
    assert(log.getvalue() == "hello")

def test_done():
    """ Check the 'done' logic """
    wsp = Workspace()
    assert(not wsp.isdone("process"))
    wsp.done("process")
    assert(wsp.isdone("process"))

def test_undone():
    """ Check the 'done' logic can be undone """
    wsp = Workspace()
    assert(not wsp.isdone("process"))
    wsp.done("process")
    assert(wsp.isdone("process"))
    wsp.done("process", status=False)
    assert(not wsp.isdone("process"))

def test_ifnone():
    wsp = Workspace(wibble=11)
    assert(wsp.ifnone("wibble", 12) == 11)
    assert(wsp.ifnone("wobble", 12) == 12)

def test_sub():
    """ Test sub-workspaces """
    wsp = Workspace()
    wsp.sub("child")
    assert(isinstance(wsp.child, Workspace))
    assert(wsp.child.wibble is None)
    assert(wsp.child.log == wsp.log)
    
def test_sub_kwargs():
    """ Test creating a sub workspace with kwargs """
    wsp = Workspace()
    wsp.sub("child", wibble="squid", pudding=4)
    assert(isinstance(wsp.child, Workspace))
    assert(wsp.child.wibble == "squid")
    assert(wsp.child.pudding == 4)

def test_savedir_created():
    """ Test save dirs are created if they don't already exist """
    tempdir = tempfile.mktemp("_oxasl")
    try:
        log = StringIO.StringIO()
        wsp = Workspace(savedir=tempdir, log=log)
        assert(wsp._savedir) == tempdir
        assert(os.path.isdir(tempdir))
        assert("WARNING" not in log.getvalue())
    finally:
        shutil.rmtree(tempdir)

def test_savedir_created_multilevel():
    """ Test multi-level save dirs are created if they don't already exist """
    tempdir = os.path.join(tempfile.mktemp("_oxasl"), "extra", "levels")
    try:
        log = StringIO.StringIO()
        wsp = Workspace(savedir=tempdir, log=log)
        assert(wsp._savedir) == tempdir
        assert(os.path.isdir(tempdir))
        assert("WARNING" not in log.getvalue())
    finally:
        shutil.rmtree(tempdir)

def test_savedir_sub():
    """ Test sub-workspace have subdirs created """
    tempdir = tempfile.mktemp("_oxasl")
    try:
        log = StringIO.StringIO()
        wsp = Workspace(savedir=tempdir, log=log)
        wsp.sub("quark")
        path = os.path.join(tempdir, "quark")
        assert(wsp.quark._savedir == path)
        assert(os.path.isdir(path))
        assert("WARNING" not in log.getvalue())
    finally:
        shutil.rmtree(tempdir)

def test_image_save():
    """ 
    Test images are saved in the savedir
    """
    tempdir = tempfile.mkdtemp("_oxasl")
    try:
        wsp = Workspace(savedir=tempdir)
        img = Image(np.random.rand(5, 5, 5))
        wsp.testimg = img
        path = os.path.join(tempdir, "testimg.nii.gz")
        assert(os.path.isfile(path))
        otherimg = Image(path)
        assert(np.all(img.data == wsp.testimg.data))
        assert(np.all(img.data == otherimg.data))
    finally:
        shutil.rmtree(tempdir)

def test_image_nosave():
    """ 
    Test setting an image without saving
    """
    tempdir = tempfile.mkdtemp("_oxasl")
    try:
        wsp = Workspace(savedir=tempdir)
        img = Image(np.random.rand(5, 5, 5))
        wsp.set_item("testimg", img, save=False)
        path = os.path.join(tempdir, "testimg.nii.gz")
        assert(not os.path.exists(path))
        assert(np.all(img.data == wsp.testimg.data))
    finally:
        shutil.rmtree(tempdir)

def test_image_save_name():
    """ 
    Test images are saved in the savedir with the specified name
    """
    tempdir = tempfile.mkdtemp("_oxasl")
    try:
        wsp = Workspace(savedir=tempdir)
        img = Image(np.random.rand(5, 5, 5))
        wsp.set_item("testimg", img, save_name="pumpkin")
        path = os.path.join(tempdir, "testimg.nii.gz")
        assert(not os.path.exists(path))
        path = os.path.join(tempdir, "pumpkin.nii.gz")
        assert(os.path.isfile(path))
        otherimg = Image(path)
        assert(np.all(img.data == wsp.testimg.data))
        assert(np.all(img.data == otherimg.data))
    finally:
        shutil.rmtree(tempdir)

def test_matrix_save():
    """ 
    Test 2D matrices are saved in the savedir
    """
    tempdir = tempfile.mkdtemp("_oxasl")
    try:
        wsp = Workspace(savedir=tempdir)
        mat = np.random.rand(4, 4)
        wsp.testmat = mat
        path = os.path.join(tempdir, "testmat")
        assert(os.path.isfile(path))
        with open(path) as matfile:
            othermat = text_to_matrix(matfile.read())
        assert(np.all(mat == wsp.testmat))
        assert(np.all(mat == othermat))
    finally:
        shutil.rmtree(tempdir)

def test_matrix_nosave():
    """ 
    Test setting an matrix without saving
    """
    tempdir = tempfile.mkdtemp("_oxasl")
    try:
        wsp = Workspace(savedir=tempdir)
        mat = np.random.rand(4, 4)
        wsp.set_item("testmat", mat, save=False)
        path = os.path.join(tempdir, "testmat")
        assert(not os.path.exists(path))
        assert(np.all(mat == wsp.testmat))
    finally:
        shutil.rmtree(tempdir)

def test_image_save_name():
    """ 
    Test matrices are saved in the savedir with the specified name
    """
    tempdir = tempfile.mkdtemp("_oxasl")
    try:
        wsp = Workspace(savedir=tempdir)
        mat = np.random.rand(4, 4)
        wsp.set_item("testmat", mat, save_name="parsnip")
        path = os.path.join(tempdir, "testmat")
        assert(not os.path.exists(path))
        path = os.path.join(tempdir, "parsnip")
        assert(os.path.isfile(path))
        with open(path) as matfile:
            othermat = text_to_matrix(matfile.read())
        assert(np.all(mat == wsp.testmat))
        assert(np.all(mat == othermat))
    finally:
        shutil.rmtree(tempdir)

def test_savedir_already_exists():
    """ 
    Test warning when save dir already exists
    """
    tempdir = tempfile.mkdtemp("_oxasl")
    try:
        log = StringIO.StringIO()
        wsp = Workspace(savedir=tempdir, log=log)
        assert("WARNING" in log.getvalue())
        assert("already exists" in log.getvalue())
    finally:
        shutil.rmtree(tempdir)

def test_fsllog_default():
    """
    Test the FSL logging context created
    """
    log = StringIO.StringIO()
    wsp = Workspace(log=log)
    assert(isinstance(wsp.fsllog, dict))
    assert(wsp.fsllog.get("stdout", None) is None)
    assert(wsp.fsllog.get("stderr", None) == log)
    assert(wsp.fsllog.get("cmd", None) is None)

def test_fsllog_debug():
    """
    Test the FSL logging context created in debug mode
    """
    log = StringIO.StringIO()
    wsp = Workspace(debug=True, log=log)
    assert(isinstance(wsp.fsllog, dict))
    assert(wsp.fsllog.get("stdout", None) == log)
    assert(wsp.fsllog.get("stderr", None) == log)
    assert(wsp.fsllog.get("cmd", None) == log)

def test_text_to_matrix_spaces():
    """
    Check that text_to_matrix works with space separated data
    """
    text = "1 2 3\n4 5 6\n"
    mat = text_to_matrix(text)
    assert(np.all(mat == [[1, 2, 3], [4, 5, 6]]))

def test_text_to_matrix_comma():
    """
    Check that text_to_matrix works with comma separated data
    """
    text = "1, 2, 3\n4,5,6\n"
    mat = text_to_matrix(text)
    assert(np.all(mat == [[1, 2, 3], [4, 5, 6]]))
    
def test_text_to_matrix_tabs():
    """
    Check that text_to_matrix works with tab separated data
    """
    text = "1\t2\t3\n4\t 5\t 6\n"
    mat = text_to_matrix(text)
    assert(np.all(mat == [[1, 2, 3], [4, 5, 6]]))
        
def test_text_to_matrix_mixed():
    """
    Check that text_to_matrix works with mixed separators
    """
    text = "1\t2   3\n4  ,  5, \t 6\n"
    mat = text_to_matrix(text)
    assert(np.all(mat == [[1, 2, 3], [4, 5, 6]]))
    
def test_text_to_matrix_not_matrix():
    text = "1 2 3\n4 5\n"
    with pytest.raises(ValueError):
        mat = text_to_matrix(text)
    
def test_text_to_matrix_not_numbers():
    text = "1 x 3\n4 5 6\n"
    with pytest.raises(ValueError):
        mat = text_to_matrix(text)
