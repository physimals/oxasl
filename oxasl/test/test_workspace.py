"""
Tests for workspace module
"""
import os
import shutil
import tempfile
from six import StringIO

import numpy as np
import pytest

from fsl.data.image import Image

from oxasl import Workspace, AslImage
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
    log = StringIO()
    wsp = Workspace(log=log)
    wsp.log.write("hello")
    assert(log.getvalue() == "hello")

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

def test_sub_inherit():
    """ Test sub workspaces can inherit values from their parent """
    wsp = Workspace()
    wsp.wibble = 7
    wsp.wobble = 6
    wsp.sub("child")
    wsp.child.wobble = 5
    assert(wsp.child.wibble == 7)
    assert(wsp.child.wobble == 5)

def test_sub_inherit_wsp():
    """ Test sub workspaces can inherit sub-workspaces from their parent """
    wsp = Workspace()
    wsp.sub("child1")
    wsp.child1.wibble = 7
    wsp.sub("child2")
    assert(wsp.child2.child1 is not None)
    assert(wsp.child2.child1.wibble == 7)

def test_input_wsp():
    """ Test putting constructor attributes in a default sub workspaces """
    wsp = Workspace(input_wsp="cakes", flapjack=4, fruit=3, defaults=[])
    assert(wsp.cakes is not None)
    assert(wsp.cakes.flapjack == 4)
    assert(wsp.cakes.fruit == 3)
    assert(wsp.flapjack is None)
    assert(wsp.fruit is None)

def test_default_wsp():
    """ Test default sub-workspaces for search """
    wsp = Workspace(defaults=["cars"])
    assert(wsp.cars is None)
    wsp.ferrari = 9
    wsp.merc = 8
    wsp.sub("cars")
    wsp.cars.porsche = 6
    wsp.cars.ferrari = 4
    assert(wsp.cars is not None)
    assert(wsp.ferrari == 9)
    assert(wsp.porsche == 6)
    assert(wsp.merc == 8)
    assert(wsp.cars.porsche == 6)
    assert(wsp.cars.ferrari == 4)
    assert(wsp.cars.merc is None)
    
def test_default_wsp_multiple():
    """ Test multiple default sub-workspaces for search """
    wsp = Workspace(defaults=["plants", "trees"])
    wsp.daffodil = 9
    wsp.larch = 1
    wsp.sub("trees")
    wsp.trees.oak = 3
    wsp.trees.larch = 2
    wsp.trees.apple = 7
    assert(wsp.daffodil == 9)
    assert(wsp.larch == 1)
    assert(wsp.oak == 3)
    assert(wsp.apple == 7)
    assert(wsp.trees.larch == 2)
    assert(wsp.trees.oak == 3)
    assert(wsp.trees.daffodil is None)
    assert(wsp.trees.apple == 7)
    wsp.sub("plants")
    wsp.plants.lily = 4
    wsp.plants.oak = 5
    assert(wsp.daffodil == 9)
    assert(wsp.larch == 1)
    assert(wsp.lily == 4)
    assert(wsp.oak == 5)
    assert(wsp.apple == 7)
    assert(wsp.trees.oak == 3)
    assert(wsp.trees.lily is None)
    assert(wsp.plants.daffodil is None)
    assert(wsp.plants.lily == 4)
    assert(wsp.plants.oak == 5)
    
def test_savedir_created():
    """ Test save dirs are created if they don't already exist """
    tempdir = tempfile.mktemp("_oxasl")
    try:
        log = StringIO()
        wsp = Workspace(savedir=tempdir, log=log)
        assert(wsp.savedir) == tempdir
        assert(os.path.isdir(tempdir))
        assert("WARNING" not in log.getvalue())
    finally:
        shutil.rmtree(tempdir)

def test_savedir_created_multilevel():
    """ Test multi-level save dirs are created if they don't already exist """
    tempdir = os.path.join(tempfile.mktemp("_oxasl"), "extra", "levels")
    try:
        log = StringIO()
        wsp = Workspace(savedir=tempdir, log=log)
        assert(wsp.savedir) == tempdir
        assert(os.path.isdir(tempdir))
        assert("WARNING" not in log.getvalue())
    finally:
        shutil.rmtree(tempdir)

def test_savedir_sub():
    """ Test sub-workspace have subdirs created """
    tempdir = tempfile.mktemp("_oxasl")
    try:
        log = StringIO()
        wsp = Workspace(savedir=tempdir, log=log)
        wsp.sub("quark")
        path = os.path.join(tempdir, "quark")
        assert(wsp.quark.savedir == path)
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

def test_matrix_save_name():
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

def _custom_save(mat):
    return "Custom Save"

def test_custom_save():
    """ 
    Test matrices are saved in the savedir with the specified name
    """
    tempdir = tempfile.mkdtemp("_oxasl")
    try:
        wsp = Workspace(savedir=tempdir)
        mat = np.random.rand(4, 4)
        wsp.set_item("testmat", mat, save_fn=_custom_save)
        path = os.path.join(tempdir, "testmat")
        assert(os.path.exists(path))
        with open(path) as sfile:
            assert("Custom Save" == sfile.read())
    finally:
        shutil.rmtree(tempdir)

def test_custom_save_name():
    """ 
    Test matrices are saved in the savedir with the specified name
    """
    tempdir = tempfile.mkdtemp("_oxasl")
    try:
        wsp = Workspace(savedir=tempdir)
        mat = np.random.rand(4, 4)
        wsp.set_item("testmat", mat, save_name="potato", save_fn=_custom_save)
        path = os.path.join(tempdir, "testmat")
        assert(not os.path.exists(path))
        path = os.path.join(tempdir, "potato")
        assert(os.path.exists(path))
        with open(path) as sfile:
            assert("Custom Save" == sfile.read())
    finally:
        shutil.rmtree(tempdir)

def test_custom_save_nosave():
    """ 
    Test matrices are saved in the savedir with the specified name
    """
    tempdir = tempfile.mkdtemp("_oxasl")
    try:
        wsp = Workspace(savedir=tempdir)
        mat = np.random.rand(4, 4)
        wsp.set_item("testmat", mat, save_fn=_custom_save, save=False)
        path = os.path.join(tempdir, "testmat")
        assert(not os.path.exists(path))
    finally:
        shutil.rmtree(tempdir)

def test_savedir_already_exists():
    """ 
    Test warning when save dir already exists
    """
    tempdir = tempfile.mkdtemp("_oxasl")
    try:
        log = StringIO()
        wsp = Workspace(savedir=tempdir, log=log)
        assert("WARNING" in log.getvalue())
        assert("already exists" in log.getvalue())
    finally:
        shutil.rmtree(tempdir)

def test_fsllog_default():
    """
    Test the FSL logging context created
    """
    log = StringIO()
    wsp = Workspace(log=log)
    assert(isinstance(wsp.fsllog, dict))
    assert(wsp.fsllog.get("stdout", None) is None)
    assert(wsp.fsllog.get("stderr", None) == log)
    assert(wsp.fsllog.get("cmd", None) is None)

def test_fsllog_debug():
    """
    Test the FSL logging context created in debug mode
    """
    log = StringIO()
    wsp = Workspace(debug=True, log=log)
    assert(isinstance(wsp.fsllog, dict))
    assert(wsp.fsllog.get("stdout", None) == log)
    assert(wsp.fsllog.get("stderr", None) == log)
    assert(wsp.fsllog.get("cmd", None) == log)

def test_aslimage():
    kwargs = {
        "asldata" : np.random.rand(5, 5, 5, 8),
        "tis" : [1, 2],
        "iaf" : "tc",
        "ibf" : "rpt",
    }
    wsp = Workspace(auto_asldata=True, **kwargs)
    assert(isinstance(wsp.asldata, AslImage))
    assert(wsp.asldata.tis == [1, 2])
    assert(wsp.asldata.iaf == "tc")
    assert(wsp.asldata.order == "ltr")
    assert(wsp.asldata.rpts == [2, 2])
    
def test_aslimage_missing():
    with pytest.raises(ValueError):
        Workspace(auto_asldata=True)
    with pytest.raises(ValueError):
        Workspace(auto_asldata=True, asldata=None)

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
