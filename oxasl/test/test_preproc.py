import numpy as np

from oxasl import Workspace, AslImage, preproc

def test_preproc_none():
    wsp = Workspace()
    d = np.random.rand(5, 5, 5, 6)
    wsp.asldata = AslImage(d, name="asldata", tis=[1.5], iaf="tc", order="lrt")

    preproc.preprocess(wsp)
    assert wsp.asldata.ntis == wsp.asldata_preproc.ntis
    assert wsp.asldata.tis == wsp.asldata_preproc.tis
    assert wsp.asldata.have_plds == wsp.asldata_preproc.have_plds
    assert wsp.asldata.rpts == wsp.asldata_preproc.rpts
    assert wsp.asldata.ntc == wsp.asldata_preproc.ntc
    assert wsp.asldata.order == wsp.asldata_preproc.order

def test_preproc_diff():
    wsp = Workspace()
    d = np.random.rand(5, 5, 5, 6)
    wsp.asldata = AslImage(d, name="asldata", tis=[1.5], iaf="tc", order="lrt")
    wsp.diff = True

    preproc.preprocess(wsp)
    assert wsp.asldata_preproc.ntis == wsp.asldata.ntis
    assert wsp.asldata_preproc.tis == wsp.asldata.tis
    assert wsp.asldata_preproc.have_plds == wsp.asldata.have_plds
    assert wsp.asldata_preproc.rpts == wsp.asldata.rpts
    assert wsp.asldata_preproc.ntc == 1
    assert wsp.asldata_preproc.order == "rt"

def test_preproc_reorder_diff():
    wsp = Workspace()
    d = np.random.rand(5, 5, 5, 6)
    wsp.asldata = AslImage(d, name="asldata", tis=[1.5], iaf="tc", order="lrt")
    wsp.diff = True
    wsp.reorder = "ltr"

    preproc.preprocess(wsp)
    assert wsp.asldata_preproc.ntis == wsp.asldata.ntis
    assert wsp.asldata_preproc.tis == wsp.asldata.tis
    assert wsp.asldata_preproc.have_plds == wsp.asldata.have_plds
    assert wsp.asldata_preproc.rpts == wsp.asldata.rpts
    assert wsp.asldata_preproc.ntc == 1
    assert wsp.asldata_preproc.order == "tr"

def test_preproc_diff_reorder():
    wsp = Workspace()
    d = np.random.rand(5, 5, 5, 6)
    wsp.asldata = AslImage(d, name="asldata", tis=[1.5], iaf="tc", order="lrt")
    wsp.diff = True
    wsp.reorder = "ltr"
    
    preproc.preprocess(wsp)
    assert wsp.asldata_preproc.ntis == wsp.asldata.ntis
    assert wsp.asldata_preproc.tis == wsp.asldata.tis
    assert wsp.asldata_preproc.have_plds == wsp.asldata.have_plds
    assert wsp.asldata_preproc.rpts == wsp.asldata.rpts
    assert wsp.asldata_preproc.ntc == 1
    assert wsp.asldata_preproc.order == "tr"

def test_preproc_smooth():
    wsp = Workspace()
    d = np.random.rand(5, 5, 5, 6)
    wsp.asldata = AslImage(d, name="asldata", tis=[1.5], iaf="tc", order="lrt")
    wsp.smooth = True
    wsp.fwhm = 1.5
    
    preproc.preprocess(wsp)
    assert wsp.asldata_preproc.ntis == wsp.asldata.ntis
    assert wsp.asldata_preproc.tis == wsp.asldata.tis
    assert wsp.asldata_preproc.have_plds == wsp.asldata.have_plds
    assert wsp.asldata_preproc.rpts == wsp.asldata.rpts
    assert wsp.asldata_preproc.ntc == wsp.asldata.ntc
    assert wsp.asldata_preproc.order == wsp.asldata.order

def test_preproc_moco():
    wsp = Workspace()
    d = np.random.rand(5, 5, 5, 6)
    wsp.asldata = AslImage(d, name="asldata", tis=[1.5], iaf="tc", order="lrt")
    wsp.mc = True
    
    preproc.preprocess(wsp)
    assert wsp.asldata_preproc.ntis == wsp.asldata.ntis
    assert wsp.asldata_preproc.tis == wsp.asldata.tis
    assert wsp.asldata_preproc.have_plds == wsp.asldata.have_plds
    assert wsp.asldata_preproc.rpts == wsp.asldata.rpts
    assert wsp.asldata_preproc.ntc == wsp.asldata.ntc
    assert wsp.asldata_preproc.order == wsp.asldata.order
