import numpy as np

from oxasl import AslImage, preproc

def test_preproc_none():
    d = np.random.rand(5, 5, 5, 6)
    orig = AslImage(d, name="asldata", tis=[1.5], order="lrt")
    orig.summary()

    proc = preprocess(orig)
    proc.summary()
    assert orig.name == proc.name
    assert orig.ntis == proc.ntis
    assert orig.tis == proc.tis
    assert orig.have_plds == proc.have_plds
    assert orig.rpts == proc.rpts
    assert orig.ntc == proc.ntc
    assert orig.order == proc.order

def test_preproc_diff():
    d = np.random.rand(5, 5, 5, 6)
    orig = AslImage(d, name="asldata", tis=[1.5], order="lrt")

    proc = preprocess(orig, diff=True)
    assert proc.ntis == orig.ntis
    assert proc.tis == orig.tis
    assert proc.have_plds == orig.have_plds
    assert proc.rpts == orig.rpts
    assert proc.name == orig.name + "_diff"
    assert proc.ntc == 1
    assert proc.order == "rt"

def test_preproc_reorder_diff():
    d = np.random.rand(5, 5, 5, 6)
    orig = AslImage(d, name="asldata", tis=[1.5], order="lrt")
    
    proc = preprocess(orig, diff=True, reorder="ltr")
    assert proc.ntis == orig.ntis
    assert proc.tis == orig.tis
    assert proc.have_plds == orig.have_plds
    assert proc.rpts == orig.rpts
    assert proc.name == orig.name + "_diff_reorder"
    assert proc.ntc == 1
    assert proc.order == "tr"

def test_preproc_diff_reorder():
    d = np.random.rand(5, 5, 5, 6)
    orig = AslImage(d, name="asldata", tis=[1.5], order="lrt")
    
    proc = preprocess(orig, diff=True, reorder="tr")
    assert proc.ntis == orig.ntis
    assert proc.tis == orig.tis
    assert proc.have_plds == orig.have_plds
    assert proc.rpts == orig.rpts
    assert proc.name == orig.name + "_diff_reorder"
    assert proc.ntc == 1
    assert proc.order == "tr"

def test_preproc_smooth():
    d = np.random.rand(5, 5, 5, 6)
    orig = AslImage(d, name="asldata", tis=[1.5], order="lrt")
    
    proc = preprocess(orig, smooth=True, fwhm=1.5)
    assert proc.ntis == orig.ntis
    assert proc.tis == orig.tis
    assert proc.have_plds == orig.have_plds
    assert proc.rpts == orig.rpts
    assert proc.name == orig.name + "_smooth"
    assert proc.ntc == orig.ntc
    assert proc.order == orig.order

def test_preproc_moco():
    d = np.random.rand(5, 5, 5, 6)
    orig = AslImage(d, name="asldata", tis=[1.5], order="lrt")
    
    proc = preprocess(orig, mc=True)
    assert proc.ntis == orig.ntis
    assert proc.tis == orig.tis
    assert proc.have_plds == orig.have_plds
    assert proc.rpts == orig.rpts
    assert proc.name == orig.name + "_mc"
    assert proc.ntc == orig.ntc
    assert proc.order == orig.order
