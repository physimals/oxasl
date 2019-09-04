"""
Tests for CALIB module

FIXME need satrecov tests
FIXME need sensitivity correction tests
FIXME need test with edge correction
"""
import math
from six import StringIO

import pytest
import numpy as np

from fsl.data.image import Image

from oxasl import Workspace, calib, AslImage

def _get_imgs(shape=(5, 5, 5)):
    """
    Create test images
    """
    perf_d = np.random.rand(5, 5, 5)
    perf_img = Image(name="perfusion", image=perf_d)

    calib_d = np.random.rand(5, 5, 5)
    calib_img = Image(name="calib", image=calib_d)

    return perf_img, calib_img

def _expected_m0(ref_mean, t1, t2, pc, gain=1.0, alpha=1.0, te=0, tr=3.2, t2b=150):
    """
    Expected M0 value from reference region mean
    """
    #print("ref mean: ", ref_mean)
    m0 =  ref_mean / (1 - math.exp(- tr / t1))
    #print("T1 correction: ", m0)
    m0 = m0 / math.exp(-te/t2)
    m0 = m0 * gain / pc
    m0 = m0 * math.exp(-te/t2b)
    #print("T2 correction: ", m0)
    m0 *= alpha
    #print("IE correction: ", m0)
    return m0

def test_nodata():
    """
    Check we get an error if there is no perfusion data
    """
    d = np.random.rand(5, 5, 5, 6)
    calib_img = Image(name="calib", image=d)

    wsp = Workspace(calib=calib_img, calib_method="voxelwise")
    with pytest.raises(ValueError):
        calib.calibrate(wsp, None)

def test_nocalib():
    """
    Check we get an error if there is no calibration data
    """
    d = np.random.rand(5, 5, 5, 6)
    perf_img = Image(name="perfusion", image=d)

    wsp = Workspace(calib_method="voxelwise")
    with pytest.raises(ValueError):
        calib.calibrate(wsp, perf_img)

def test_bad_method():
    """
    Check we get an error on an invalid method
    """
    d = np.random.rand(5, 5, 5, 6)
    perf_img = Image(name="perfusion", image=d)

    d = np.random.rand(5, 5, 5, 6)
    calib_img = Image(name="calib", image=d)

    wsp = Workspace(calib=calib_img, calib_method="random")
    with pytest.raises(ValueError):
        calib.calibrate(wsp, perf_img)

def test_defaults():
    """
    Check default calibration works
    """
    perf_d = np.random.rand(5, 5, 5)
    perf_img = Image(name="perfusion", image=perf_d)

    calib_d = np.random.rand(5, 5, 5)
    calib_img = Image(name="calib", image=calib_d)

    wsp = Workspace(calib=calib_img, calib_method="voxelwise")
    perf_calib = calib.calibrate(wsp, perf_img)
    calibrated_d = perf_calib.data
    assert(perf_calib.name == "perfusion_calib")
    assert(perf_calib.shape == perf_img.shape)

    # Default partition coefficient is 0.9
    np.testing.assert_allclose(calibrated_d, 0.9 * perf_d / calib_d)

def test_cgain():
    """
    Check calibration gain
    """
    GAIN = 1.23
    perf_d = np.random.rand(5, 5, 5)
    perf_img = Image(name="perfusion", image=perf_d)

    calib_d = np.random.rand(5, 5, 5)
    calib_img = Image(name="calib", image=calib_d)

    wsp = Workspace(calib=calib_img, calib_method="voxelwise", calib_gain=GAIN)
    perf_calib = calib.calibrate(wsp, perf_img)
    calibrated_d = perf_calib.data
    assert(perf_calib.name == "perfusion_calib")
    assert(perf_calib.shape == perf_img.shape)

    # Default partition coefficient is 0.9
    np.testing.assert_allclose(calibrated_d, 0.9 / GAIN * perf_d / calib_d)

def test_alpha():
    """
    Check inversion efficiency
    """
    ALPHA = 0.74
    perf_d = np.random.rand(5, 5, 5)
    perf_img = Image(name="perfusion", image=perf_d)

    calib_d = np.random.rand(5, 5, 5)
    calib_img = Image(name="calib", image=calib_d)

    wsp = Workspace(calib=calib_img, calib_method="voxelwise")
    perf_calib = calib.calibrate(wsp, perf_img, alpha=ALPHA)
    calibrated_d = perf_calib.data
    assert(perf_calib.name == "perfusion_calib")
    assert(perf_calib.shape == perf_img.shape)

    # Default partition coefficient is 0.9
    np.testing.assert_allclose(calibrated_d, 0.9 / ALPHA * perf_d / calib_d)

def test_multiplier():
    """
    Check end multiplier for physical units
    """
    MULTIPLIER = 7654
    perf_d = np.random.rand(5, 5, 5)
    perf_img = Image(name="perfusion", image=perf_d)

    calib_d = np.random.rand(5, 5, 5)
    calib_img = Image(name="calib", image=calib_d)

    wsp = Workspace(calib=calib_img, calib_method="voxelwise")
    perf_calib = calib.calibrate(wsp, perf_img, multiplier=MULTIPLIER)
    calibrated_d = perf_calib.data
    assert(perf_calib.name == "perfusion_calib")
    assert(perf_calib.shape == perf_img.shape)

    # Default partition coefficient is 0.9
    np.testing.assert_allclose(calibrated_d, 0.9 * MULTIPLIER * perf_d / calib_d)

def test_partition_coeff():
    """
    Check tissue/blood partition coefficient
    """
    PC = 0.67
    perf_d = np.random.rand(5, 5, 5)
    perf_img = Image(name="perfusion", image=perf_d)

    calib_d = np.random.rand(5, 5, 5)
    calib_img = Image(name="calib", image=calib_d)

    wsp = Workspace(calib=calib_img, calib_method="voxelwise", pct=PC)
    perf_calib = calib.calibrate(wsp, perf_img)
    calibrated_d = perf_calib.data
    assert(perf_calib.name == "perfusion_calib")
    assert(perf_calib.shape == perf_img.shape)

    # Default partition coefficient is 0.9
    np.testing.assert_allclose(calibrated_d, PC * perf_d / calib_d)

def test_shorttr_corr():
    """
    Check correction for short TR
    """
    TR, T1 = 3.0, 1.1
    perf_d = np.random.rand(5, 5, 5)
    perf_img = Image(name="perfusion", image=perf_d)

    calib_d = np.random.rand(5, 5, 5)
    calib_img = Image(name="calib", image=calib_d)

    wsp = Workspace(calib=calib_img, calib_method="voxelwise", tr=TR, t1t=T1)
    perf_calib = calib.calibrate(wsp, perf_img)
    calibrated_d = perf_calib.data
    assert(perf_calib.name == "perfusion_calib")
    assert(perf_calib.shape == perf_img.shape)

    factor =  1 - math.exp(-TR / T1)
    # Default partition coefficient is 0.9
    np.testing.assert_allclose(calibrated_d, 0.9 * factor * perf_d / calib_d)

def test_shorttr_corr_not1():
    """
    Check correction for short TR is not done (and generates warning) if T1 not given
    """
    TR = 3.0
    perf_d = np.random.rand(5, 5, 5)
    perf_img = Image(name="perfusion", image=perf_d)

    calib_d = np.random.rand(5, 5, 5)
    calib_img = Image(name="calib", image=calib_d)

    log = StringIO()
    wsp = Workspace(calib=calib_img, calib_method="voxelwise", tr=TR, log=log)
    perf_calib = calib.calibrate(wsp, perf_img)
    calibrated_d = perf_calib.data
    assert(perf_calib.name == "perfusion_calib")
    assert(perf_calib.shape == perf_img.shape)
    assert("WARNING" in log.getvalue())

    # Default partition coefficient is 0.9
    np.testing.assert_allclose(calibrated_d, 0.9 * perf_d / calib_d)

def test_defaults_var():
    """
    Check default calibration works on variances
    """
    perf_d = np.random.rand(5, 5, 5)
    perf_img = Image(name="perfusion", image=perf_d)

    calib_d = np.random.rand(5, 5, 5)
    calib_img = Image(name="calib", image=calib_d)

    wsp = Workspace(calib=calib_img, calib_method="voxelwise")
    perf_calib = calib.calibrate(wsp, perf_img, var=True)
    calibrated_d = perf_calib.data
    assert(perf_calib.name == "perfusion_calib")
    assert(perf_calib.shape == perf_img.shape)

    # Default partition coefficient is 0.9
    np.testing.assert_allclose(calibrated_d, 0.9 * 0.9 * perf_d / np.square(calib_d))

def test_cgain_var():
    """
    Check calibration gain on variances
    """
    GAIN = 1.23
    perf_d = np.random.rand(5, 5, 5)
    perf_img = Image(name="perfusion", image=perf_d)

    calib_d = np.random.rand(5, 5, 5)
    calib_img = Image(name="calib", image=calib_d)

    wsp = Workspace(calib=calib_img, calib_method="voxelwise", calib_gain=GAIN)
    perf_calib = calib.calibrate(wsp, perf_img, var=True)
    calibrated_d = perf_calib.data
    assert(perf_calib.name == "perfusion_calib")
    assert(perf_calib.shape == perf_img.shape)

    # Default partition coefficient is 0.9
    np.testing.assert_allclose(calibrated_d, 0.9 * 0.9 / GAIN / GAIN * perf_d / calib_d / calib_d)

def test_alpha_var():
    """
    Check inversion efficiency on variances
    """
    ALPHA = 0.74
    perf_d = np.random.rand(5, 5, 5)
    perf_img = Image(name="perfusion", image=perf_d)

    calib_d = np.random.rand(5, 5, 5)
    calib_img = Image(name="calib", image=calib_d)

    wsp = Workspace(calib=calib_img, calib_method="voxelwise",)
    perf_calib = calib.calibrate(wsp, perf_img, var=True, alpha=ALPHA)
    calibrated_d = perf_calib.data
    assert(perf_calib.name == "perfusion_calib")
    assert(perf_calib.shape == perf_img.shape)

    # Default partition coefficient is 0.9
    np.testing.assert_allclose(calibrated_d, 0.9 * 0.9 / ALPHA /ALPHA * perf_d / calib_d / calib_d)

def test_multiplier_var():
    """
    Check end multiplier for physical units on variances
    """
    MULTIPLIER = 7654
    perf_d = np.random.rand(5, 5, 5)
    perf_img = Image(name="perfusion", image=perf_d)

    calib_d = np.random.rand(5, 5, 5)
    calib_img = Image(name="calib", image=calib_d)

    wsp = Workspace(calib=calib_img, calib_method="voxelwise")
    perf_calib = calib.calibrate(wsp, perf_img, var=True, multiplier=MULTIPLIER)
    calibrated_d = perf_calib.data
    assert(perf_calib.name == "perfusion_calib")
    assert(perf_calib.shape == perf_img.shape)

    # Default partition coefficient is 0.9
    np.testing.assert_allclose(calibrated_d, 0.9 * 0.9 * MULTIPLIER * MULTIPLIER * perf_d / calib_d / calib_d)

def test_partition_coeff_var():
    """
    Check tissue/blood partition coefficient on variances
    """
    PC = 0.67
    perf_d = np.random.rand(5, 5, 5)
    perf_img = Image(name="perfusion", image=perf_d)

    calib_d = np.random.rand(5, 5, 5)
    calib_img = Image(name="calib", image=calib_d)

    wsp = Workspace(calib=calib_img, calib_method="voxelwise", pct=PC)
    perf_calib = calib.calibrate(wsp, perf_img, var=True)
    calibrated_d = perf_calib.data
    assert(perf_calib.name == "perfusion_calib")
    assert(perf_calib.shape == perf_img.shape)

    # Default partition coefficient is 0.9
    np.testing.assert_allclose(calibrated_d, PC * PC * perf_d / calib_d / calib_d)

def test_shorttr_corr_var():
    """
    Check correction for short TR on variances
    """
    TR, T1 = 3.0, 1.1
    perf_d = np.random.rand(5, 5, 5)
    perf_img = Image(name="perfusion", image=perf_d)

    calib_d = np.random.rand(5, 5, 5)
    calib_img = Image(name="calib", image=calib_d)

    wsp = Workspace(calib=calib_img, calib_method="voxelwise", tr=TR, t1t=T1)
    perf_calib = calib.calibrate(wsp, perf_img, var=True)
    calibrated_d = perf_calib.data
    assert(perf_calib.name == "perfusion_calib")
    assert(perf_calib.shape == perf_img.shape)

    factor =  1 - math.exp(-TR / T1)
    # Default partition coefficient is 0.9
    np.testing.assert_allclose(calibrated_d, 0.9 * 0.9 * factor * factor * perf_d / calib_d / calib_d)

def test_refregion_defaults():
    """
    Reference region defaults to CSF
    """
    perf_img, calib_img = _get_imgs()

    ref_d = np.ones((5, 5, 5))
    ref_img = Image(name="refmask", image=ref_d)

    wsp = Workspace(calib=calib_img, calib_method="refregion", refmask=ref_img)
    perf_calib = calib.calibrate(wsp, perf_img)
    calibrated_d = perf_calib.data
    assert(perf_calib.name == "perfusion_calib")
    assert(perf_calib.shape == perf_img.shape)

    # CSF defaults
    m0_expected =  _expected_m0(np.mean(calib_img.data), 4.3, 750, 1.15)
    np.testing.assert_allclose(calibrated_d, perf_img.data / m0_expected)

def test_refregion_unknown():
    """
    Reference region unknown type
    """
    perf_img, calib_img = _get_imgs()
    
    ref_d = np.ones((5, 5, 5))
    ref_img = Image(name="refmask", image=ref_d)

    wsp = Workspace(calib=calib_img, calib_method="refregion", refmask=ref_img, tissref="kryptonite")
    with pytest.raises(ValueError):
        calib.calibrate(wsp, perf_img)
    
def test_refregion_csf_all():
    """
    Reference region mask is whole image
    """
    perf_img, calib_img = _get_imgs()
    
    ref_d = np.ones((5, 5, 5))
    ref_img = Image(name="refmask", image=ref_d)

    wsp = Workspace(calib=calib_img, calib_method="refregion", refmask=ref_img, tissref="csf")
    perf_calib = calib.calibrate(wsp, perf_img)
    calibrated_d = perf_calib.data

    # CSF defaults
    m0_expected =  _expected_m0(np.mean(calib_img.data), 4.3, 750, 1.15)
    np.testing.assert_allclose(calibrated_d, perf_img.data / m0_expected)

def test_refregion_wm_all():
    """
    Reference region mask is whole image
    """
    perf_img, calib_img = _get_imgs()

    ref_d = np.ones((5, 5, 5))
    ref_img = Image(name="refmask", image=ref_d)

    wsp = Workspace(calib=calib_img, calib_method="refregion", refmask=ref_img, tissref="wm")
    perf_calib = calib.calibrate(wsp, perf_img)
    calibrated_d = perf_calib.data

    # WM defaults
    m0_expected =  _expected_m0(np.mean(calib_img.data), 1.0, 50, 0.82)
    np.testing.assert_allclose(calibrated_d, perf_img.data / m0_expected)

def test_refregion_gm_all():
    """
    Reference region mask is whole image
    """
    perf_img, calib_img = _get_imgs()

    ref_d = np.ones((5, 5, 5))
    ref_img = Image(name="refmask", image=ref_d)

    wsp = Workspace(calib=calib_img, calib_method="refregion", refmask=ref_img, tissref="gm")
    perf_calib = calib.calibrate(wsp, perf_img)
    calibrated_d = perf_calib.data

    # GM defaults
    m0_expected =  _expected_m0(np.mean(calib_img.data), 1.3, 100, 0.98)
    np.testing.assert_allclose(calibrated_d, perf_img.data / m0_expected)

def test_refregion_tr():
    """
    Reference region override TR
    """
    TR = 7.1
    perf_img, calib_img = _get_imgs()

    ref_d = np.ones((5, 5, 5))
    ref_img = Image(name="refmask", image=ref_d)

    wsp = Workspace(calib=calib_img, calib_method="refregion", refmask=ref_img, tissref="wm", tr=TR)
    perf_calib = calib.calibrate(wsp, perf_img)
    calibrated_d = perf_calib.data

    # WM defaults
    m0_expected =  _expected_m0(np.mean(calib_img.data), 1.0, 50, 0.82, tr=TR)
    np.testing.assert_allclose(calibrated_d, perf_img.data / m0_expected)

def test_refregion_te():
    """
    Reference region override TE
    """
    TE = 1.2
    perf_img, calib_img = _get_imgs()

    ref_d = np.ones((5, 5, 5))
    ref_img = Image(name="refmask", image=ref_d)

    wsp = Workspace(calib=calib_img, calib_method="refregion", refmask=ref_img, tissref="wm", te=TE)
    perf_calib = calib.calibrate(wsp, perf_img)
    calibrated_d = perf_calib.data

    m0_expected =  _expected_m0(np.mean(calib_img.data), 1.0, 50, 0.82, te=TE)
    np.testing.assert_allclose(calibrated_d, perf_img.data / m0_expected)

def test_refregion_t2b():
    """
    Reference region override T2b
    """
    T2b = 782
    perf_img, calib_img = _get_imgs()

    ref_d = np.ones((5, 5, 5))
    ref_img = Image(name="refmask", image=ref_d)

    wsp = Workspace(calib=calib_img, calib_method="refregion", refmask=ref_img, tissref="wm", t2b=T2b)
    perf_calib = calib.calibrate(wsp, perf_img)
    calibrated_d = perf_calib.data

    m0_expected =  _expected_m0(np.mean(calib_img.data), 1.0, 50, 0.82, t2b=T2b)
    np.testing.assert_allclose(calibrated_d, perf_img.data / m0_expected)

def test_refregion_t1r():
    """
    Reference region override T1 of ref tissue
    """
    T1 = 1.03
    perf_img, calib_img = _get_imgs()

    ref_d = np.ones((5, 5, 5))
    ref_img = Image(name="refmask", image=ref_d)

    wsp = Workspace(calib=calib_img, calib_method="refregion", refmask=ref_img, tissref="wm", t1r=T1)
    perf_calib = calib.calibrate(wsp, perf_img)
    calibrated_d = perf_calib.data

    m0_expected =  _expected_m0(np.mean(calib_img.data), T1, 50, 0.82)
    np.testing.assert_allclose(calibrated_d, perf_img.data / m0_expected)

def test_refregion_t2r():
    """
    Reference region override T2 of ref tissue
    """
    T2 = 635
    perf_img, calib_img = _get_imgs()

    ref_d = np.ones((5, 5, 5))
    ref_img = Image(name="refmask", image=ref_d)

    wsp = Workspace(calib=calib_img, calib_method="refregion", refmask=ref_img, tissref="wm", t2r=T2)
    perf_calib = calib.calibrate(wsp, perf_img)
    calibrated_d = perf_calib.data

    m0_expected =  _expected_m0(np.mean(calib_img.data), 1.0, T2, 0.82)
    np.testing.assert_allclose(calibrated_d, perf_img.data / m0_expected)

def test_refregion_pc():
    """
    Reference region override partition coeff of ref tissue
    """
    PC = 0.67
    perf_img, calib_img = _get_imgs()

    ref_d = np.ones((5, 5, 5))
    ref_img = Image(name="refmask", image=ref_d)

    wsp = Workspace(calib=calib_img, calib_method="refregion", refmask=ref_img, tissref="wm", pcr=PC)
    perf_calib = calib.calibrate(wsp, perf_img)
    calibrated_d = perf_calib.data

    m0_expected =  _expected_m0(np.mean(calib_img.data), 1.0, 50, PC)
    np.testing.assert_allclose(calibrated_d, perf_img.data / m0_expected)

def test_refregion_csf_t2star():
    """
    Using T2* instead of T2 for CSF
    """
    perf_img, calib_img = _get_imgs()
    
    ref_d = np.ones((5, 5, 5))
    ref_img = Image(name="refmask", image=ref_d)

    wsp = Workspace(calib=calib_img, calib_method="refregion", refmask=ref_img, tissref="csf", t2star=True)
    perf_calib = calib.calibrate(wsp, perf_img)
    calibrated_d = perf_calib.data

    m0_expected =  _expected_m0(np.mean(calib_img.data), 4.3, 400, 1.15)
    np.testing.assert_allclose(calibrated_d, perf_img.data / m0_expected)

def test_refregion_wm_t2star():
    """
    Using T2* instead of T2 for WM
    """
    perf_img, calib_img = _get_imgs()

    ref_d = np.ones((5, 5, 5))
    ref_img = Image(name="refmask", image=ref_d)

    wsp = Workspace(calib=calib_img, calib_method="refregion", refmask=ref_img, tissref="wm", t2star=True)
    perf_calib = calib.calibrate(wsp, perf_img)
    calibrated_d = perf_calib.data

    # WM defaults
    m0_expected =  _expected_m0(np.mean(calib_img.data), 1.0, 50, 0.82)
    np.testing.assert_allclose(calibrated_d, perf_img.data / m0_expected)

def test_refregion_gm_t2star():
    """
    Using T2* instead of T2 for GM
    """
    perf_img, calib_img = _get_imgs()

    ref_d = np.ones((5, 5, 5))
    ref_img = Image(name="refmask", image=ref_d)

    wsp = Workspace(calib=calib_img, calib_method="refregion", refmask=ref_img, tissref="gm", t2star=True)
    perf_calib = calib.calibrate(wsp, perf_img)
    calibrated_d = perf_calib.data

    # GM defaults
    m0_expected =  _expected_m0(np.mean(calib_img.data), 1.3, 60, 0.98)
    np.testing.assert_allclose(calibrated_d, perf_img.data / m0_expected)

def test_refregion_gain():
    """
    Modify calibration gain
    """
    GAIN = 3.4
    perf_img, calib_img = _get_imgs()

    ref_d = np.ones((5, 5, 5))
    ref_img = Image(name="refmask", image=ref_d)

    wsp = Workspace(calib=calib_img, calib_method="refregion", refmask=ref_img, tissref="wm", calib_gain=GAIN)
    perf_calib = calib.calibrate(wsp, perf_img)
    calibrated_d = perf_calib.data

    m0_expected =  _expected_m0(np.mean(calib_img.data), 1.0, 50, 0.82, gain=GAIN)
    np.testing.assert_allclose(calibrated_d, perf_img.data / m0_expected)

def test_refregion_alpha():
    """
    Modify calibration gain
    """
    ALPHA = 0.46
    perf_img, calib_img = _get_imgs()

    ref_d = np.ones((5, 5, 5))
    ref_img = Image(name="refmask", image=ref_d)

    wsp = Workspace(calib=calib_img, calib_method="refregion", refmask=ref_img, tissref="wm")
    perf_calib = calib.calibrate(wsp, perf_img, alpha=ALPHA)
    calibrated_d = perf_calib.data

    m0_expected =  _expected_m0(np.mean(calib_img.data), 1.0, 50, 0.82, alpha=ALPHA)
    np.testing.assert_allclose(calibrated_d, perf_img.data / m0_expected)
