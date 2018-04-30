"""
Tests for CALIB module


calib(perf_data, calib_data, method, output_name=None, multiplier=1.0, var=False, log=sys.stdout, **kwargs)
get_m0_voxelwise(calib_data, gain=1.0, alpha=1.0, tr=None, t1t=None, part_coeff=0.9, mask=None, edgecorr=False, log=sys.stdout)
"""
import math
import StringIO

import pytest
import numpy as np

from oxasl import calib, AslImage, fsl

def test_nodata():
    """
    Check we get an error if there is no perfusion data
    """
    d = np.random.rand(5, 5, 5, 6)
    calib_img = fsl.Image("calib", data=d)

    log = StringIO.StringIO()
    with pytest.raises(ValueError):
        calib(None, calib_img, "voxelwise", log=log)

def test_nocalib():
    """
    Check we get an error if there is no calibration data
    """
    d = np.random.rand(5, 5, 5, 6)
    perf_img = fsl.Image("perfusion", data=d)

    log = StringIO.StringIO()
    with pytest.raises(ValueError):
        calib(perf_img, None, "voxelwise", log=log)

def test_bad_method():
    """
    Check we get an error on an invalid method
    """
    d = np.random.rand(5, 5, 5, 6)
    perf_img = fsl.Image("perfusion", data=d)

    d = np.random.rand(5, 5, 5, 6)
    calib_img = fsl.Image("calib", data=d)

    log = StringIO.StringIO()
    with pytest.raises(ValueError):
        calib(perf_img, calib_img, "random", log=log)

def test_defaults():
    """
    Check default calibration works
    """
    perf_d = np.random.rand(5, 5, 5)
    perf_img = fsl.Image("perfusion", data=perf_d)

    calib_d = np.random.rand(5, 5, 5)
    calib_img = fsl.Image("calib", data=calib_d)

    log = StringIO.StringIO()
    perf_calib = calib(perf_img, calib_img, "voxelwise", log=log)
    calibrated_d = perf_calib.data()
    assert(perf_calib.iname == "perfusion_calib")
    assert(perf_calib.shape == perf_img.shape)

    # Default partition coefficient is 0.9
    np.testing.assert_allclose(calibrated_d, 0.9 * perf_d / calib_d)

def test_cgain():
    """
    Check calibration gain
    """
    GAIN = 1.23
    perf_d = np.random.rand(5, 5, 5)
    perf_img = fsl.Image("perfusion", data=perf_d)

    calib_d = np.random.rand(5, 5, 5)
    calib_img = fsl.Image("calib", data=calib_d)

    log = StringIO.StringIO()
    perf_calib = calib(perf_img, calib_img, "voxelwise", gain=GAIN, log=log)
    calibrated_d = perf_calib.data()
    assert(perf_calib.iname == "perfusion_calib")
    assert(perf_calib.shape == perf_img.shape)

    # Default partition coefficient is 0.9
    np.testing.assert_allclose(calibrated_d, 0.9 / GAIN * perf_d / calib_d)

def test_alpha():
    """
    Check inversion efficiency
    """
    ALPHA = 0.74
    perf_d = np.random.rand(5, 5, 5)
    perf_img = fsl.Image("perfusion", data=perf_d)

    calib_d = np.random.rand(5, 5, 5)
    calib_img = fsl.Image("calib", data=calib_d)

    log = StringIO.StringIO()
    perf_calib = calib(perf_img, calib_img, "voxelwise", alpha=ALPHA, log=log)
    calibrated_d = perf_calib.data()
    assert(perf_calib.iname == "perfusion_calib")
    assert(perf_calib.shape == perf_img.shape)

    # Default partition coefficient is 0.9
    np.testing.assert_allclose(calibrated_d, 0.9 / ALPHA * perf_d / calib_d)

def test_multiplier():
    """
    Check end multiplier for physical units
    """
    MULTIPLIER = 7654
    perf_d = np.random.rand(5, 5, 5)
    perf_img = fsl.Image("perfusion", data=perf_d)

    calib_d = np.random.rand(5, 5, 5)
    calib_img = fsl.Image("calib", data=calib_d)

    log = StringIO.StringIO()
    perf_calib = calib(perf_img, calib_img, "voxelwise", multiplier=MULTIPLIER, log=log)
    calibrated_d = perf_calib.data()
    assert(perf_calib.iname == "perfusion_calib")
    assert(perf_calib.shape == perf_img.shape)

    # Default partition coefficient is 0.9
    np.testing.assert_allclose(calibrated_d, 0.9 * MULTIPLIER * perf_d / calib_d)

def test_partition_coeff():
    """
    Check tissue/blood partition coefficient
    """
    PC = 0.67
    perf_d = np.random.rand(5, 5, 5)
    perf_img = fsl.Image("perfusion", data=perf_d)

    calib_d = np.random.rand(5, 5, 5)
    calib_img = fsl.Image("calib", data=calib_d)

    log = StringIO.StringIO()
    perf_calib = calib(perf_img, calib_img, "voxelwise", part_coeff=PC, log=log)
    calibrated_d = perf_calib.data()
    assert(perf_calib.iname == "perfusion_calib")
    assert(perf_calib.shape == perf_img.shape)

    # Default partition coefficient is 0.9
    np.testing.assert_allclose(calibrated_d, PC * perf_d / calib_d)

def test_shorttr_corr():
    """
    Check correction for short TR
    """
    TR, T1 = 3.0, 1.1
    perf_d = np.random.rand(5, 5, 5)
    perf_img = fsl.Image("perfusion", data=perf_d)

    calib_d = np.random.rand(5, 5, 5)
    calib_img = fsl.Image("calib", data=calib_d)

    log = StringIO.StringIO()
    perf_calib = calib(perf_img, calib_img, "voxelwise", tr=TR, t1t=T1, log=log)
    calibrated_d = perf_calib.data()
    assert(perf_calib.iname == "perfusion_calib")
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
    perf_img = fsl.Image("perfusion", data=perf_d)

    calib_d = np.random.rand(5, 5, 5)
    calib_img = fsl.Image("calib", data=calib_d)

    log = StringIO.StringIO()
    perf_calib = calib(perf_img, calib_img, "voxelwise", tr=TR, log=log)
    calibrated_d = perf_calib.data()
    assert(perf_calib.iname == "perfusion_calib")
    assert(perf_calib.shape == perf_img.shape)
    assert("WARNING" in log.getvalue())

    # Default partition coefficient is 0.9
    np.testing.assert_allclose(calibrated_d, 0.9 * perf_d / calib_d)

def test_defaults_var():
    """
    Check default calibration works on variances
    """
    perf_d = np.random.rand(5, 5, 5)
    perf_img = fsl.Image("perfusion", data=perf_d)

    calib_d = np.random.rand(5, 5, 5)
    calib_img = fsl.Image("calib", data=calib_d)

    log = StringIO.StringIO()
    perf_calib = calib(perf_img, calib_img, "voxelwise", var=True, log=log)
    calibrated_d = perf_calib.data()
    assert(perf_calib.iname == "perfusion_calib")
    assert(perf_calib.shape == perf_img.shape)

    # Default partition coefficient is 0.9
    np.testing.assert_allclose(calibrated_d, 0.9 * 0.9 * perf_d / calib_d / calib_d)

def test_cgain_var():
    """
    Check calibration gain on variances
    """
    GAIN = 1.23
    perf_d = np.random.rand(5, 5, 5)
    perf_img = fsl.Image("perfusion", data=perf_d)

    calib_d = np.random.rand(5, 5, 5)
    calib_img = fsl.Image("calib", data=calib_d)

    log = StringIO.StringIO()
    perf_calib = calib(perf_img, calib_img, "voxelwise", gain=GAIN, var=True, log=log)
    calibrated_d = perf_calib.data()
    assert(perf_calib.iname == "perfusion_calib")
    assert(perf_calib.shape == perf_img.shape)

    # Default partition coefficient is 0.9
    np.testing.assert_allclose(calibrated_d, 0.9 * 0.9 / GAIN / GAIN * perf_d / calib_d / calib_d)

def test_alpha_var():
    """
    Check inversion efficiency on variances
    """
    ALPHA = 0.74
    perf_d = np.random.rand(5, 5, 5)
    perf_img = fsl.Image("perfusion", data=perf_d)

    calib_d = np.random.rand(5, 5, 5)
    calib_img = fsl.Image("calib", data=calib_d)

    log = StringIO.StringIO()
    perf_calib = calib(perf_img, calib_img, "voxelwise", alpha=ALPHA, var=True, log=log)
    calibrated_d = perf_calib.data()
    assert(perf_calib.iname == "perfusion_calib")
    assert(perf_calib.shape == perf_img.shape)

    # Default partition coefficient is 0.9
    np.testing.assert_allclose(calibrated_d, 0.9 * 0.9 / ALPHA /ALPHA * perf_d / calib_d / calib_d)

def test_multiplier_var():
    """
    Check end multiplier for physical units on variances
    """
    MULTIPLIER = 7654
    perf_d = np.random.rand(5, 5, 5)
    perf_img = fsl.Image("perfusion", data=perf_d)

    calib_d = np.random.rand(5, 5, 5)
    calib_img = fsl.Image("calib", data=calib_d)

    log = StringIO.StringIO()
    perf_calib = calib(perf_img, calib_img, "voxelwise", multiplier=MULTIPLIER, var=True, log=log)
    calibrated_d = perf_calib.data()
    assert(perf_calib.iname == "perfusion_calib")
    assert(perf_calib.shape == perf_img.shape)

    # Default partition coefficient is 0.9
    np.testing.assert_allclose(calibrated_d, 0.9 * 0.9 * MULTIPLIER * MULTIPLIER * perf_d / calib_d / calib_d)

def test_partition_coeff_var():
    """
    Check tissue/blood partition coefficient on variances
    """
    PC = 0.67
    perf_d = np.random.rand(5, 5, 5)
    perf_img = fsl.Image("perfusion", data=perf_d)

    calib_d = np.random.rand(5, 5, 5)
    calib_img = fsl.Image("calib", data=calib_d)

    log = StringIO.StringIO()
    perf_calib = calib(perf_img, calib_img, "voxelwise", part_coeff=PC, var=True, log=log)
    calibrated_d = perf_calib.data()
    assert(perf_calib.iname == "perfusion_calib")
    assert(perf_calib.shape == perf_img.shape)

    # Default partition coefficient is 0.9
    np.testing.assert_allclose(calibrated_d, PC * PC * perf_d / calib_d / calib_d)

def test_shorttr_corr_var():
    """
    Check correction for short TR on variances
    """
    TR, T1 = 3.0, 1.1
    perf_d = np.random.rand(5, 5, 5)
    perf_img = fsl.Image("perfusion", data=perf_d)

    calib_d = np.random.rand(5, 5, 5)
    calib_img = fsl.Image("calib", data=calib_d)

    log = StringIO.StringIO()
    perf_calib = calib(perf_img, calib_img, "voxelwise", tr=TR, t1t=T1, var=True, log=log)
    calibrated_d = perf_calib.data()
    assert(perf_calib.iname == "perfusion_calib")
    assert(perf_calib.shape == perf_img.shape)

    factor =  1 - math.exp(-TR / T1)
    # Default partition coefficient is 0.9
    np.testing.assert_allclose(calibrated_d, 0.9 * 0.9 * factor * factor * perf_d / calib_d / calib_d)

def test_refregion_not_implemented():
    TR, T1 = 3.0, 1.1
    perf_d = np.random.rand(5, 5, 5)
    perf_img = fsl.Image("perfusion", data=perf_d)

    calib_d = np.random.rand(5, 5, 5)
    calib_img = fsl.Image("calib", data=calib_d)

    ref_d = np.ones((5, 5, 5))
    ref_img = fsl.Image("ref_mask", data=ref_d)

    log = StringIO.StringIO()
    with pytest.raises(NotImplementedError):
        calib(perf_img, calib_img, "refregion", ref_mask=ref_img, log=log)