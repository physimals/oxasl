"""
Tests for BASIL module

FIXME need more multi-step tests
FIXME need to test more error conditions
"""
import StringIO

import pytest
import numpy as np

from fsl.data.image import Image

from oxasl import basil, AslImage

DEFAULTS = {
    "method" : "vb",
    "noise" : "white",
    "allow-bad-voxels" : True,
    "convergence" : "trialmode",
    "max-iterations" : 20,
    "max-trials" : 10,
    "model" : "aslrest",
    "disp" : "none",
    "exch" : "mix",
    "inctiss" : True,
    "incbat" : True,
    "infertiss" : True,
    "inferbat" : True,
}

def _get_defaults(img):
    options = dict(DEFAULTS)
    for idx, ti in enumerate(img.tis):
        options["ti%i" % (idx+1)] = ti
        options["rpt%i" % (idx+1)] = img.rpts[idx]
    return options


def _check_step(step, step_num=None, desc_text=None, data_name=None, options=None, prev_step=None):
    if step_num:
        assert(step[0] == step_num)
    if desc_text:
        assert(desc_text.lower().strip() in step[1].lower())
    if data_name:
        assert(step[2]["data"].name == data_name)
    if options:
        for k, v in options.items():
            assert(step[2][k] == v)
    if prev_step:
        assert(step[3] == prev_step)

def test_nodata():
    """
    Check we get an error if there is nothing to infer
    """
    log = StringIO.StringIO()

    with pytest.raises(ValueError):
        steps = basil.get_steps(None, log=log)

def test_infer_nothing():
    """
    Check we get an error if there is nothing to infer
    """
    d = np.random.rand(5, 5, 5, 6)
    img = AslImage(name="asldata", image=d, tis=[1.5], order="prt")
    log = StringIO.StringIO()

    with pytest.raises(ValueError):
        steps = basil.get_steps(img, infertiss=False, inferbat=False, log=log)

def test_defaults():
    """
    Check the basic defaults (infer tissue perfusion and bolus arrival time)
    """
    d = np.random.rand(5, 5, 5, 6)
    img = AslImage(name="asldata", image=d, tis=[1.5], order="prt")
    log = StringIO.StringIO()

    options = _get_defaults(img)
    steps = basil.get_steps(img, log=log)
    assert(len(steps) == 1)

    _check_step(steps[0], step_num=1, desc_text="tissue", 
                options=options)
                
def test_fix_bat():
    """
    Check fixing the arrival time, which is normally inferred
    """
    d = np.random.rand(5, 5, 5, 6)
    img = AslImage(name="asldata", image=d, tis=[1.5], order="prt")
    log = StringIO.StringIO()

    options = _get_defaults(img)
    steps = basil.get_steps(img, infer_bat=False, log=log)
    assert(len(steps) == 1)

    options.pop("incbat")
    options.pop("inferbat")
    _check_step(steps[0], step_num=1, desc_text="tissue", 
                options=options)
 
def test_inferart():
    """
    Check inference of arterial component
    """
    d = np.random.rand(5, 5, 5, 6)
    img = AslImage(name="asldata", image=d, tis=[1.5], order="prt")
    log = StringIO.StringIO()

    options = _get_defaults(img)
    steps = basil.get_steps(img, inferart=True, log=log)
    assert(len(steps) == 2)

    _check_step(steps[0], step_num=1, desc_text="tissue", 
                options=dict(options, **{
                    "incart" : True,
                }))

    _check_step(steps[1], step_num=2, desc_text="arterial", 
                options=dict(options, **{
                    "incart" : True,
                    "inferart" : True,
                }))

def test_infertau():
    """
    Check inference of bolus duration (tau)
    """
    d = np.random.rand(5, 5, 5, 6)
    img = AslImage(name="asldata", image=d, tis=[1.5], order="prt")
    log = StringIO.StringIO()

    options = _get_defaults(img)
    steps = basil.get_steps(img, infertau=True, log=log)
    assert(len(steps) == 2)

    _check_step(steps[0], step_num=1, desc_text="tissue", 
                options=dict(options, **{
                    "inctau" : True,
                }))

    _check_step(steps[1], step_num=2, desc_text="bolus", 
                options=dict(options, **{
                    "inctau" : True,
                    "infertau" : True,
                }))

def test_inferarttau():
    """
    Check inference of bolus duration (tau) and arterial component
    """
    d = np.random.rand(5, 5, 5, 6)
    img = AslImage(name="asldata", image=d, tis=[1.5], order="prt")
    log = StringIO.StringIO()

    options = _get_defaults(img)
    steps = basil.get_steps(img, infertau=True, inferart=True, log=log)
    assert(len(steps) == 3)

    _check_step(steps[0], step_num=1, desc_text="tissue", 
                options=dict(options, **{
                    "inctau" : True,
                    "incart" : True,
                }))

    _check_step(steps[1], step_num=2, desc_text="arterial", 
                options=dict(options, **{
                    "inctau" : True,
                    "incart" : True,
                    "inferart" : True,
                }))

    _check_step(steps[2], step_num=3, desc_text="bolus", 
                options=dict(options, **{
                    "inctau" : True,
                    "incart" : True,
                    "infertau" : True,
                    "inferart" : True,
                }))

def test_infert1():
    """
    Check inference of T1
    """
    d = np.random.rand(5, 5, 5, 6)
    img = AslImage(name="asldata", image=d, tis=[1.5], order="prt")
    log = StringIO.StringIO()

    options = _get_defaults(img)
    steps = basil.get_steps(img, infert1=True, log=log)
    assert(len(steps) == 2)

    _check_step(steps[0], step_num=1, desc_text="tissue", 
                options=dict(options, **{
                    "inct1" : True,
                }))

    _check_step(steps[1], step_num=2, desc_text="T1", 
                options=dict(options, **{
                    "inct1" : True,
                    "infert1" : True,
                }))

def test_t1im():
    """
    Check T1 image priors are correctly handled
    """
    d = np.random.rand(5, 5, 5, 6)
    img = AslImage(name="asldata", image=d, tis=[1.5], order="prt")
    log = StringIO.StringIO()

    t1d = np.random.rand(5, 5, 5)
    t1im = Image(name="t1file", image=t1d)
    steps = basil.get_steps(img, infert1=True, t1im=t1im, log=log)
    assert(len(steps) == 2)

    options = _get_defaults(img)
    _check_step(steps[0], step_num=1, desc_text="tissue", 
                options=dict(options, **{
                    "inct1" : True,
                }))

    _check_step(steps[1], step_num=2, desc_text="T1", 
                options=dict(options, **{
                    "inct1" : True,
                    "infert1" : True,
                    "PSP_byname1" : "T_1",
                    "PSP_byname1_type" : "I",
                    #"PSP_byname1_image" : "t1file",
                }))

def test_inferpc():
    """
    Check the pre-capiliary component
    """
    d = np.random.rand(5, 5, 5, 6)
    img = AslImage(name="asldata", image=d, tis=[1.5], order="prt")
    log = StringIO.StringIO()

    options = _get_defaults(img)
    steps = basil.get_steps(img, inferpc=True, log=log)
    assert(len(steps) == 2)

    _check_step(steps[0], step_num=1, desc_text="tissue", 
                options=dict(options, **{
                    "incpc" : True,
                }))

    _check_step(steps[1], step_num=2, desc_text="pre-capiliary", 
                options=dict(options, **{
                    "incpc" : True,
                    "inferpc" : True,
                }))

def test_artonly():
    """
    Check we can infer arterial component without tissue step
    """
    d = np.random.rand(5, 5, 5, 6)
    img = AslImage(name="asldata", image=d, tis=[1.5], order="prt")
    log = StringIO.StringIO()

    steps = basil.get_steps(img, infertiss=False, inferart=True, log=log)
    assert(len(steps) == 1)

    options = _get_defaults(img)
    options.update({
        "incart" : True,
        "inferart" : True,
    })
    options.pop("inctiss")
    options.pop("infertiss")
    _check_step(steps[0], step_num=1, desc_text="arterial", 
                options=options)

def test_initmvn():
    """
    Check the supply of an initialization MVN
    """
    d = np.random.rand(5, 5, 5, 6)
    img = AslImage(name="asldata", image=d, tis=[1.5], order="prt")
    log = StringIO.StringIO()

    mvnd = np.random.rand(5, 5, 5, 6)
    initmvn = Image(name="mvnfile", image=mvnd)
    steps = basil.get_steps(img, initmvn=initmvn, log=log)
    assert(len(steps) == 1)

    options = _get_defaults(img)
    options.update({
        "continue-from-mvn" : "mvnfile"
    })
    _check_step(steps[0], step_num=1, desc_text="tissue", 
                options=options)

def test_spatial():
    """
    Check final spatial step
    """
    d = np.random.rand(5, 5, 5, 6)
    img = AslImage(name="asldata", image=d, tis=[1.5], order="prt")
    log = StringIO.StringIO()

    options = _get_defaults(img)
    steps = basil.get_steps(img, spatial=True, log=log)
    assert(len(steps) == 2)

    _check_step(steps[0], step_num=1, desc_text="tissue", 
                options=options)

    options.update({
        "method" : "spatialvb",
        "param-spatial-priors" : "N+",
        "PSP_byname1" : "ftiss",
        "PSP_byname1_type" : "M",
#        "convergence" : "maxiters", FIXME
    })
    options.pop("max-trials")
    _check_step(steps[1], step_num=2, desc_text="spatial", options=options)
    
def test_onestep():
    """
    Check that single step mode works when you would normally get multiple steps
    """
    d = np.random.rand(5, 5, 5, 6)
    img = AslImage(name="asldata", image=d, tis=[1.5], order="prt")
    log = StringIO.StringIO()

    steps = basil.get_steps(img, infertau=True, inferart=True, spatial=True, onestep=True, log=log)
    assert(len(steps) == 1)

    options = _get_defaults(img)
    options.update({
        "method" : "spatialvb",
        "param-spatial-priors" : "N+",
        "PSP_byname1" : "ftiss",
        "PSP_byname1_type" : "M",
        "inctau" : True,
        "incart" : True,
        "inferart" : True,
        "infertau" : True,
#        "convergence" : "maxiters", FIXME
    })
    options.pop("max-trials")
    _check_step(steps[0], step_num=1, desc_text="spatial", options=options)

def test_max_iterations():
    """
    Check that max iterations can be overridden
    """
    d = np.random.rand(5, 5, 5, 6)
    img = AslImage(name="asldata", image=d, tis=[1.5], order="prt")
    log = StringIO.StringIO()

    kwargs = {
        "max-iterations" : 123,
    }
    steps = basil.get_steps(img, log=log, **kwargs)
    assert(len(steps) == 1)

    options = _get_defaults(img)
    options.update({
        "max-iterations" : 123,
    })
    _check_step(steps[0], step_num=1, desc_text="tissue", options=options)

def test_random_extra_options():
    """
    Check that any additional keyword arguments are passed to Fabber
    """
    d = np.random.rand(5, 5, 5, 6)
    img = AslImage(name="asldata", image=d, tis=[1.5], order="prt")
    log = StringIO.StringIO()

    kwargs = {
        "phase-of-moon-correction-factor" : 7,
        "random-output-proportion-percent" : 36,
    }
    steps = basil.get_steps(img, log=log, **kwargs)
    assert(len(steps) == 1)

    options = _get_defaults(img)
    options.update(kwargs)
    _check_step(steps[0], step_num=1, desc_text="tissue", options=options)

def test_pvc_only_one_map_given1():
    """
    Check that PVC correction fails if you only give the GM map
    """
    d = np.random.rand(5, 5, 5, 6)
    img = AslImage(name="asldata", image=d, tis=[1.5], order="prt")
    log = StringIO.StringIO()

    pgmd = np.random.rand(5, 5, 5)
    pgm = Image(name="pgm_map", image=pgmd)
    with pytest.raises(ValueError):
        basil.get_steps(img, pgm=pgm, log=log)
    
def test_pvc_only_one_map_given2():
    """
    Check that PVC correction fails if you only give the WM map
    """
    d = np.random.rand(5, 5, 5, 6)
    img = AslImage(name="asldata", image=d, tis=[1.5], order="prt")
    log = StringIO.StringIO()

    pwmd = np.random.rand(5, 5, 5)
    pwm = Image(name="pwm_map", image=pwmd)
    with pytest.raises(ValueError):
        basil.get_steps(img, pwm=pwm, log=log)
    
def test_pvc_no_tissue():
    """
    Check that PVC correction fails if you do not infer the tissue component
    """
    d = np.random.rand(5, 5, 5, 6)
    img = AslImage(name="asldata", image=d, tis=[1.5], order="prt")
    log = StringIO.StringIO()

    pgmd = np.random.rand(5, 5, 5)
    pgm = Image(name="pgm_map", image=pgmd)

    pwmd = np.random.rand(5, 5, 5)
    pwm = Image(name="pwm_map", image=pwmd)

    with pytest.raises(ValueError):
        basil.get_steps(img, pgm=pgm, pwm=pwm, infertiss=False, log=log)
    
def test_pvc():
    """
    FIXME we need to test the PVC initialization step
    and how to do this is not finalized
    """
    d = np.random.rand(5, 5, 5, 6)
    img = AslImage(name="asldata", image=d, tis=[1.5], order="prt")
    log = StringIO.StringIO()

    pgmd = np.random.rand(5, 5, 5)
    pgm = Image(name="pgm_map", image=pgmd)

    pwmd = np.random.rand(5, 5, 5)
    pwm = Image(name="pwm_map", image=pwmd)

    steps = basil.get_steps(img, pgm=pgm, pwm=pwm, log=log)
    assert(len(steps) == 3)

    options = _get_defaults(img)
    options.update({
        "incpve" : True,
    })
    _check_step(steps[0], step_num=1, desc_text="tissue", options=options)

    options.update({
        "method" : "spatialvb",
        "param-spatial-priors" : "N+",
        "PSP_byname1" : "ftiss",
        "PSP_byname1_type" : "M",
        "max-iterations" : 200,
#        "convergence" : "maxiters", FIXME
    })
    options.pop("max-trials")
    _check_step(steps[2], step_num=3, desc_text="PVE", options=options)
    _check_step(steps[2], step_num=3, desc_text="spatial")

