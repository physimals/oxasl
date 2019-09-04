"""
Tests for BASIL module

FIXME need more multi-step tests
FIXME need to test more error conditions
"""
import pytest
import numpy as np

from fsl.data.image import Image

from oxasl import AslImage, Workspace
import oxasl.basil as basil

pytest.skip("skipping basil tests as neeed rewrite", allow_module_level=True)

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

def _check_step(step, desc_text=None, options=None):
    print(step.options)
    if desc_text:
        assert(desc_text.lower().strip() in step.desc.lower())
    if options:
        for k, v in options.items():
            assert(step.options[k] == v)

def test_nodata():
    """
    Check we get an error if there is no data
    """
    wsp = Workspace()

    with pytest.raises(ValueError):
        steps = basil.basil_steps(wsp, None)

def test_infer_nothing():
    """
    Check we get an error if there is nothing to infer
    """
    d = np.random.rand(5, 5, 5, 6)
    img = AslImage(name="asldata", image=d, tis=[1.5], order="prt")
    wsp = Workspace()

    with pytest.raises(ValueError):
        steps = basil.basil_steps(wsp, img)

def test_defaults():
    """
    Check the basic defaults (infer tissue perfusion and bolus arrival time)
    """
    d = np.random.rand(5, 5, 5, 6)
    img = AslImage(name="asldata", image=d, tis=[1.5], order="prt")
    wsp = Workspace()
    wsp.infertiss = True
    wsp.inferbat = True

    steps = basil.basil_steps(wsp, img)
    assert(len(steps) == 1)

    options = _get_defaults(img)
    _check_step(steps[0], desc_text="tissue", options=options)
                
def test_fix_bat():
    """
    Check fixing the arrival time, which is normally inferred
    """
    d = np.random.rand(5, 5, 5, 6)
    img = AslImage(name="asldata", image=d, tis=[1.5], order="prt")
    wsp = Workspace(infertiss=True, inferbat=False)

    steps = basil.basil_steps(wsp, img)
    assert(len(steps) == 1)

    options = _get_defaults(img)
    options.pop("incbat")
    options.pop("inferbat")
    _check_step(steps[0], desc_text="tissue", options=options)
 
def test_inferart():
    """
    Check inference of arterial component
    """
    d = np.random.rand(5, 5, 5, 6)
    img = AslImage(name="asldata", image=d, tis=[1.5], order="prt")
    wsp = Workspace(infertiss=True, inferbat=True, inferart=True)

    options = _get_defaults(img)
    steps = basil.basil_steps(wsp, img)
    assert(len(steps) == 2)

    _check_step(steps[0], desc_text="tissue", 
                options=dict(options, **{
                    "incart" : True,
                }))

    _check_step(steps[1], desc_text="arterial", 
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
    wsp = Workspace(infertiss=True, inferbat=True, infertau=True)

    steps = basil.basil_steps(wsp, img)
    assert(len(steps) == 2)

    options = _get_defaults(img)
    _check_step(steps[0], desc_text="tissue", 
                options=dict(options, **{
                    "inctau" : True,
                }))

    _check_step(steps[1], desc_text="bolus", 
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
    wsp = Workspace(infertiss=True, inferbat=True, infertau=True, inferart=True)

    steps = basil.basil_steps(wsp, img)
    assert(len(steps) == 3)

    options = _get_defaults(img)
    _check_step(steps[0], desc_text="tissue", 
                options=dict(options, **{
                    "inctau" : True,
                    "incart" : True,
                }))

    _check_step(steps[1], desc_text="arterial", 
                options=dict(options, **{
                    "inctau" : True,
                    "incart" : True,
                    "inferart" : True,
                }))

    _check_step(steps[2], desc_text="bolus", 
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
    wsp = Workspace(infertiss=True, inferbat=True, infert1=True)

    steps = basil.basil_steps(wsp, img)
    assert(len(steps) == 2)

    options = _get_defaults(img)
    _check_step(steps[0], desc_text="tissue", 
                options=dict(options, **{
                    "inct1" : True,
                }))

    _check_step(steps[1], desc_text="T1", 
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
    wsp = Workspace(infertiss=True, inferbat=True, infert1=True)

    wsp.t1im = Image(np.random.rand(5, 5, 5))
    steps = basil.basil_steps(wsp, img)
    assert(len(steps) == 2)

    options = _get_defaults(img)
    _check_step(steps[0], desc_text="tissue", 
                options=dict(options, **{
                    "inct1" : True,
                }))

    _check_step(steps[1], desc_text="T1", 
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
    wsp = Workspace(infertiss=True, inferbat=True, inferpc=True)

    steps = basil.basil_steps(wsp, img)
    assert(len(steps) == 2)

    options = _get_defaults(img)
    _check_step(steps[0], desc_text="tissue", 
                options=dict(options, **{
                    "incpc" : True,
                }))

    _check_step(steps[1], desc_text="pre-capiliary", 
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
    wsp = Workspace(infertiss=False, inferbat=True, inferart=True)

    steps = basil.basil_steps(wsp, img)
    assert(len(steps) == 1)

    options = _get_defaults(img)
    options.update({
        "incart" : True,
        "inferart" : True,
    })
    options.pop("inctiss")
    options.pop("infertiss")
    _check_step(steps[0], desc_text="arterial", options=options)

def test_initmvn():
    """
    Check the supply of an initialization MVN
    """
    d = np.random.rand(5, 5, 5, 6)
    img = AslImage(name="asldata", image=d, tis=[1.5], order="prt")
    wsp = Workspace(infertiss=True, inferbat=True)

    wsp.initmvn = Image(np.random.rand(5, 5, 5, 6))
    steps = basil.basil_steps(wsp, img)
    assert(len(steps) == 1)

    options = _get_defaults(img)
    options.update({
        "continue-from-mvn" : wsp.initmvn
    })
    _check_step(steps[0], desc_text="tissue", options=options)

def test_spatial():
    """
    Check final spatial step
    """
    d = np.random.rand(5, 5, 5, 6)
    img = AslImage(name="asldata", image=d, tis=[1.5], order="prt")
    wsp = Workspace(infertiss=True, inferbat=True, spatial=True)

    steps = basil.basil_steps(wsp, img)
    assert(len(steps) == 2)
    
    options = _get_defaults(img)
    _check_step(steps[0], desc_text="tissue", options=options)

    options.update({
        "method" : "spatialvb",
        "param-spatial-priors" : "N+",
        "PSP_byname1" : "ftiss",
        "PSP_byname1_type" : "M",
        "convergence" : "maxits",
    })
    options.pop("max-trials")
    _check_step(steps[1], desc_text="spatial", options=options)
    
def test_onestep():
    """
    Check that single step mode works when you would normally get multiple steps
    """
    d = np.random.rand(5, 5, 5, 6)
    img = AslImage(name="asldata", image=d, tis=[1.5], order="prt")
    wsp = Workspace(infertiss=True, inferbat=True, infertau=True, inferart=True, spatial=True, onestep=True)

    steps = basil.basil_steps(wsp, img)
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
        "convergence" : "maxits",
    })
    options.pop("max-trials")
    _check_step(steps[0], desc_text="spatial", options=options)

def test_max_iterations():
    """
    Check that max iterations can be overridden
    """
    d = np.random.rand(5, 5, 5, 6)
    img = AslImage(name="asldata", image=d, tis=[1.5], order="prt")
    wsp = Workspace(infertiss=True, inferbat=True)

    kwargs = {
        "max-iterations" : 123,
    }
    steps = basil.basil_steps(wsp, img, **kwargs)
    assert(len(steps) == 1)

    options = _get_defaults(img)
    options.update({
        "max-iterations" : 123,
    })
    _check_step(steps[0], desc_text="tissue", options=options)

def test_random_extra_options():
    """
    Check that any additional keyword arguments are passed to Fabber
    """
    d = np.random.rand(5, 5, 5, 6)
    img = AslImage(name="asldata", image=d, tis=[1.5], order="prt")
    wsp = Workspace(infertiss=True, inferbat=True)

    kwargs = {
        "phase-of-moon-correction-factor" : 7,
        "random-output-proportion-percent" : 36,
    }
    steps = basil.basil_steps(wsp, img, **kwargs)
    assert(len(steps) == 1)

    options = _get_defaults(img)
    options.update(kwargs)
    _check_step(steps[0], desc_text="tissue", options=options)

def test_pvc_only_one_map_given1():
    """
    Check that PVC correction fails if you only give the GM map
    """
    d = np.random.rand(5, 5, 5, 6)
    img = AslImage(name="asldata", image=d, tis=[1.5], order="prt")
    wsp = Workspace(infertiss=True, inferbat=True)

    wsp.pgm = Image(np.random.rand(5, 5, 5))
    with pytest.raises(ValueError):
        basil.basil_steps(wsp, img)
    
def test_pvc_only_one_map_given2():
    """
    Check that PVC correction fails if you only give the WM map
    """
    d = np.random.rand(5, 5, 5, 6)
    img = AslImage(name="asldata", image=d, tis=[1.5], order="prt")
    wsp = Workspace(infertiss=True, inferbat=True)

    wsp.pwm = Image(np.random.rand(5, 5, 5))
    with pytest.raises(ValueError):
        basil.basil_steps(wsp, img)
    
def test_pvc_no_tissue():
    """
    Check that PVC correction fails if you do not infer the tissue component
    """
    d = np.random.rand(5, 5, 5, 6)
    img = AslImage(name="asldata", image=d, tis=[1.5], order="prt")
    wsp = Workspace(infertiss=False, inferbat=True)

    wsp.pgm = Image(np.random.rand(5, 5, 5))
    wsp.pwm = Image(np.random.rand(5, 5, 5))

    with pytest.raises(ValueError):
        basil.basil_steps(wsp, img)
    
def test_pvc():
    """
    FIXME we need to test the PVC initialization step
    and how to do this is not finalized
    """
    d = np.random.rand(5, 5, 5, 6)
    img = AslImage(name="asldata", image=d, tis=[1.5], order="prt")
    wsp = Workspace(infertiss=True, inferbat=True)

    wsp.pgm = Image(np.random.rand(5, 5, 5))
    wsp.pwm = Image(np.random.rand(5, 5, 5))
    
    steps = basil.basil_steps(wsp, img)
    assert(len(steps) == 3)

    options = _get_defaults(img)
#    options.update({
#        "incpve" : True,
#    })
    _check_step(steps[0], desc_text="tissue", options=options)
#    _check_step(steps[1], desc_text="PVE", options=options)

    options.update({
        "method" : "spatialvb",
        "param-spatial-priors" : "N+",
        "PSP_byname1" : "ftiss",
        "PSP_byname1_type" : "M",
        "max-iterations" : 200,
        "convergence" : "maxits",
    })
    options.pop("max-trials")
    _check_step(steps[2], desc_text="spatial")

