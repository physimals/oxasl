import StringIO

import numpy as np

from oxasl import basil, AslImage

"""
def get_steps(asldata, mask, 
              infertau=False, inferart=False, infert1=False, inferpc=False,
              artonly=False, fixbat=False, spatial=False, onestep=False,
              t1im=None, pgm=None, pwm=None,
              initmvn=None,
              log=sys.stdout, **kwargs):

def _do_step(wsp, step, step_desc, infile, mask, options, prev_step=None, log=sys.stdout):
"""

def _check_step(step, step_num=None, desc_text=None, data_name=None, options=None, prev_step=None):
    if step_num:
        assert(step[0] == step_num)
    if desc_text:
        assert(desc_text.lower().strip() in step[1].lower())
    if data_name:
        assert(step[2].iname == data_name)
    if options:
        for k, v in options.items():
            assert(step[4][k] == v)
    if prev_step:
        assert(step[5] == prev_step)

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

def test_defaults():
    d = np.random.rand(5, 5, 5, 6)
    img = AslImage("asldata", data=d, tis=[1.5], order="prt")
    log = StringIO.StringIO()

    options = _get_defaults(img)
    steps = basil.get_steps(img, log=log)
    assert(len(steps) == 1)

    _check_step(steps[0], step_num=1, desc_text="tissue", 
                options=options)
                
def test_fixbat():
    d = np.random.rand(5, 5, 5, 6)
    img = AslImage("asldata", data=d, tis=[1.5], order="prt")
    log = StringIO.StringIO()

    options = _get_defaults(img)
    steps = basil.get_steps(img, fixbat=True, log=log)
    assert(len(steps) == 1)

    options.pop("incbat")
    options.pop("inferbat")
    _check_step(steps[0], step_num=1, desc_text="tissue", 
                options=options)
 
def test_inferart():
    d = np.random.rand(5, 5, 5, 6)
    img = AslImage("asldata", data=d, tis=[1.5], order="prt")
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
    d = np.random.rand(5, 5, 5, 6)
    img = AslImage("asldata", data=d, tis=[1.5], order="prt")
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
    d = np.random.rand(5, 5, 5, 6)
    img = AslImage("asldata", data=d, tis=[1.5], order="prt")
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
    d = np.random.rand(5, 5, 5, 6)
    img = AslImage("asldata", data=d, tis=[1.5], order="prt")
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

def test_inferpc():
    d = np.random.rand(5, 5, 5, 6)
    img = AslImage("asldata", data=d, tis=[1.5], order="prt")
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

def test_spatial():
    d = np.random.rand(5, 5, 5, 6)
    img = AslImage("asldata", data=d, tis=[1.5], order="prt")
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
#        "convergence" : "maxiters",
    })
    options.pop("max-trials")
    _check_step(steps[1], step_num=2, desc_text="spatial", options=options)
    