"""
Wrapper for epi_reg command
"""

import fsl.utils.assertions as asrt
from fsl.wrappers import wrapperutils  as wutils

@wutils.fileOrImage('epi', 't1', 't1brain', 'fmap', 'fmapmag', 'fmapmagbrain', 'gdc', 'wmseg', 'weight', outprefix='out')
@wutils.fileOrArray('init', 'out.mat')
@wutils.fslwrapper
def epi_reg(**kwargs):
    """Wrapper for the ``epi_reg`` command.
    
    Required options:
    
    :arg epi: EPI image
    :arg t1: T1 weighted image (wholehead)
    :arg t1_brain: T1 weighted image (brain extracted)
    :arg out: Output image name or LOAD to return in-memory image
    
    Additional options:

    :arg fmap:          fieldmap image (in rad/s)
    :arg fmapmag:       fieldmap magnitude image - wholehead extracted
    :arg fmapmagbrain:  fieldmap magnitude image - brain extracted
    :arg gdc:           Gradient-distortion corection warpfield
    :arg wmseg:         white matter segmentation of T1 image
    :arg echospacing:   Effective EPI echo spacing (sometimes called dwell time) - in seconds
    :arg pedir:         Phase encoding direction, dir = x/y/z/-x/-y/-z
    :arg weight:        weighting image (in T1 space)
    :arg nofmapreg:     do not perform registration of fmap to T1 (use if fmap already registered) 
    """

    valmap = {
        'nofmapreg' : wutils.SHOW_IF_TRUE,
    }

    cmd = ['epi_reg', ]
    cmd += wutils.applyArgStyle('--=', valmap=valmap, singlechar_args=True, **kwargs)
    
    return cmd
