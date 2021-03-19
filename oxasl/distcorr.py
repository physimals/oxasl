#!/bin/env python
"""
OXASL - Distortion correction module

The functions in this module calculate non-linear distortion correction transformations
to apply to the ASL and calibration images.

Currently two sources of transformation exist:

 - Fieldmap-based distortion correction. This generates a nonlinear warp image
   in structural space which is then transformed to ASL space

 - Phase encoding reversed (CBLIP) distortion correction using TOPUP. This generates
   a nonlinear warp image in ASL space

Copyright (c) 2008-2020 Univerisity of Oxford
"""
from __future__ import unicode_literals

import os
import sys
import tempfile
import shutil

import numpy as np

import fsl.wrappers as fsl
from fsl.data.image import Image

from oxasl import reg, struc
from oxasl.options import OptionCategory, OptionGroup
from oxasl.reporting import LightboxImage, LineGraph
from oxasl.wrappers import epi_reg, fnirtfileutils

def add_options(parser):
    g = OptionGroup(parser, "Distortion correction using fieldmap")
    g.add_option("--fmap", help="fieldmap image (in rad/s)", type="image")
    g.add_option("--fmapmag", help="fieldmap magnitude image - wholehead extracted", type="image")
    g.add_option("--fmapmagbrain", help="fieldmap magnitude image - brain extracted", type="image")
    g.add_option("--nofmapreg", help="Do not perform registration of fmap to T1 (use if fmap already in T1-space)", action="store_true", default=False)
    parser.add_option_group(g)
    
    g = OptionGroup(parser, "Distortion correction using phase-encode-reversed calibration image (TOPUP)")
    g.add_option("--cblip", help="phase-encode-reversed (blipped) calibration image", type="image")
    parser.add_option_group(g)
    
    g = OptionGroup(parser, "General distortion correction options")
    g.add_option("--echospacing", help="Effective EPI echo spacing (sometimes called dwell time) - in seconds", type=float)
    g.add_option("--pedir", help="Phase encoding direction, dir = x/y/z/-x/-y/-z")
    g.add_option("--gdc-warp", "--gdcwarp", help="Additional warp image for gradient distortion correction - will be combined with fieldmap or TOPUP distortion correction", type="image")
    parser.add_option_group(g)
        
def run(wsp):
    wsp.sub("distcorr")
    wsp.log.write("\nCalculating distortion corrections\n")
    get_fieldmap_correction(wsp.distcorr)
    get_cblip_correction(wsp.distcorr)
    get_gdc_correction(wsp.distcorr)

def get_cblip_correction(wsp):
    """
    Get the cblip based distortion correction warp

    Required workspace attributes
    -----------------------------

     - ``calib`` : Calibration image
     - ``cblip`` : Phase-encode-reversed calibration image
     - ``echospacing`` :
     - ``pedir`` :

    Updated workspace attributes
    ----------------------------

     - ``cblip_warp``    : CBLIP Distortion correction warp image
    """
    if wsp.topup is not None:
        return
    elif wsp.cblip is None:
        wsp.log.write(" - No CBLIP images provided for distortion correction\n")
        return

    wsp.sub("topup")
    wsp.log.write(" - Calculating distortion Correction using TOPUP\n")

    topup_params = {
        "x"  : [[1, 0, 0, -99], [-1, 0, 0, -99]],
        "-x" : [[-1, 0, 0, -99], [1, 0, 0, -99]],
        "y"  : [[0, 1, 0, -99], [0, -1, 0, -99]],
        "-y" : [[0, -1, 0, -99], [0, 1, 0, -99]],
        "z"  : [[0, 0, 1, -99], [0, 0, -1, -99]],
        "-z" : [[0, 0, -1, -99], [0, 0, 1, -99]],
    }
    dim_idx = {
        "x"  : 0, "-x" : 0,
        "y"  : 1, "-y" : 1,
        "z"  : 2, "-z" : 2,
    }
    my_topup_params = np.array(topup_params[wsp.pedir], dtype=np.float)
    dimsize = wsp.asldata.shape[dim_idx[wsp.pedir]]
    my_topup_params[:, 3] = wsp.echospacing * (dimsize - 1)
    wsp.topup.params = my_topup_params

    # Run TOPUP to calculate correction
    wsp.topup.calib_blipped = reg.change_space(wsp, Image(np.stack((wsp.calib.data, wsp.cblip.data), axis=-1), header=wsp.calib.header), 'native')
    topup_result = fsl.topup(imain=wsp.topup.calib_blipped, datain=wsp.topup.params, config="b02b0.cnf", out=fsl.LOAD, iout=fsl.LOAD, fout=fsl.LOAD, log=wsp.fsllog)
    wsp.topup.fieldcoef, wsp.topup.movpar = topup_result["out_fieldcoef"], topup_result["out_movpar"]
    wsp.topup.iout = topup_result["iout"]
    wsp.topup.fout = topup_result["fout"]

    page = wsp.report.page("topup")
    page.heading("TOPUP distortion correction", level=0)
    page.text("PE direction: %s" % wsp.pedir)
    page.text("Echo spacing: %f s" % wsp.echospacing)
    page.heading("Correction image", level=1)
    for dim in range(3):
        img = Image(wsp.topup.fieldcoef.data[..., dim], header=wsp.topup.fieldcoef.header)
        page.text("Dimension %i" % dim)
        page.image("fmap_warp%i" % dim, LightboxImage(img))

def get_fieldmap_correction(wsp):
    """
    Get the fieldmap based distortion correction warp

    Required workspace attributes
    -----------------------------

     - ``pwi``          : Perfusion weighted image (generated by preproc_asl)
     - ``fmap``         : Fieldmap image
     - ``fmapmag``      : Fieldmap magnitude image
     - ``fmapmagbrain`` : Fieldmap magnitude brain image
     - ``echospacing``  :
     - ``pedir``        :

    Optional workspace attributes
    -----------------------------

     - ``nofmapreg``       : If True assume fieldmap in structural space

    Updated workspace attributes
    ----------------------------

     - ``fmap_warp``    : Fieldmap distortion correction warp image in ASL space
    """
    if wsp.fmap is None or wsp.fmapmag is None or wsp.fmapmagbrain is None:
        wsp.log.write(" - No fieldmap images for distortion correction\n")
        return
    elif wsp.pedir is None or wsp.echospacing is None:
        wsp.log.write(" -WARNING: Fieldmap images supplied but pedir and echospacing required for distortion correction\n")
        return

    wsp.sub("fieldmap")
    wsp.log.write(" - Calculating distortion correction from fieldmap images using EPI_REG\n")

    epi_reg_opts = {
        "inweight" : wsp.inweight,
        "init" : wsp.reg.asl2struc,
        "fmap" : wsp.fmap,
        "fmapmag" : wsp.fmapmag,
        "fmapmagbrain" : wsp.fmapmagbrain,
        "pedir" : wsp.pedir,
        "echospacing" : wsp.echospacing,
        "nofmapreg" : wsp.ifnone("nofmapreg", False),
    }

    # Windows can't run epi_reg as it's a batch script. Use our experimental python
    # implementation but use the standard epi_reg on other platforms until the python
    # version is better tested
    if sys.platform.startswith("win"):
        import oxasl.epi_reg as pyepi
        result = pyepi.epi_reg(wsp, epi=wsp.reg.nativeref, **epi_reg_opts)
    else:
        result = epi_reg(epi=wsp.asldata.perf_weighted(), t1=wsp.structural.struc, t1brain=wsp.structural.brain, out=fsl.LOAD, wmseg=wsp.structural.wm_seg, log=wsp.fsllog, **epi_reg_opts)

    # Occasionally we end up with NaN in the output of epi_reg and this will ruin the entire distcorr warp.
    # So remove any NaN values and replace with zero.
    warp_struc = result["out_warp"]
    wsp.fieldmap.warp_struc = Image(np.nan_to_num(warp_struc.data, nan=0), header=warp_struc.header)
    wsp.fieldmap.asl2struc = result["out"]
    wsp.fieldmap.struc2asl = np.linalg.inv(wsp.fieldmap.asl2struc)

    result = fsl.convertwarp(out=fsl.LOAD, ref=wsp.reg.nativeref, warp1=wsp.fieldmap.warp_struc, postmat=wsp.fieldmap.struc2asl, rel=True, log=wsp.fsllog)
    wsp.fieldmap.warp = result["out"]

    page = wsp.report.page("fmap")
    page.heading("Fieldmap distortion correction", level=0)
    page.text("PE direction: %s" % wsp.pedir)
    page.text("Echo spacing: %f s" % wsp.echospacing)
    page.heading("Correction warps", level=1)
    for dim in range(3):
        img = Image(wsp.fieldmap.warp.data[..., dim], header=wsp.fieldmap.warp.header)
        page.text("Dimension %i" % dim)
        page.image("fmap_warp%i" % dim, LightboxImage(img))

def get_gdc_correction(wsp):
    """
    User-specified gradient distortion correction warp assumed to be in ASL space
    """
    wsp.gdc_warp = wsp.gdc_warp
