#!/bin/env python
"""
OXASL - Module to apply Moco/Distortion/sensitivity corrections

This module generates corrected ASL/calibration data from previously calculated
corrections with the minimum of interpolation

Currently the following sources of transformation exist:

 - Motion correction of the ASL data. This generates a series of linear (rigid body)
   transformations in ASL space, one for each ASL volume. If calibration data is also
   present a calibration->ASL transform is also generated as part of this process

 - Fieldmap-based distortion correction. This generates a nonlinear warp image
   in structural space which is then transformed to ASL space

 - Phase encoding reversed (CBLIP) distortion correction using TOPUP. This generates
   a nonlinear warp image in ASL space FIXME calibration space?

 - User-supplied nonlinear warp image for gradient distortion corection

 - Sensitivity correction

Except for the TOPUP correction, all of the above can be combined in a single
transformation to minimise interpolation of the ASL data

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
from oxasl.options import OptionCategory
from oxasl.reporting import LightboxImage, LineGraph
from oxasl.wrappers import epi_reg, fnirtfileutils

class Options(OptionCategory):
    """
    Options for corrections of the input data
    """

    def __init__(self):
        OptionCategory.__init__(self, "corrections")

    def groups(self, parser):
        ret = []
        return ret

def run(wsp):
    """
    Apply distortion and motion corrections to ASL and calibration data

    Required workspace attributes
    -----------------------------

     - ``asldata_orig`` : Uncorrected ASL data image

    Optional workspace attributes
    -----------------------------

     - ``calib_orig``      : Calibration image
     - ``cref_orig``       : Calibration reference image
     - ``cblip_orig``      : Calibration BLIP image
     - ``asldata_mc_mats`` : ASL motion correction matrices
     - ``calib2asl``       : Calibration -> ASL transformation matrix
     - ``distcorr_warp``   : Distortion correction warp image
     - ``gdc_warp``        : Gradient distortion correction warp image

    Updated workspace attributes
    ----------------------------

     - ``asldata``    : Corrected ASL data
     - ``calib``      : Corrected calibration image
     - ``cref``       : Corrected calibration reference image
     - ``cblip``      : Corrected calibration BLIP image
    """
    wsp.sub("corrected")
    wsp.log.write("\nApplying data corrections\n")

    warps, moco_mats = [], None
    
    if wsp.moco is not None:
        wsp.log.write(" - Using motion correction\n")
        moco_mats = wsp.moco.mc_mats

    if wsp.distcorr is not None:
        if wsp.distcorr.fieldmap is not None:
            wsp.log.write(" - Using fieldmap distortion correction\n")
            warps.append(wsp.distcorr.fieldmap.warp)

        if wsp.distcorr.gdc_warp:
            wsp.log.write(" - Using user-supplied GDC warp\n")
            warps.distcorr.append(wsp.gdc_warp)

    if warps:
        kwargs = {}
        for idx, warp in enumerate(warps):
            kwargs["warp%i" % (idx+1)] = warp

        wsp.log.write(" - Converting all warps to single transform and extracting Jacobian\n")
        result = fsl.convertwarp(ref=wsp.reg.nativeref, out=fsl.LOAD, rel=True, jacobian=fsl.LOAD, log=wsp.fsllog, **kwargs)
        wsp.corrected.total_warp = result["out"]

        # Calculation of the jacobian for the warp - method suggested in:
        # https://www.jiscmail.ac.uk/cgi-bin/webadmin?A2=FSL;d3fee1e5.0908
        wsp.corrected.warp_coef = fnirtfileutils(wsp.corrected.total_warp, outformat="spline", out=fsl.LOAD, log=wsp.fsllog)["out"]
        jacobian = fnirtfileutils(wsp.corrected.warp_coef, jac=fsl.LOAD, log=wsp.fsllog)["jac"]
        wsp.corrected.jacobian = Image(jacobian.data, header=wsp.corrected.total_warp.header)

    if not warps and moco_mats is None and wsp.senscorr is None:
        wsp.log.write(" - No corrections to apply to ASL data\n")
        wsp.corrected.asldata = wsp.preproc.asldata
    else:
        # Apply all corrections to ASL data - note that we make sure the output keeps all the ASL metadata
        wsp.log.write(" - Applying corrections to ASL data\n")
        asldata_corr = correct_img(wsp, wsp.preproc.asldata, moco_mats)
        wsp.corrected.asldata = wsp.preproc.asldata.derived(asldata_corr.data)

    if wsp.preproc.calib is not None:
        # Apply corrections to calibration images if we have calib2asl registration or any other correction
        if not warps and (wsp.reg is None or wsp.reg.calib2asl is None) and wsp.senscorr is None:
            wsp.log.write(" - No corrections to apply to calibration data\n")
            wsp.corrected.calib = wsp.preproc.calib
            if wsp.cref is not None:
                wsp.corrected.cref = wsp.preproc.cref
            if wsp.cblip is not None:
                wsp.corrected.cblip = wsp.preproc.cblip
        else:
            wsp.log.write(" - Applying corrections to calibration data\n")
            wsp.corrected.calib = correct_img(wsp, wsp.preproc.calib, wsp.reg.calib2asl)

            if wsp.cref is not None:
                wsp.corrected.cref = correct_img(wsp, wsp.preproc.cref, wsp.reg.calib2asl)
            if wsp.cblip is not None:
                wsp.corrected.cblip = correct_img(wsp, wsp.preproc.cblip, wsp.reg.calib2asl)

    if wsp.distcorr is not None and wsp.distcorr.topup is not None:
        wsp.log.write(" - Adding TOPUP distortion correction\n")
        # This can't currently be done using the FSL wrappers - we need the TOPUP output as two prefixed files
        # Only workaround currently is to create a temp directory to store appropriately named input files
        topup_input = tempfile.mkdtemp(prefix="topup_input")
        try:
            wsp.distcorr.topup.fieldcoef.save("%s/topup_fieldcoef" % topup_input)
            movpar_file = open("%s/topup_movpar.txt" % topup_input, "w")
            for row in wsp.distcorr.topup.movpar:
                movpar_file.write("\t".join([str(val) for val in row]) + "\n")
            movpar_file.close()
            # TOPUP does not do the jacobian magntiude correction - so only okay if using voxelwise calibration
            wsp.corrected.calib = fsl.applytopup(wsp.corrected.calib, datain=wsp.distcorr.topup.params, index=1, topup="%s/topup" % topup_input, out=fsl.LOAD, method="jac", log=wsp.fsllog)["out"]
            if wsp.cref is not None:
                wsp.corrected.cref = fsl.applytopup(wsp.corrected.cref, datain=wsp.distcorr.topup.params, index=1, topup="%s/topup" % topup_input, out=fsl.LOAD, method="jac", log=wsp.fsllog)["out"]
            if wsp.cblip is not None:
                wsp.corrected.cblip = fsl.applytopup(wsp.corrected.cblip, datain=wsp.distcorr.topup.params, index=2, topup="%s/topup" % topup_input, out=fsl.LOAD, method="jac", log=wsp.fsllog)["out"]
            post_topup = fsl.applytopup(wsp.corrected.asldata, datain=wsp.distcorr.topup.params, index=1, topup="%s/topup" % topup_input, out=fsl.LOAD, method="jac", log=wsp.fsllog)["out"]
            wsp.corrected.asldata = wsp.corrected.asldata.derived(post_topup.data)
            # FIXME warning below
            # FIXME do we need to correct anything else in ASL or calibration space, e.g. mask, reference region mask
            #if wsp.calib_method != "voxel":
            #    wsp.log.write("WARNING: Using TOPUP does not correct for magntiude using the jocbian in distortion correction")
            #    wsp.log.write("         This is not optimal when not using voxelwise calibration\n")
            #    wsp.log.write("         To avoid this supply structural image(s)\n")
        finally:
            shutil.rmtree(topup_input)

    try:
        wsp.corrected.pwi = wsp.corrected.asldata.perf_weighted()
    except:
        # Ignore - not all data can generate a PWI
        pass

def correct_img(wsp, img, linear_mat):
    """
    Apply combined warp/linear transformations to an image in ASL space

    :param img: fsl.data.image.Image to correct
    :param linear_mat: img->ASL space linear transformation matrix.
    :return: Corrected Image

    If a jacobian is present, also corrects for quantitative signal magnitude as volume has been locally scaled

    FIXME there are slight differences to oxford_asl here due to use of spline interpolation rather than
    applyxfm4D which uses sinc interpolation.

    Required workspace attributesd
    -----------------------------

     - ``asldata_mean`` : Mean ASL image used as reference space

    Optional workspace attributes
    -----------------------------

     - ``total_warp``      : Combined warp image
     - ``jacobian``        : Jacobian associated with warp image
     - ``senscorr``        : Sensitivity correction
    """
    if wsp.corrected.total_warp is not None:
        img = reg.transform(wsp, img, trans=wsp.corrected.total_warp, ref=wsp.preproc.aslspace, premat=linear_mat)
    elif linear_mat is not None:
        img = reg.transform(wsp, img, trans=linear_mat, ref=wsp.preproc.aslspace)

    if wsp.corrected.jacobian is not None:
        wsp.log.write(" - Correcting for local volume scaling using Jacobian\n")
        jdata = wsp.corrected.jacobian.data
        if img.data.ndim == 4:
            # Required to make broadcasting work
            jdata = jdata[..., np.newaxis]
        img = Image(img.data * jdata, header=img.header)

    if wsp.senscorr is not None:
        wsp.log.write(" - Applying sensitivity correction\n")
        sens_data = reg.change_space(wsp, wsp.senscorr.sensitivity, img).data
        if img.ndim == 4:
            sens_data = sens_data[..., np.newaxis]
        img = Image(img.data / sens_data, header=img.header)

    return img
