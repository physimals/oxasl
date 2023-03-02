#!/bin/env python
"""
OXASL - Registration for ASL data

Michael Chappell, IBME QuBIc & FMRIB Image Analysis Groups

Copyright (c) 2008-2020 Univerisity of Oxford
"""
from __future__ import unicode_literals

import os

import numpy as np

import fsl.wrappers as fsl

from oxasl import reg
from oxasl.reporting import LineGraph

def run(wsp):
    """
    Calculate motion correction transforms for ASL data

    Note simple motion correction of multi-volume calibration data is done in preprocessing.

    The reference volume for motion correction is the calibration image, if supplied, or
    otherwise the middle volume of the ASL data is used.

    If the calibration image is used, the inverse of the middle ASL volume -> calibration
    transform is applied to each transform matrix. This ensures that the middle volume of
    the ASL data is unchanged and interpolation on the other volumes is also minimised.
    In this case, all calibration images are also themselves transformed to bring them in
    to ASL middle volume space.

    Required workspace attributes
    -----------------------------

     - ``asldata`` : ASL data image

    Optional workspace attributes
    -----------------------------

     - ``calib``    : Calibration image

    Updated workspace attributes
    ----------------------------

     - ``asldata_mc_mats`` : Sequence of matrices giving motion correction transform for each ASL volume
     - ``asl2calib``       : ASL->calibration image transformation
     - ``calib2asl``       : Calibration->ASL image transformation
    """
    if not wsp.mc:
        return

    wsp.sub("moco")
    wsp.log.write("\nCalculating Motion Correction\n")
    nvols = wsp.preproc.asldata.shape[3]
    # If available, use the calibration image as reference since this will be most consistent if the data has a range
    # of different TIs and background suppression etc. This also removes motion effects between asldata and calibration image
    if wsp.input.aslref is not None:
        wsp.log.write(" - Using user-specified aslref as reference\n")
        ref_source = "User specified: %s" % (wsp.input.aslref.name)
        wsp.moco.ref = wsp.input.aslref
        wsp.moco.input = wsp.preproc.asldata
        mcflirt_result = fsl.mcflirt(wsp.moco.input, reffile=wsp.moco.ref, out=fsl.LOAD, mats=fsl.LOAD, log=wsp.fsllog)
        mats = [mcflirt_result[os.path.join("out.mat", "MAT_%04i" % vol)] for vol in range(nvols)]
    elif wsp.preproc.calib is not None:
        wsp.log.write(" - Using calibration image as reference\n")
        ref_source = "Calibration image"
        wsp.moco.ref = wsp.preproc.calib
        wsp.moco.input = wsp.preproc.asldata
        mcflirt_result = fsl.mcflirt(wsp.moco.input, reffile=wsp.moco.ref, out=fsl.LOAD, mats=fsl.LOAD, log=wsp.fsllog)
        mats = [mcflirt_result[os.path.join("out.mat", "MAT_%04i" % vol)] for vol in range(nvols)]

        # To reduce interpolation of the ASL data change the transformations so that we end up in the space of the central volume of asldata
        wsp.moco.asl2calib = mats[int(float(len(mats))/2)]
        wsp.moco.calib2asl = np.linalg.inv(wsp.moco.asl2calib)
        mats = [np.dot(wsp.moco.calib2asl, mat) for mat in mats]
        wsp.log.write("   ASL middle volume->Calib:\n%s\n" % str(wsp.moco.asl2calib))
        wsp.log.write("   Calib->ASL middle volume:\n%s\n" % str(wsp.moco.calib2asl))
    else:
        wsp.log.write(" - Using ASL data middle volume as reference\n")
        ref_source = "ASL data middle volume: %i" % int(float(wsp.preproc.asldata.shape[3])/2)
        wsp.moco.input = wsp.preproc.asldata
        mcflirt_result = fsl.mcflirt(wsp.moco.input, out=fsl.LOAD, mats=fsl.LOAD, log=wsp.fsllog)
        mats = [mcflirt_result[os.path.join("out.mat", "MAT_%04i" % vol)] for vol in range(nvols)]

    # Convert motion correction matrices into single (4*nvols, 4) matrix - convenient for writing
    # to file, and same form that applywarp expects. Also save output ASL data + PWI
    wsp.moco.mc_mats = np.concatenate(mats, axis=0)
    output = reg.transform(wsp, wsp.moco.input, trans=wsp.moco.mc_mats, ref=wsp.preproc.aslspace)
    wsp.moco.output = wsp.asldata.derived(output)
    try:
        wsp.moco.pwi = wsp.moco.output.perf_weighted()
    except:
        # Ignore - not all data can generate a PWI
        pass

    page = wsp.report.page("moco")
    page.heading("Motion correction", level=0)
    page.heading("Reference volume", level=1)
    page.text(ref_source)
    page.heading("Motion parameters", level=1)
    moco_params = [reg.get_transform_params(mat) for mat in mats]
    trans = [p[0] for p in moco_params]
    abstrans = np.fabs(trans)
    rot = [p[1] for p in moco_params]
    absrot = np.fabs(rot)
    page.table([
        ["Mean translation", "%.3g mm" % np.mean(trans)],
        ["Translation std.dev.", "%.3g mm" % np.std(trans)],
        ["Absolute maximum translation", "%.3g mm (volume %i)" % (np.max(abstrans), np.argmax(abstrans))],
        ["Mean rotation", "%.3g \N{DEGREE SIGN}" % np.mean(rot)],
        ["Rotation std.dev.", "%.3g \N{DEGREE SIGN}" % np.std(rot)],
        ["Absolute maximum rotation", "%.3g \N{DEGREE SIGN} (volume %i)" % (np.max(absrot), np.argmax(absrot))],
    ])
    page.image("moco_trans", LineGraph(trans, "Volume number", "Translation (mm)"))
    page.image("moco_rot", LineGraph(rot, "Volume number", "Rotation relative to reference (\N{DEGREE SIGN})"))
