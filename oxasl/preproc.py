"""
OXASL - ASL data preprocessing

Copyright (c) 2008-2020 Univerisity of Oxford
"""
from __future__ import print_function

import sys

import fsl.wrappers as fsl
from fsl.data.image import Image

import scipy
import numpy as np

from oxasl import Workspace, image
from oxasl.options import AslOptionParser, OptionCategory, OptionGroup, GenericOptions

def run(wsp):
    wsp.log.write("\nPre-processing input images\n")
    wsp.sub("preproc")
    asldata_reorient = _reorient(wsp, wsp.input.asldata)
    wsp.preproc.asldata = wsp.input.asldata.derived(asldata_reorient.data)
    wsp.preproc.aslspace = wsp.preproc.asldata.mean()
    try:
        wsp.preproc.pwi = wsp.asldata.perf_weighted()
    except:
        # Ignore - not all data can generate a PWI
        pass

    if wsp.calib_first_vol and wsp.input.calib is None:
        wsp.input.calib = wsp.asldata.calib

    wsp.preproc.calib = _single_volume(wsp, wsp.input.calib)
    wsp.preproc.cref = _single_volume(wsp, wsp.input.cref)
    wsp.preproc.cact = _single_volume(wsp, wsp.input.cact)
    wsp.preproc.cblip = _single_volume(wsp, wsp.input.cblip)

def _reorient(wsp, img):
    if wsp.noreorient:
        wsp.log.write(" - NOT reorienting input data to standard orientation - need to check registration is OK\n")
        return img
    else:
        wsp.log.write(" - Reorienting to standard orientation: %s\n" % img.name)
        ret = fsl.fslreorient2std(img, output=fsl.LOAD)
        return ret["output"]

def _single_volume(wsp, img, moco=True, discard_first=True):
    """
    Convert a potentially 4D image into a single 3D volume

    :param moco: If True, perform basic motion correction
    :param discard_first: If True, discard first volume if nvols > 1

    """
    if img is not None:
        wsp.log.write(" - Pre-processing image: %s\n" % img.name)
        if img.ndim == 4:
            if discard_first and img.shape[3] > 1:
                wsp.log.write("   - Removing first volume to ensure data is in steady state\n")
                img = Image(img.data[..., :-1], header=img.header)

            if moco and img.shape[3] > 1:
                if moco:
                    wsp.log.write("   - Motion correcting\n")
                    img = fsl.mcflirt(img, out=fsl.LOAD, log=wsp.fsllog)["out"]

            wsp.log.write("   - Taking mean across volumes\n")
            img = Image(np.mean(img.data, axis=-1), header=img.header)

        return _reorient(wsp, img)
    else:
        return None

def preprocess(wsp):
    """
    Basic preprocessing of ASL data

    :param wsp: Workspace object

    Required workspace attributes
    -----------------------------

     - ``asldata`` : AslImage containing data to preprocess


    Optional workspace attributes
    -----------------------------

     - ``diff`` : If True, perform tag/control subtraction
     - ``reorder`` : If specified, reorder the data using this ordering string
     - ``mc`` : If True, perform motion correction using mcflirt
     - ``smooth`` : If True, perform smoothing
     - ``fwhm`` : If smooth=True, the full-width-half-maximum for the smoothing kernel
     - ``ref`` : If mc=True, a reference image for motion correction (default is to use the middle volume)

    Workspace attributes updated
    -----------------------------

     - ``asldata_preproc`` : preprocessed AslImage instance
    """
    wsp.log.write("ASL preprocessing...\n")

    wsp.asldata_preproc = wsp.asldata

    if wsp.diff:
        if wsp.asldata.iaf == "hadamard":
            wsp.log.write("  - Hadamard decoding\n")
        else:
            wsp.log.write("  - Tag-control subtraction\n")
        wsp.asldata_preproc = wsp.asldata_preproc.diff()

    if wsp.reorder:
        wsp.log.write("  - Re-ordering to %s\n" % wsp.reorder)
        if "l" in wsp.reorder.lower() and wsp.diff:
            wsp.reorder = wsp.reorder.replace("l", "")
        wsp.asldata_preproc = wsp.asldata_preproc.reorder(wsp.reorder)

    if wsp.mc:
        wsp.log.write("  - Motion correction\n")
        mcimg = fsl.mcflirt(wsp.asldata_preproc, cost="mutualinfo", out=fsl.LOAD, reffile=wsp.ref)["out"]
        wsp.asldata_preproc = wsp.asldata_preproc.derived(mcimg.data, suffix="_mc")

    if wsp.smooth:
        wsp.sigma = round(wsp.fwhm / 2.355, 2)
        wsp.log.write("  - Spatial smoothing with FWHM: %f (sigma=%f)\n" % (wsp.fwhm, wsp.sigma))
        smoothed = scipy.ndimage.gaussian_filter(wsp.asldata_preproc.data, sigma=wsp.sigma)
        wsp.asldata_preproc = wsp.asldata_preproc.derived(smoothed, suffix="_smooth")

    wsp.log.write("DONE\n\n")

class AslPreprocOptions(OptionCategory):
    """
    OptionCategory which contains options for preprocessing ASL data
    """

    def __init__(self, **kwargs):
        OptionCategory.__init__(self, "preproc")

    def groups(self, parser):
        group = OptionGroup(parser, "Preprocessing")
        group.add_option("--diff", help="Perform tag-control subtraction / Hadamard decoding", action="store_true", default=False)
        group.add_option("--smooth", help="Spatially smooth data", action="store_true", default=False)
        group.add_option("--fwhm", help="FWHM for spatial filter kernel", type="float", default=6)
        group.add_option("--mc", help="Motion correct data", action="store_true", default=False)
        group.add_option("--ref", help="Optional reference image for motion correction", type="image", default=None)
        group.add_option("--reorder", help="Re-order data in specified order")
        return [group, ]

def main():
    """
    Entry point for command line tool
    """
    try:
        parser = AslOptionParser(usage="asl_preproc -i <filename> [options]")
        parser.add_category(GenericOptions(output_type="file"))
        parser.add_category(image.Options())
        parser.add_category(AslPreprocOptions())

        options, _ = parser.parse_args()
        if not options.output:
            options.output = "oxasl_preproc_out"
        wsp = Workspace(auto_asldata=True, **vars(options))

        print("")
        wsp.asldata.summary()
        preprocess(wsp)
        wsp.asldata_preproc.save(options.output)

    except ValueError as exc:
        sys.stderr.write("ERROR: " + str(exc) + "\n")
        sys.exit(1)

if __name__ == "__main__":
    main()
