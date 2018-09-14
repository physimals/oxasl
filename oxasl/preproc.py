"""
ASL data preprocessing command line tool
"""
import sys

import numpy as np
import scipy

import fsl.wrappers as fsl

from .options import AslOptionParser, OptionCategory, IgnorableOptionGroup
from .image import AslImage, AslImageOptions
from .reporting import ReportPage

def motion_correct(wsp):
    """
    Motion Correction of ASL data
    
    Note motion correction of multi-volume calibration data is done in preprocessing.

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
     - ``cref``     : Calibration reference image
     - ``cblip``     : Calibration BLIP image

    Updated workspace attributes
    ----------------------------

     - ``asldata_mc_mats`` : Sequence of matrices giving motion correction transform for each ASL volume
     - ``asl2calib``       : ASL->calibration image transformation
     - ``calib2asl``       : Calibration->ASL image transformation
    """
    wsp.log.write("\nMotion Correction\n")
    if wsp.calib:
        # use supplied image as our reference for motion correction
        # Normally the calibration image since this will be most consistent if the data has a range of different TIs and background suppression etc
        # this also removes motion effects between asldata and calibration image
        wsp.log.write(" - Using calibration image as reference\n")
        ref_source = "calibration image: %s" % wsp.calib.name
        mcflirt_result = fsl.mcflirt(wsp.asldata, reffile=wsp.calib, out=fsl.LOAD, mats=fsl.LOAD, log=wsp.fsllog)
        wsp.asldata_mc_mats = [mcflirt_result["out.mat/MAT_%04i" % vol] for vol in range(wsp.asldata.shape[3])]

        # To reduce interpolation of the ASL data change the transformations so that we end up in the space of the central volume of asldata
        wsp.asl2calib = wsp.asldata_mc_mats[len(wsp.asldata_mc_mats)/2+1]
        wsp.calib2asl = np.linalg.inv(wsp.asl2calib)
        wsp.asldata_mc_mats = [np.dot(mat, wsp.calib2asl) for mat in wsp.asldata_mc_mats]
        wsp.log.write("   ASL middle volume->Calib:\n%s\n" % str(wsp.asl2calib))
        wsp.log.write("   Calib->ASL middle volume:\n%s\n" % str(wsp.calib2asl))
    else:
        wsp.log.write(" - Using ASL data middle volume as reference\n")
        ref_source = "ASL data %s middle volume: %s" % wsp.asldata.name
        mcflirt_result = fsl.mcflirt(wsp.asldata, out=fsl.LOAD, mats=fsl.LOAD, log=wsp.fsllog)
        wsp.asldata_mc_mats = [mcflirt_result["out.mat/MAT_%04i" % vol] for vol in range(wsp.asldata.shape[3])]
        
    page = ReportPage()
    page.heading("Motion correction", level=0)
    page.text("Reference volume: %s" % ref_source)
    page.heading("Motion parameters", level=1)
    for vol, mat in enumerate(wsp.asldata_mc_mats):
        rst_math = "\\begin{bmatrix}\n"
        for row in mat:
            rst_math += "    " + " & ".join([str(v) for v in row]) + " \\"
        rst_math += "\\end{bmatrix}\n"
        page.maths(rst_math)
    wsp.report.add("moco", page)

def preproc_asl(wsp):
    """
    Preprocessing on the main ASL data - calculate averaged images and run brain extraction
    """
    if wsp.asldata is not None and wsp.diffasl is None:
        wsp.log.write("\nPre-processing ASL data: %s\n" % wsp.asldata.name)
        wsp.asldata_orig = wsp.asldata
        wsp.asldata_mean_across_repeats = wsp.asldata.mean_across_repeats()

        if "p" in wsp.asldata.order.lower():
            # Data contains tag-control pairs to can extract mean raw signal
            wsp.asldata_mean = wsp.asldata.mean()
            bet_result = fsl.bet(wsp.asldata_mean, seg=True, mask=True, fracintensity=0.2, output=fsl.LOAD, log=wsp.fsllog)
            wsp.asldata_mean_brain = bet_result["output"]
            wsp.asldata_mean_brain_mask = bet_result["output_mask"]

        wsp.diffasl = wsp.asldata.diff()
        wsp.diffasl_mean = wsp.diffasl.mean()
        bet_result = fsl.bet(wsp.diffasl_mean, seg=True, mask=True, fracintensity=0.2, output=fsl.LOAD, log=wsp.fsllog)
        wsp.diffasl_mean_brain = bet_result["output"]

        wsp.pwi = wsp.asldata.perf_weighted()
        bet_result = fsl.bet(wsp.pwi, seg=True, mask=True, fracintensity=0.2, output=fsl.LOAD, log=wsp.fsllog)
        wsp.pwi_brain = bet_result["output"]
        wsp.pwi_brain_mask = bet_result["output_mask"]
        wsp.log.write(" - Generated subtracted and averaged copies of input data\n")

def preprocess(asldata, diff=False, reorder=None, mc=False, smooth=False, fwhm=None, ref=None, log=sys.stdout):
    """
    Basic preprocessing of ASL data

    :param asldata: AslImage instance
    :param diff: If True, perform tag/control subtraction
    :param reorder: If specified, reorder the data using this ordering string
    :param mc: If True, perform motion correction using mcflirt
    :param smooth: If True, perform smoothing
    :param fwhm: If smooth=True, the full-width-half-maximum for the smoothing kernel
    :param ref: If mc=True, a reference image for motion correction (default is to use the middle volume)
    :param log: Output log stream (default: stdout)

    :return: AslImage instance
    """
    log.write("ASL preprocessing...\n")

    # Keep original AslImage with info about TIs, repeats, etc
    orig = asldata

    if diff: 
        log.write("  - Tag-control subtraction\n")
        asldata = asldata.diff()
        
    if reorder:
        log.write("  - Re-ordering to %s\n" % reorder)
        if "l" in reorder.lower() and diff:
            reorder = reorder.replace("l", "")
        asldata = asldata.reorder(reorder)

    if mc: 
        log.write("  - Motion correction\n")
        mcimg = fsl.mcflirt(asldata, cost="mutualinfo", out=fsl.LOAD, reffile=ref)["out"]
        asldata = orig.derived(mcimg.data, suffix="_mc")

    if smooth:
        sigma = round(fwhm/2.355, 2)
        log.write("  - Spatial smoothing with FWHM: %f (sigma=%f)\n" % (fwhm, sigma))
        smoothed = scipy.ndimage.gaussian_filter(asldata.nibImage.get_data(), sigma=sigma)
        asldata = asldata.derived(smoothed, suffix="_smooth")

    log.write("DONE\n\n")
    return asldata

class AslPreprocOptions(OptionCategory):
    """
    OptionCategory which contains options for preprocessing ASL data
    """

    def __init__(self, **kwargs):
        OptionCategory.__init__(self, "preproc", **kwargs)

    def groups(self, parser):
        group = IgnorableOptionGroup(parser, "Preprocessing", ignore=self.ignore)
        group.add_option("--diff", help="Perform tag-control subtraction", action="store_true", default=False)
        group.add_option("--smooth", help="Spatially smooth data", action="store_true", default=False)
        group.add_option("--fwhm", help="FWHM for spatial filter kernel", type="float", default=6)
        group.add_option("--mc", help="Motion correct data", action="store_true", default=False)
        group.add_option("--ref", help="Optional reference image for motion correction", default=None)
        group.add_option("--reorder", help="Re-order data in specified order")
        return [group, ]

def main():
    """
    Entry point for command line tool
    """
    try:
        parser = AslOptionParser(usage="asl_preproc -i <filename> [options]")
        parser.add_option("-o", dest="output", help="Output file", default=None)
        parser.add_option("--debug", dest="debug", help="Debug mode", action="store_true", default=False)
        
        parser.add_category(AslImageOptions())
        parser.add_category(AslPreprocOptions())

        options, _ = parser.parse_args(sys.argv)
        kwopts = vars(options)

        if not options.asldata:
            sys.stderr.write("Input file not specified\n")
            parser.print_help()
            sys.exit(1)

        options.asldata = AslImage(options.asldata, role="Input data", **parser.filter(kwopts, "image"))
        if options.output is None:
            options.output = options.asldata.name + "_out"
        options.asldata.summary()
        print("")

        data_preproc = preprocess(options.asldata, **parser.filter(kwopts, "preproc"))
        data_preproc.save(options.output)

    except ValueError as exc:
        sys.stderr.write("ERROR: " + str(exc) + "\n")
        sys.exit(1)

if __name__ == "__main__":
    main()
