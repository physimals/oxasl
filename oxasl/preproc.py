"""
ASL data preprocessing command line tool
"""
from __future__ import print_function

import sys

import scipy

import fsl.wrappers as fsl

from .options import AslOptionParser, OptionCategory, IgnorableOptionGroup
from .image import AslImage, AslImageOptions

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
        smoothed = scipy.ndimage.gaussian_filter(asldata.data, sigma=sigma)
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
        group.add_option("--ref", help="Optional reference image for motion correction", type="image", default=None)
        group.add_option("--reorder", help="Re-order data in specified order")
        return [group, ]

def main():
    """
    Entry point for command line tool
    """
    try:
        parser = AslOptionParser(usage="asl_preproc -i <filename> [options]")
        parser.add_option("--output", "-o", help="Output file", default=None)
        parser.add_option("--debug", help="Debug mode", action="store_true", default=False)
        
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
