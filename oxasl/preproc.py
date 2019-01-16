"""
ASL data preprocessing command line tool
"""
from __future__ import print_function

import sys

import scipy

import fsl.wrappers as fsl

from oxasl import Workspace, image
from oxasl.options import AslOptionParser, OptionCategory, IgnorableOptionGroup, GenericOptions

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
        parser.add_category(GenericOptions(output_type="file", ignore=["mask"]))
        parser.add_category(image.AslImageOptions())
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
