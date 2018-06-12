"""
ASL data preprocessing command line tool
"""
import sys
from optparse import OptionParser

import scipy

from fsl.wrappers import mcflirt, LOAD

from .image import AslImage, add_data_options, AslOptionGroup

def add_preproc_options(parser, ignore=()):
    g = AslOptionGroup(parser, "Preprocessing", ignore=ignore)
    g.add_option("--diff", dest="diff", help="Perform tag-control subtraction", action="store_true", default=False)
    g.add_option("--smooth", dest="smooth", help="Spatially smooth data", action="store_true", default=False)
    g.add_option("--fwhm", dest="fwhm", help="FWHM for spatial filter kernel", type="float", default=6)
    g.add_option("--mc", dest="mc", help="Motion correct data", action="store_true", default=False)
    g.add_option("--reorder", dest="reorder", help="Re-order data in specified order")
    parser.add_option_group(g)

def preprocess(asldata, diff=False, reorder=None, mc=False, smooth=False, fwhm=None, ref=None, log=sys.stdout, **kwargs):
    log.write("ASL preprocessing...\n")

    # Keep original AslImage with info about TIs, repeats, etc
    orig = asldata

    if diff: 
        log.write("  - Tag-control subtraction\n")
        asldata = asldata.diff()
        
    if reorder:
        log.write("  - Re-ordering to %s\n" % reorder)
        if "p" in reorder.lower() and diff:
            reorder = reorder.replace("p", "").replace("P", "") 
        asldata = asldata.reorder(reorder)

    if mc: 
        log.write("  - Motion correction\n")
        output = mcflirt(asldata, cost="mutualinfo", out=LOAD)
        asldata = asldata.derived(output["out"].nibImage.get_data(), suffix="_mc")

    if smooth:
        sigma = round(fwhm/2.355, 2)
        log.write("Spatial smoothing with FWHM: %f (sigma=%f)\n" % (fwhm, sigma))
        smoothed = scipy.ndimage.gaussian_filter(asldata.nibImage.get_data(), sigma=sigma)
        asldata = asldata.derived(smoothed, suffix="_smooth")

    log.write("DONE\n\n")
    return asldata

def main():
    """
    Entry point for command line tool
    """
    try:
        p = OptionParser(usage="asl_preproc -i <filename> [options]", version="@VERSION@")
        add_data_options(p, output_type="file")
        add_preproc_options(p)
        options, _ = p.parse_args(sys.argv)

        if not options.asldata:
            sys.stderr.write("Input file not specified\n")
            p.print_help()
            sys.exit(1)

        options.asldata = AslImage(options.asldata, role="Input data", **vars(options))
        if options.output is None:
            options.output = options.asldata.name + "_out"
        options.asldata.summary()
        print("")

        data_preproc = preprocess(**vars(options))
        data_preproc.save(options.output)

    except Exception as e:
        sys.stderr.write("ERROR: " + str(e) + "\n")
        sys.exit(1)

if __name__ == "__main__":
    main()
