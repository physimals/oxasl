"""
ASL data preprocessing command line tool
"""
import sys

from optparse import OptionParser

from .image import AslImage, add_data_options, AslOptionGroup, AslWorkspace

def add_preproc_options(parser, ignore=()):
    g = AslOptionGroup(parser, "Preprocessing", ignore=ignore)
    g.add_option("--diff", dest="diff", help="Perform tag-control subtraction", action="store_true", default=False)
    g.add_option("--smooth", dest="smooth", help="Spatially smooth data", action="store_true", default=False)
    g.add_option("--fwhm", dest="fwhm", help="FWHM for spatial filter kernel", type="float", default=6)
    g.add_option("--mc", dest="mc", help="Motion correct data", action="store_true", default=False)
    g.add_option("--reorder", dest="reorder", help="Re-order data in specified order")
    parser.add_option_group(g)

def main():
    """
    Entry point for command line tool
    """
    try:
        p = OptionParser(usage="asl_preproc", version="@VERSION@")
        add_data_options(p, output_type="file")
        add_preproc_options(p)
        options, _ = p.parse_args(sys.argv)

        options.asldata = AslImage(options.asldata, role="Input data", **vars(options))
        if options.output is None:
            options.output = options.asldata.iname + "_out"
        options.asldata.summary()
        
        wsp = AslWorkspace(echo=options.debug)
        data_preproc = wsp.preprocess(**vars(options))
        data_preproc.save(options.output)

    except Exception as e:
        sys.stderr.write("ERROR: " + str(e) + "\n")
        sys.exit(1)

if __name__ == "__main__":
    main()
