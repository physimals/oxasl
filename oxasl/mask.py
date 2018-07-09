"""
Functions for generating a suitable mask for ASL data
"""
import sys

from .options import AslOptionParser, GenericOptions
from .image import AslImage, AslImageOptions
from .struc import StructuralImageOptions
from .workspace import AslWorkspace

from .reporting import Report
from ._version import __version__

def main():
    """
    Entry point for command line tool
    """
    try:
        parser = AslOptionParser(usage="asl_mask -i <asl_image> [options...]", version=__version__)
        parser.add_option("-c", "--calib", dest="calib", help="Calibration image", default=None)
        parser.add_option("--use-pwi", help="Use the perfusion weighted average rather than the timeseries mean", action="store_true", default=False)
        parser.add_category(AslImageOptions())
        parser.add_category(StructuralImageOptions())
        parser.add_category(GenericOptions())

        options, _ = parser.parse_args(sys.argv)
        options.mask = None # No point in using command line tool if you already have a mask!
        wsp = AslWorkspace(report=Report("mask_report"), **vars(options))

        if not options.asldata:
            sys.stderr.write("Input file not specified\n")
            parser.print_help()
            sys.exit(1)

        wsp.asldata = AslImage(wsp.asldata, **parser.filter(vars(options), "image"))
        wsp.asldata.summary()

        wsp.generate_mask()
        
        if wsp.output is None:
            wsp.output = wsp.asldata.name + "_mask"
        wsp.mask.save(wsp.output)

    except ValueError as exc:
        sys.stderr.write("ERROR: " + str(exc) + "\n")
        sys.exit(1)

if __name__ == "__main__":
    main()
