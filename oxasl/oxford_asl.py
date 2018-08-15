#!/bin/env python
"""
OXFORD_ASL: Converts ASL images in to perfusion maps

Michael Chappell, FMRIB Image Analysis Group & IBME

Copyright (c) 2008-2013 Univerisity of Oxford
"""

import sys
import os
import traceback

from fsl.data.image import Image

from ._version import __version__
from .image import AslImage, AslImageOptions
from . import calib, struc, basil, preproc, mask, corrections

from .calib import CalibOptions
from .struc import StructuralImageOptions
from .basil import BasilOptions
from .workspace import Workspace
from .options import AslOptionParser, GenericOptions, OptionCategory, IgnorableOptionGroup
from .reg import reg_asl2struc

class OxfordAslOptions(OptionCategory):
    """
    OptionCategory which contains options for preprocessing ASL data
    """

    def __init__(self, **kwargs):
        OptionCategory.__init__(self, "oxford_asl", **kwargs)

    def groups(self, parser):
        ret = []
        g = IgnorableOptionGroup(parser, "Main Options")
        g.add_option("--wp", dest="wp", help="Analysis which conforms to the 'white papers' (Alsop et al 2014)", action="store_true", default=False)
        g.add_option("--mc", dest="mc", help="Motion correct data", action="store_true", default=False)
        ret.append(g)
        g = IgnorableOptionGroup(parser, "Acquisition/Data specific")
        g.add_option("--casl", dest="casl", help="ASL acquisition is  pseudo cASL (pcASL) rather than pASL", action="store_true", default=False)
        g.add_option("--bolus", dest="bolus", help="Bolus duration", type=float, default=1.8)
        g.add_option("--bat", dest="bat", help="Bolus arrival time (default=0.7 (pASL), 1.3 (cASL)", type=float)
        g.add_option("--t1", dest="t1", help="Tissue T1 value", type=float, default=1.3)
        g.add_option("--t1b", dest="t1b", help="Blood T1 value", type=float, default=1.65)
        g.add_option("--slicedt", dest="slicedt", help="Timing difference between slices", type=float, default=0.0)
        g.add_option("--sliceband", dest="sliceband", help="Number of slices per pand in multi-band setup", type=int)
        g.add_option("--artsupp", dest="artsupp", help="Arterial suppression (vascular crushing) was used", action="store_true", default=False)
        ret.append(g)
        return ret

"""
    echo "Partial volume correction"
    echo " --pvcorr    : Do partial volume correction"
    echo "  PV estimates will be taken from:"
    echo "   fsl_anat dir (--fslanat), if supplied"
    echo "   exising fast segmentation (--fastsrc), if supplied"
    echo "   FAST segmenation of structural (if using -s and --sbet)"
    echo "   User supplied PV estimates (--pvgm, --pvwm)"   
    echo "   --pvgm    : Partial volume estimates for GM"
    echo "   --pvwm    : Partial volume estimates for WM"

    echo " Epochs "
    echo " --elen      : Length of each epoch in TIs"
    echo " --eol       : Overlap of each epoch in TIs- {deafult: 0}"
"""

def main():
    debug = True
    try:
        parser = AslOptionParser(usage="oxford_asl -i <asl_image> [options]", version=__version__)
        parser.add_category(AslImageOptions())
        parser.add_category(StructuralImageOptions())
        parser.add_category(OxfordAslOptions())
        parser.add_category(CalibOptions(ignore=["perf", "tis"]))
        parser.add_category(BasilOptions(ignore=["spatial"]))
        parser.add_category(corrections.DistcorrOptions())
        parser.add_category(GenericOptions())

        options, _ = parser.parse_args(sys.argv)
        if not options.output:
            options.output = "oxasl"

        if not options.asldata:
            sys.stderr.write("Input file not specified\n")
            parser.print_help()
            sys.exit(1)

        # Some oxford_asl command-line specific defaults
        options.spatial = True
        if options.calib is not None and options.calib_method is None:
            if options.struc is not None:
                options.calib_method = "refregion"
            else:
                options.calib_method = "voxelwise"

        #basil_options = parser.filter(vars(options), "basil", consume=True)
        wsp = Workspace(savedir=options.output, **vars(options))
        #wsp.basil_options = basil_options
        wsp.asldata = AslImage(options.asldata, **parser.filter(vars(wsp), "image"))

        oxasl(wsp)
        report_build_dir = None
        if wsp.debug:
            report_build_dir = os.path.join(options.output, "report_build")
        wsp.report.generate_html(os.path.join(options.output, "report"), report_build_dir)

    except Exception as e:
        sys.stderr.write("ERROR: " + str(e) + "\n")
        if debug:
            traceback.print_exc()
        sys.exit(1)

def oxasl(wsp):
    """
    Run standard processing on ASL data

    This method requires wsp to be a Workspace containing certain standard information.
    As a minimum, the attribute ``asldata`` must contain an AslImage object.
    """
 
    wsp.log.write("OXFORD_ASL - running\n")
    wsp.log.write("Version: %s\n" % __version__)
    wsp.log.write("Input ASL data: %s\n" % wsp.asldata.name)
    wsp.asldata.summary(wsp.log)

    # Create output workspace. Always output in native space
    wsp.sub("output").sub("native")
    wsp.do_flirt, wsp.do_bbr = True, False # FIXME

    preproc.preproc_asl(wsp)
    calib.preproc_calib(wsp)
    struc.preproc_struc(wsp)

    if wsp.mc: 
        corrections.get_motion_correction(wsp)

    if wsp.fmap:
        corrections.get_fieldmap_correction(wsp)

    corrections.apply_corrections(wsp)
    mask.generate_mask(wsp)
    basil.basil(wsp, output_wsp=wsp.sub("basil"))

    wsp.do_flirt, wsp.do_bbr = False, True # FIXME

    new_regfrom = wsp.basil.main.finalstep.mean_ftiss.data
    new_regfrom[new_regfrom < 0] = 0
    wsp.regfrom = Image(new_regfrom, header=wsp.regfrom.header)
    wsp.asl2struc_initial = wsp.asl2struc
    wsp.struc2asl_initial = wsp.struc2asl
    wsp.done("reg_asl2struc", status=False)
    reg_asl2struc(wsp)

    do_output(wsp)
    if wsp.calib is not None:
        wsp.output.native.perfusion_calib = calib.calibrate(wsp, wsp.output.native.perfusion)

    wsp.log.write("\nOutput is %s\n" % wsp.savedir)
    wsp.log.write("OXFORD_ASL - done\n")

def do_output(wsp):
    for space in dir(wsp.output):
        output_wsp = getattr(wsp.output, space)
        if isinstance(output_wsp, Workspace):
            # Make negative values = 0
            img = wsp.basil.main.finalstep.mean_ftiss
            img.data[img.data < 0] = 0
            output_wsp.perfusion = img

def set_output_spaces(wsp):
    """
    Determine the output spaces we should be using and the transformations into them
    """
    if not wsp.output_spaces:
        

        if wsp.struc2std_warp:
            wsp.output_spaces["std"] = {"warp" : wsp.struc2std_warp}
        elif wsp.struc2std_mat:
            wsp.output_spaces["std"] = {"mat" : wsp.struc2std_mat}

        if "std" in wsp.output_spaces:
            if wsp.std_brain:
                wsp.output_spaces["std"]["brain"] = wsp.std_brain
            else:
                wsp.output_spaces["std"]["brain"] = os.path.join(os.environ["FSLDIR"], "data", "standard", "MNI152_T1_2mm")

            wsp.log.write("Standard brain is: %s" % wsp.output_spaces["std"]["brain"])
        else:
            wsp.log.write("No standard space transform found - output will be in native/structural space only\n")
