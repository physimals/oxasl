"""
Functions for generating a suitable mask for ASL data
"""
import sys

import numpy as np

import fsl.wrappers as fsl
from fsl.data.image import Image

from oxasl import __version__, AslImage, Workspace, image, reg, struc
from oxasl.options import AslOptionParser, GenericOptions
from oxasl.reporting import LightboxImage

def generate_mask(wsp):
    """
    Generate mask for ASL data

    - If a ready-made mask image is provided or has already been generated, this is returned
    - If a structural image is provided this will be used. Brain extraction and registration 
      will be performed if required
    - If a calibration image is provided, this is used. It is assumed to be in the same space
      as the ASL data
    - If none of the above are present, the ASL data itself is averaged and brain extracted 
      to produce the mask

    Required workspace attributes
    -----------------------------

    Formally there are no required attributes, however at least one image must be provided
    which enables a mask to be generated.

    Optional workspace attributes
    -----------------------------

     - ``asldata`` : ASL data image
     - ``mask``    : Existing brain mask
     - ``struc``   : Structural image (wholehead)
     - ``struc_brain``: Already brain-extracted structural image
     - ``asl2struc`` : Existring ASL->Structural space transformation matrix
     - ``calib``   : Calibration image
     - ``regfrom`` : ASL registration source image
    """
    if wsp.rois is not None and wsp.rois.mask is not None:
        return
        
    wsp.sub("rois")

    # Reporting
    page = wsp.report.page("mask")
    page.heading("Mask generation", level=0)

    if wsp.mask is not None:
        wsp.rois.mask_src = "user"
        mask_source = "provided by user (assumed to be ASL space): %s" % wsp.mask.name
        wsp.rois.mask = wsp.mask
    elif wsp.structural.struc is not None:
        # Preferred option is to use brain extracted structural
        wsp.rois.mask_src = "struc"
        struc.init(wsp)
        page.heading("Brain extracted structural image", level=1)
        page.image("struc_brain", LightboxImage(wsp.structural.brain, bgimage=wsp.structural.struc))
        brain_mask_asl = reg.struc2asl(wsp, wsp.structural.brain_mask)
        wsp.rois.mask = fsl.fslmaths(brain_mask_asl).thr(0.25).bin().fillh().run()
        mask_source = "generated from brain extracting structural image and registering to ASL space"
    else:
        # Alternatively, use registration image (which will be BETed calibration or mean ASL image)
        reg.get_regfrom(wsp)
        wsp.rois.mask_src = "regfrom"
        wsp.rois.mask = Image((wsp.reg.regfrom.data != 0).astype(np.int), header=wsp.reg.regfrom.header)
        mask_source = "generated from brain extracted registration ASL image"
    
    wsp.log.write("\nGenerated ASL data mask\n")
    wsp.log.write(" - Mask %s\n" % mask_source)
    
    page.heading("Masked ASL brain image", level=1)
    page.text("Mask was %s" % mask_source)
    page.text("PW ASL image masked by ASL-space mask")

    if wsp.asldata.iaf in ("diff", "tc", "ct"):
        page.image("mask_outline", LightboxImage(wsp.rois.mask, bgimage=wsp.asldata.perf_weighted(), outline=True))
    else:
        page.image("mask_outline", LightboxImage(wsp.rois.mask, bgimage=wsp.asldata.mean(), outline=True))

def main():
    """
    Entry point for command line tool
    """
    try:
        parser = AslOptionParser(usage="asl_mask -i <asl_image> [options...]", version=__version__)
        parser.add_option("--calib", "-c", help="Calibration image", default=None)
        parser.add_option("--use-pwi", help="Use the perfusion weighted average rather than the timeseries mean", action="store_true", default=False)
        parser.add_category(image.AslImageOptions())
        parser.add_category(struc.StructuralImageOptions())
        parser.add_category(GenericOptions())

        options, _ = parser.parse_args(sys.argv)
        options.mask = None # No point in using command line tool if you already have a mask!
        wsp = Workspace(**vars(options))

        if not options.asldata:
            sys.stderr.write("Input file not specified\n")
            parser.print_help()
            sys.exit(1)

        wsp.asldata = AslImage(wsp.asldata, **parser.filter(vars(options), "image"))
        wsp.asldata.summary()

        wsp.generate_mask()
        
        if wsp.output is None:
            wsp.output = wsp.asldata.name + "_mask"
        wsp.rois.mask.save(wsp.output)

    except ValueError as exc:
        sys.stderr.write("ERROR: " + str(exc) + "\n")
        sys.exit(1)

if __name__ == "__main__":
    main()
