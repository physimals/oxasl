"""
OXASL - Module to generate a suitable mask for ASL data

Copyright (c) 2008-2020 Univerisity of Oxford
"""
import numpy as np
import scipy as sp

from fsl.data.image import Image

from oxasl import reg
from oxasl.reporting import LightboxImage

def generate_mask(wsp):
    """
    For compatibility
    """
    run(wsp)

def run(wsp):
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
     - ``aslref`` : ASL registration source image
    """
    if wsp.rois is not None and wsp.rois.mask is not None:
        return

    wsp.sub("rois")
    wsp.log.write("\nGenerating ASL data mask\n")

    # Reporting
    page = wsp.report.page("mask")
    page.heading("Mask generation", level=0)

    if wsp.mask is not None:
        wsp.rois.mask_src = "user"
        mask_source = "provided by user (assumed to be ASL space): %s" % wsp.mask.name
        wsp.rois.mask = wsp.mask
    elif wsp.structural is not None and wsp.structural.struc is not None:
        # Preferred option is to use brain extracted structural
        wsp.rois.mask_src = "struc"
        page.heading("Brain extracted structural image", level=1)
        page.image("struc_brain", LightboxImage(wsp.structural.brain, bgimage=wsp.structural.struc))
        wsp.rois.mask_struc = wsp.structural.brain_mask
        wsp.rois.mask_asl = reg.change_space(wsp, wsp.structural.brain_mask, "asl")
        wsp.rois.mask = Image(sp.ndimage.morphology.binary_fill_holes((wsp.rois.mask_asl.data > 0.25)).astype(np.int), header=wsp.rois.mask_asl.header)
        mask_source = "generated from brain extracting structural image and registering to ASL space"
    else:
        # Alternatively, use registration image (which will be BETed calibration or mean ASL image)
        wsp.rois.mask_src = "aslref"
        wsp.rois.mask = Image((wsp.reg.aslref.data != 0).astype(np.int), header=wsp.reg.aslref.header)
        mask_source = "generated from brain extracted registration ASL image"

    wsp.log.write(" - Mask %s\n" % mask_source)

    page.heading("Masked ASL brain image", level=1)
    page.text("Mask was %s" % mask_source)
    page.text("PW ASL image masked by ASL-space mask")

    if wsp.asldata.iaf in ("diff", "tc", "ct"):
        page.image("mask_outline", LightboxImage(wsp.rois.mask, bgimage=wsp.asldata.perf_weighted(), outline=True))
    else:
        page.image("mask_outline", LightboxImage(wsp.rois.mask, bgimage=wsp.asldata.mean(), outline=True))
