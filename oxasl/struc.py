from fsl.wrappers import bet, LOAD

from .wrappers import fast
from .options import OptionCategory, IgnorableOptionGroup

class StructuralImageOptions(OptionCategory):
    """
    OptionGroup which contains options for describing a structural image
    """

    def __init__(self, title="Structural image", **kwargs):
        OptionCategory.__init__(self, "struc", **kwargs)
        self.title = title

    def groups(self, parser):
        group = IgnorableOptionGroup(parser, self.title, ignore=self.ignore)
        group.add_option("-s", "--struc", dest="struc", help="Structural image", type="image", default=None)
        group.add_option("--sbet", "--struc-brain", "--struc-bet", dest="struc_brain", type="image", help="Structural image (brain extracted)", default=None)
        group.add_option("--struc2asl", help="Structural->ASL transformation matrix", default=None)
        group.add_option("--asl2struc", help="ASL->Structural transformation matrix", default=None)
        group.add_option("--wm-seg", help="White matter segmentation of structural image", type="image", default=None)
        group.add_option("--gm-seg", help="Grey matter segmentation of structural image", type="image", default=None)
        group.add_option("--csf-seg", help="CSF segmentation of structural image", type="image", default=None)
        group.add_option("--fslanat", help="FSL_ANAT output directory for structural information", default=None)
        group.add_option("--fastsrc", dest="fastsrc", help="Images from a FAST segmentation - if not set FAST will be run on structural image")
        group.add_option("--senscorr", dest="senscorr", help="Use bias field (from segmentation) for sensitivity correction", action="store_true", default=False)
        
        return [group, ]

def brain_extract(options):
    """
    Brain extract the structural image if required
    """
    if options.struc_brain is None:
        if options.struc is None:
            raise ValueError("No structural data provided")
        bet_result = bet(options.struc, output=LOAD, seg=True, mask=True)
        options.struc_brain = bet_result["output"]
        options.struc_brain_mask = bet_result["output_mask"]

def seg(options):
    """
    Segment the structural image if required
    """
    if None in (options.wm_seg, options.gm_seg, options.csf_seg):
        if options.fslanat:
            raise NotImplementedError("Getting segmentation from FSL_ANAT output")
        else:
            brain_extract(options)
            fast_result = fast(options.struc_brain, out=LOAD)
            if options.wm_seg is None:
                options.wm_seg = fast_result["wm_pve"].data > 0.5
            if options.gm_seg is None:
                options.gm_seg = fast_result["gm_pve"].data > 0.5
            if options.csf_seg is None:
                options.csf_seg = fast_result["csf_pve"].data > 0.5
        