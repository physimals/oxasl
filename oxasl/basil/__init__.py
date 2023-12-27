"""
OXASL - BASIL: Quantification using Bayesian model fitting

The BASIL module is a little more complex than the other Workspace based
modules because of the number of options available and the need for flexibility
in how the modelling steps are run.

The main function is ``basil`` which performs model fitting on ASL data
in the Workspace ``asldata`` attribute.

    wsp = Workspace()
    wsp.asldata = AslImage("asldata.nii.gz", tis=[1.6,])
    wsp.infertiss = True
    basil.run(wsp.sub("basil"))
    wsp.basil.finalstep.mean_ftiss.save("mean_ftiss.nii.gz")

Because of the number of options possible for the modelling process, the
workspace attribute ``basil_options`` can be set as a dictionary of extra
options relevant only to Basil:

    wsp = Workspace()
    wsp.asldata = AslImage("asldata.nii.gz", tis=[1.6,])
    wsp.basil_options = {"infertiss" : True, "spatial" : True}
    basil.run(wsp.sub("basil"))
    wsp.basil.finalstep.mean_ftiss.save("mean_ftiss.nii.gz")

All options specified in basil_options are either consumed by Basil, or
if not passed directly to the model.

Copyright (c) 2008-2020 Univerisity of Oxford
"""
"""
OXASL - BASIL: Quantification using Bayesian model fitting

Copyright (c) 2008-2020 Univerisity of Oxford
"""
from oxasl.options import OptionCategory, OptionGroup

from . import fabber_method, svb_method, vaby_method

class Options(OptionCategory):
    """
    Options for corrections of the input data
    """

    def __init__(self, **kwargs):
        OptionCategory.__init__(self, "corrections")

    def groups(self, parser):
        ret = []

        g = OptionGroup(parser, "BASIL quantification")
        g.add_option("--basil-mask", help="Masking policy to use for model fitting. Does not affect analysis mask used in rest of pipeline. 'dilate' means dilate the default analysis mask. 'none' means use no masking",
                     type="choice", choices=["default", "dilated", "none"])
        g.add_option("--basil-options", "--fit-options", help="File containing additional options for Basil model fitting", type="optfile", default=None)
        g.add_option("--basil-method", help="Model fitting method", type="choice", choices=["fabber", "vaby", "svb"])
        ret.append(g)

        g = OptionGroup(parser, "BASIL partial volume correction (PVEc)")
        g.add_option("--pvcorr", help="Apply PVEc using FAST estimates taken from --fslanat dir", action="store_true", default=False)
        g.add_option("--pvgm", help="GM PV estimates in ASL space (apply PVEc only, don't estimate PVs)", type="image", default=None)
        g.add_option("--pvwm", help="As above, WM PV estimates in ASL space", type="image", default=None)
        g.add_option("--pvcsf", help="As above, CSF PV estimates in ASL space", type="image", default=None)
        ret.append(g)

        return ret

def run(wsp):
    method = wsp.ifnone("basil_method", "fabber")
    wsp.log.write(" - Running Bayesian model fitting using method '%s'\n" % method)
    if method == "fabber":
        fabber_method.run(wsp)
    elif method == "svb":
        svb_method.run(wsp)
    elif method == "vaby":
        vaby_method.run(wsp)
    else:
        raise ValueError("Unknown BASIL method: %s" % method)
