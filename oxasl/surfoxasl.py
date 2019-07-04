import os.path as op 

import toblerone

import workspace
import oxford_asl
from fsl.data.image import Image

def prepare_surf_pvs(wsp):

    # Pipeline: 
    # Do an initial run of oxford_asl using just ASL data to get a perfusion image
    # Register (epi_reg, BBR) the structural to this perfusion image 
    # Run Toblerone in the space of the perfusion image
    # Motion correct the ASL data to the perfusion image 
    # Run oxford_asl with PVEc from surface estimates
    # Optional: run oxford_asl with PVEc from FAST estimates

    # Estimate WM and GM PVs
    struct2asl = wsp.reg.struc2asl
    ref = wsp.reg.regfrom.dataSource
    struct = wsp.reg.regto.dataSource
    pvdir = wsp.sub('surf_pvs')

    if not toblerone.utils.check_anat_dir(wsp.fslanat):
        raise RuntimeError("fsl_anat dir not complete with surfaces")

    if True: 
        pvs, _ = toblerone.estimate_all(ref=ref, anat=wsp.fslanat, struct2ref=struct2asl, flirt=True)

        spc = toblerone.classes.ImageSpace(ref)
        for k,v in pvs.items():
            # setattr(pvdir, k, v)
            spc.saveImage(v, op.join(pvdir.savedir, k+'.nii.gz'))
        
    wm,gm = [ Image(op.join(pvdir.savedir, 'all_%s.nii.gz' % t)) 
        for t in ['WM', 'GM'] ]
    wsp.basil_options.update({"pwm" : wm, "pgm" : gm})

    return 
        