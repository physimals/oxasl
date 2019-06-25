"""
Python implementation of EPI registration
"""
import os
import sys

import numpy as np

import fsl.wrappers as fsl

from oxasl import struc

def epi_reg(wsp, epi_img, use_fmap=False, **kwargs):
    """
    Do EPI registration
    """
    struc.segment(wsp)
    wmseg = wsp.structural.wm_seg
    bbr_sch = os.path.join(os.environ["FSLDIR"], "etc/flirtsch/bbr.sch")

    if wsp.asl2struc is None:
        # Do pre-alignment in the same was as epi_reg
        wsp.asl2struc = fsl.flirt(epi_img, ref=wsp.structural.brain, omat=fsl.LOAD, dof=6, log=wsp.fsllog)["omat"]

    if not use_fmap:
        wsp.log.write(" - Running BBR\n")
        trans = fsl.flirt(epi_img, ref=wsp.structural.struc, dof=6, cost="bbr", wmseg=wmseg, init=wsp.asl2struc, omat=fsl.LOAD, out=fsl.LOAD, schedule=bbr_sch, log=wsp.fsllog)["omat"]
        out_img = fsl.applywarp(epi_img, ref=wsp.structural.struc, out=fsl.LOAD, premat=trans, interp="spline", log=wsp.fsllog)["out"]
        return {"out.nii.gz" : out_img, "out" : trans}
    else:
        dirmap = {
            "x" : (1, "x"),
            "y" : (2, "y"),
            "z" : (3, "z"),
            "-x" : (-1, "x-"),
            "-y" : (-2, "y-"),
            "-z" : (-3, "z-"),
            "x-" : (-1, "x-"),
            "y-" : (-2, "y-"),
            "z-" : (-3, "z-"),
        }
        pedir, fdir = dirmap.get(wsp.pedir, (None, None))
        if pedir is None:
            raise ValueError("Invalid phase encode direction specified: %s" % wsp.pedir)

        if wsp.nofmapreg:
            wsp.fmap2struc = np.identity(4)
            wsp.fmapmag_struc = wsp.fmapmag
        else:
            # Register fmap to structural image
            wsp.log.write(" - Registering fieldmap to structural\n")
            wsp.fmap2struc = fsl.flirt(wsp.fmapmagbrain, ref=wsp.structural.brain, dof=6, omat=fsl.LOAD)["omat"]
            flirt_result = fsl.flirt(wsp.fmapmag, ref=wsp.structural.struc, dof=6, init=wsp.fmap2struc, omat=fsl.LOAD, out=fsl.LOAD, nosearch=True)
            wsp.fmap2struc = flirt_result["omat"]
            wsp.fmapmag_struc = flirt_result["out"]

        # Unmask the fieldmap (necessary to avoid edge effects)
        wsp.fmap_mask = fsl.fslmaths(wsp.fmapmagbrain).abs().bin().run()
        wsp.fmap_mask = fsl.fslmaths(wsp.fmap).abs().bin().mul(wsp.fmap_mask).run()

        # The direction here should take into account the initial affine (it needs to be the direction in the EPI)
        wsp.fmap_unmasked = fsl.fugue(loadfmap=wsp.fmap, mask=wsp.fmap_mask, unmaskfmap=True, savefmap=fsl.LOAD, unwarpdir=fdir)["out"]

        # The following is a NEW HACK to fix extrapolation when fieldmap is too small
        wsp.fmap_struc_pad0 = fsl.applywarp(wsp.fmap_unmasked, ref=wsp.structural.struc, premat=wsp.fmap2struc, out=fsl.LOAD)["out"]
        wsp.fmap_struc_innermask = fsl.fslmaths(wsp.fmap_struc_pad0).abs().bin().run()
        wsp.fmap_struc_dilated = fsl.fugue(loadfmap=wsp.fmap_struc_pad0, mask=wsp.fmap_struc_innermask, unmaskfmap=True, unwarpdir=fdir, savefmap=fsl.LOAD)["savefmap"]
        wsp.fmap_struc = wsp.fmap_struc_dilated

        # Run bbr with fieldmap
        wsp.log.write("Running BBR with fieldmap\n")
        if not wsp.epi_reg_use_weighting:
            refweight = None

        wsp.reg.epi2struc = fsl.flirt(epi_img, ref=wsp.structural.struc, dof=6, cost="bbr", wmseg=wmseg, init=wsp.reg.asl2struc, omat=fsl.LOAD, schedule=bbr_sch, echospacing=wsp.echospacing, pedir=pedir, fieldmap=wsp.fmap_struc, refweight=refweight)["omat"]

        # Make equivalent warp fields
        wsp.log.write("Making warp fields and applying registration to EPI series\n")
        wsp.reg.struc2epi = np.linalg.inv(wsp.reg.epi2struc)
        fsl.concatxfm(wsp.reg.struc2epi, wsp.fmap2struc, outmat=fsl.LOAD)
        fsl.applywarp(wsp.fmap, ref=epi_img, premat=wsp.fmap2epi, out=fsl.LOAD)
        wsp.fmap2epi_mask = fsl.fslmaths(wsp.fmap2epi).abs().bin().run()
        # ${vout}_fieldmaprads2epi -abs -bin ${vout}_fieldmaprads2epi_mask
        ret = fsl.fugue(loadfmap=wsp.fmap2epi, mask=wsp.fmap2epi_mask, saveshift=fsl.LOAD, unmaskshift=True, dwell=wsp.dwell, unwarpdir=fdir)
        print(ret)
        wsp.fmap2epi_shift = ret["fieldmaprads2epi_shift"]
        ret = fsl.convertwarp(wsp.structural.struc, s=wsp.fmap2epi_shift, postmat=wsp.epi2struc, out=fsl.LOAD, shiftdir=fdir, relout=True)
        print(ret)
        warp = ret["out"]
        ret = fsl.applywarp(epi_img, ref=wsp.structural.struc, out=fsl.LOAD, warp1=warp, interp="spline", rel=True)
        print(ret)
        return ret["out"], warp
        # if gdcwarp:
        #     ${FSLDIR}/bin/convertwarp --ref=${vepi} --warp1=${vout}_warp --midmat=${vout}_inv.mat --warp2=${gdcwarp}  --out=${vepi}2${vepi}_undistorted_warp --relout
        #     ${FSLDIR}/bin/applywarp -i ${vepi} -r ${vepi} -o ${vepi}_fully_undistorted -w ${vepi}2${vepi}_undistorted_warp --interp=spline --rel
        #     ${FSLDIR}/bin/flirt -ref ${vrefhead} -in ${vepi}_fully_undistorted -dof 6 -cost bbr -wmseg ${vout}_fast_wmseg -init ${vout}.mat -omat ${vout}.mat -out ${vout} -schedule ${FSLDIR}/etc/flirtsch/bbr.sch
        #     ${FSLDIR}/bin/convertwarp --ref=${vrefhead} --warp1=${vout}_warp --midmat=${vout}_inv.mat --warp2=${gdcwarp} --postmat=${vout}.mat --out=${vout}_warp --relout
        #     ${FSLDIR}/bin/applywarp -i ${vepi} -r ${vrefhead} -o ${vout} -w ${vout}_warp --interp=spline --rel

def main():
    """ Simple test case """
    from fsl.data.image import Image
    from oxasl import Workspace
    wsp = Workspace(savedir="epi_py_out", debug=True)
    wsp.sub("reg")
    wsp.regfrom = Image(sys.argv[1])
    wsp.fslanat = sys.argv[2]
    struc.init(wsp)

    result = epi_reg(wsp, wsp.regfrom)
    wsp.imgreg, wsp.img2struc = result["out.nii.gz"], result["out"]

if __name__ == "__main__":
    main()
