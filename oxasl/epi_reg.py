def epi_reg(wsp, epi_img, use_fmap):
    """ 
    Do EPI registration
    """
    struc.segment(wsp)
    wmseg = wsp.wm_seg_struc
    
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

    bbr_sch = os.path.join(os.environ["FSLDIR"], "etc/flirtsch/bbr.sch")

    # make a WM edge map for visualisation (good to overlay in FSLView)
    fsl.maths(wmseg).edge().bin().mask(wmseg).run(out=fsl.LOAD)

    if not use_fmap:
        # NO FIELDMAP
        wsp.log.write("Running BBR\n")
        result = fsl.flirt(ref=wsp.struc, in=epi_img, dof=6, cost="bbr", wmseg=wmseg, init=wsp.asl2struc. omat=fsl.LOAD, out=fsl.LOAD, schedule=bbr_sch)
        fsl.applywarp(epi_img, ref=wsp.struc, out=fsl.LOAD, premat=result["omat"], interp="spline")
    else:
        # WITH FIELDMAP
        if not wsp.nofmapreg:
            # Register fmap to structural image
            echo "Registering fieldmap to structural"
            wsp.fmap2struc = fsl.flirt(in=wsp.fmapmagbrain, ref=wsp.struc_brain, dof=6, omat=fsl.LOAD)["omat"]
            flirt_result = fsl.flirt(in=wsp.fmapmag, ref=wsp.struc, dof=6, init=wsp.fmap2struc, omat=fsl.LOAD, out=LOAD, nosearch=True)
            wsp.fmap2struc = flirt_result["omat"]
            wsp.fmapmag_struc = flirt_result["out"]
        else:
            wsp.fmap2struc = np.identity(4)
            wsp.fmapmag_struc = wsp.fmapmag
        
        # Unmask the fieldmap (necessary to avoid edge effects)
        wsp.fmap_mask = fsl.fslmaths(wsp.fmapmagbrain).abs().bin().run(out=LOAD)
        wsp.fmap_mask = fsl.fslmaths(wsp.fmap).abs().bin().mul(wsp.fmap_mask).run(out=fsl.LOAD)
        # The direction here should take into account the initial affine (it needs to be the direction in the EPI)
        fsl.fugue(loadfmap=wsp.fmap, mask=fmap_mask, unmaskfmap=True, savefmap=fsl.LOAD, unwarpdir=fdir)
        
        # The following is a NEW HACK to fix extrapolation when fieldmap is too small
        wsp.fmap_struc_pad0 = fsl.applywarp(wsp.fmap_unmasked, ref=wsp.struc, premat=wsp.fmap2struc, out=fsl.LOAD)["out"]
        wsp.fmap_struc_innermask = fsl.fslmaths(wsp.fmap_struc_pad0).abs().bin().run(out=fsl.LOAD)
        wsp.fmap_struc_dilated = fsl.fugue(loadfmap=wsp.fmap_struc_pad0, mask=fmap_struc_innermask, unmaskfmap=True, unwarpdir=fdir, savefmap=fsl.LOAD)["savefmap"]
        wsp.fmap_struc = wsp.fmap_struc_dilated

        # Run bbr with fieldmap
        wsp.log.write("Running BBR with fieldmap\n")
        if not use_weighting:
            refweight = None
            
        wsp.epi2struc = fsl.flirt(ref=wsp.struc, in=epi_img, dof=6, cost="bbr", wmseg=wmseg, init=wsp.asl2struc, omat=fsl.LOAD, schedule=bbr_sch echospacing=wsp.echospacing, pedir=pedir, fieldmap=wsp.fmap_struc, refweight=refweight)["omat"]
        
        # Make equivalent warp fields
        wsp.log.write("Making warp fields and applying registration to EPI series\n")
        wsp.struc2epi = np.linalg.inv(wsp.epi2struc)
        fsl.convert_xfm(omat=fsl.LOAD, concat=True, wsp.struc2epi, fmap2struc)
        fsl.applywarp(wsp.fmap, ref=epi_img, premat=wsp.fmap2epi, out=fsl.LOAD)
        wsp.fmap2epi_mask = fsl.fslmaths(wsp.fmap2epi).abs().bin().run(out=fsl.LOAD) ${vout}_fieldmaprads2epi -abs -bin ${vout}_fieldmaprads2epi_mask
        fsl.fugue(loadfmap=wsp.fmap2epi, mask=wsp.fmap2epi_mask, saveshift=fsl.LOAD, unmaskshift=True, dwell=wsp.dwell, unwarpdir=fdir)
        fsl.convertwarp(wsp.struc, s=wsp.fmap2epi_shift, postmat=..., out=fsl.LOAD, shiftdir=fdir, relout=True)
        fsl.applywarp(epi_img, ref=wsp.struc, out=fsl.LOAD, warp1=_warp, interp="spline", rel=True)
        
        # if gdcwarp:
        #     ${FSLDIR}/bin/convertwarp --ref=${vepi} --warp1=${vout}_warp --midmat=${vout}_inv.mat --warp2=${gdcwarp}  --out=${vepi}2${vepi}_undistorted_warp --relout
        #     ${FSLDIR}/bin/applywarp -i ${vepi} -r ${vepi} -o ${vepi}_fully_undistorted -w ${vepi}2${vepi}_undistorted_warp --interp=spline --rel
        #     ${FSLDIR}/bin/flirt -ref ${vrefhead} -in ${vepi}_fully_undistorted -dof 6 -cost bbr -wmseg ${vout}_fast_wmseg -init ${vout}.mat -omat ${vout}.mat -out ${vout} -schedule ${FSLDIR}/etc/flirtsch/bbr.sch
        #     ${FSLDIR}/bin/convertwarp --ref=${vrefhead} --warp1=${vout}_warp --midmat=${vout}_inv.mat --warp2=${gdcwarp} --postmat=${vout}.mat --out=${vout}_warp --relout
        #     ${FSLDIR}/bin/applywarp -i ${vepi} -r ${vrefhead} -o ${vout} -w ${vout}_warp --interp=spline --rel
        


