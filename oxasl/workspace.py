"""
Workspace for executing FSL commands
"""

import os
import sys
import errno
import tempfile

import numpy as np

from fsl.data.image import Image
import fsl.wrappers as fsl

from .reporting import RstContent

class Workspace(object):
    """
    A workspace for data processing

    The contents of a workspace are modelled as attributes on the workspace
    object. An initial set of contents can be specified using keyword arguments.
    For command line tools, typically these are provided directly from the
    OptionParser results, e.g.

        options, args = parser.parse_args(sys.argv)
        wsp = Workspace(**vars(options))

    A workspace may optionally be associated with a physical directory. In this
    case certain types of objects will automatically be saved to the workspace.
    Supported types are currently:

         - ``fsl.data.image.Image`` - Saved as Nifti
         - 2D Numpy array - Saved as ASCII matrix

    To avoid saving a particular item, use the ``add`` method rather than
    directly setting an attribute, as it supports a ``save`` option.
    """

    def __init__(self, savedir=None, **kwargs):
        """
        Create workspace

        :param savedir: If specified, use this path to save data. Will be created
                        if it does not not already exist
        :param log:     File stream to write log output to (default: sys.stdout)
        """
        if savedir is not None:
            self._savedir = os.path.abspath(savedir)
            mkdir(savedir, log=kwargs.get("log", sys.stdout))
        else:
            self._savedir = None

        # Defaults might be overridden by kwargs
        self.log = sys.stdout
        self.debug = False
        self.fsllog = None
        for key, value in kwargs.items():
            setattr(self, key, value)
        
        # Default log configuration for FSL wrapper commands
        if not self.fsllog:
            self.fsllog = {"stderr" : self.log}
            if self.debug:
                self.fsllog.update({"stdout" : self.log, "cmd" : self.log})

    def __getattr__(self, name):
        return None

    def __setattr__(self, name, value):
        self.set_item(name, value)

    def set_item(self, name, value, save=True, save_name=None):
        """
        Add an item to the workspace

        The item will be set as an attribute on the workspace
        If a save directory is configured and the value is a supported
        type it will be saved there. This can be disabled by setting save=False

        :param name: Name, must be a valid Python identifier
        :param value: Value to set
        :param save: If False, do not save item even if savedir defined
        :param save: If specified, alternative name to use for saving this item
        """
        if value is not None and save and self._savedir:
            if not save_name:
                save_name = name
            if isinstance(value, Image):
                # Save as Nifti file
                value.save(os.path.join(self._savedir, save_name))
                value.name = save_name
            elif isinstance(value, np.ndarray) and value.ndim == 2:
                # Save as ASCII matrix
                with open(os.path.join(self._savedir, save_name), "w") as tfile:
                    tfile.write(matrix_to_text(value))
        super(Workspace, self).__setattr__(name, value)

    def sub(self, name, **kwargs):
        """
        Create a sub-workspace, (i.e. a subdir of this workspace)

        This inherits the log configuration from the parent workspace. If
        a savedir is defined, it will be created with a savedir which is
        a subdirectory of the original workspace. Additional data may be 
        set using keyword arguments.
        
        :param name: Name of subdir
        """
        if self._savedir:
            savedir = os.path.join(self._savedir, name)
        else:
            savedir = None
        setattr(self, name, Workspace(savedir=savedir, log=self.log, debug=self.debug, **kwargs))
        return getattr(self, name)

class AslWorkspace(Workspace):
    """
    Workspace with methods used by the ``oxasl`` tool

    Known attributes are as follows. As a bare minimum, ``asldata`` must be
    set before calling any methods

    asldata : AslImage containing main ASL data

    diffasl : Tag-control subtracted AslImage
    diffasl_mean : Mean across repeats at each TI/PLD of diffasl
    
    meanasl : 3D Image containing the timeseries mean of asldata
    meanasl_brain : Brain extracted copy of meanasl
    meanasl_brain_mask : Brain mask from meanasl

    pwi : 3D Perfusion weighted image (mean of mean across repeats)
    pwi_brain: Brain extracted copy of pwi
    pwi_brain_mask : Brain mask from pwi

    struc : Structural image
    struc_brain : Brain extracted structural image
    struc_brain_mask : Brain mask from structural image

    calib : Calibration image - motion corrected and averaged to 3D
    calib_brain : Brain-extracted calibration image
    calib_brain_mask : Brain mask from calibration image

    cref : Calibration reference image - motion corrected and averaged to 3D
    cblip : cblip image  - motion corrected and averaged to 3D

    asl2struc : ASL -> structural rigid body transformation matrix
    struc2asl : Structural -> ASL rigid body transformation matrix

    regfrom : 3D brain-extracted image in ASL space used as ASL-> structural registration target
    """

    def preproc_asl(self):
        """
        Preprocessing on the main ASL data - calculate averaged images and run brain extraction
        """
        if self.asldata is not None and self.meanasl is None:
            self.log.write("\nPre-processing ASL data: %s\n" % self.asldata.name)
            self.diffasl = self.asldata.diff()
            self.diffasl_mean = self.diffasl.mean_across_repeats()
            self.meanasl = self.asldata.mean()
            bet_result = fsl.bet(self.meanasl, seg=True, mask=True, fracintensity=0.2, output=fsl.LOAD, log=self.fsllog)
            self.meanasl_brain = bet_result["output"]
            self.meanasl_brain_mask = bet_result["output_mask"]
            self.pwi = self.asldata.perf_weighted()
            bet_result = fsl.bet(self.pwi, seg=True, mask=True, fracintensity=0.2, output=fsl.LOAD, log=self.fsllog)
            self.pwi_brain = bet_result["output"]
            self.pwi_brain_mask = bet_result["output_mask"]
            self.log.write(" - Generated subtracted and averaged copies of input data\n")

    def preproc_calib(self):
        """
        Preprocess calibration images
        """
        for img_type in ("calib", "cblip", "cref"):
            if hasattr(self, img_type) and not hasattr(self, img_type + "_brain"):
                img = getattr(self, img_type)
                self.log.write("\nPre-processing calibration image: %s (%s)\n" % (img.name, img_type))
                data = img.data
                if img.ndim == 4:
                    if img.shape[3] > 1:
                        self.log.write(" - Removing first volume to ensure data is in steady state\n")
                        data = data[..., :-1]
                    
                if data.shape[3] > 1:
                    if self.moco:
                        self.log.write(" - Motion correcting\n")
                        data = fsl.mcflirt(data, out=fsl.LOAD, log=self.fsllog)
            
                    self.log.write(" - Taking mean across time axis\n")
                    data = np.mean(data, axis=-1)

                img = Image(data, name=img_type, header=img.header)
                setattr(self, img_type, img)
                if img_type == "calib":
                    self.log.write(" - Doing brain extraction\n")
                    bet_result = fsl.bet(img, seg=True, mask=True, output=fsl.LOAD, log=self.fsllog)
                    setattr(self, img_type + "_brain", bet_result["output"])
                    setattr(self, img_type + "_brain_mask", bet_result["output_mask"])

    def preproc_struc(self):
        """
        Do preprocessing on supplied structural data - copy relevant image and do brain extraction
        """
        if self.fslanat:
            self.log.write("\nUsing FSL_ANAT output directory for structural data: %s\n" % self.fslanat)
            biascorr = os.path.join(self.fslanat, "T1_biascorr")
            biascorr_brain = os.path.join(self.fslanat, "T1_biascorr_brain")
            if os.path.isfile(biascorr) and os.path.isfile(biascorr_brain):
                self.log.write(" - Using bias-corrected structural images")
                self.struc = fsl.Image(biascorr)
                self.struc_brain = fsl.Image(biascorr_brain)
            else:
                self.log.write(" - Using non bias-corrected structural images")
                self.struc = fsl.Image(os.path.join(self.fslanat, "T1"))
                self.struc_brain = fsl.Image(os.path.join(self.fslanat, "T1_brain"))
                
            warp = os.path.join(self.fslanat, "T1_to_MNI_nonlin_coeff")
            mat = os.path.join(self.fslanat, "T1_to_MNI_lin.mat")
            if os.path.isfile(warp):
                self.struc2std_warp = warp
            elif os.path.isfile(mat):
                self.struc2std_mat = mat

        elif self.struc:
            self.log.write("\nUsing structural image provided: %s\n" % self.struc.name)
        #elif self.struc_lores
        #    self.log.write("Low-resolution tructural image: %s\n" % self.struc_lores.name)
        else:
            self.log.write("\nWARNING: No structural data supplied - output will be native space only\n")

        if self.struc is not None and self.struc_brain is None:
            self.log.write(" - Brain-extracting structural image\n")
            bet_result = fsl.bet(self.struc, output=fsl.LOAD, seg=True, mask=True, log=self.fsllog)
            self.struc_brain = bet_result["output"]
            self.struc_brain_mask = bet_result["output_mask"]
        elif self.struc_brain is not None and self.struc_brain_mask is None:
            self.struc_brain_mask = fsl.fslmaths(self.struc_brain).bin().run()

    def get_regfrom(self):
        """
        Set the 3D image to be used as the ASL registration target for structural->ASL registration
        """
        self.log.write("\nGetting image to use for ASL->structural registration)\n")
        self.preproc_calib()
        self.preproc_asl()
        if self.regfrom is None:
            if self.calib_brain is not None:
                self.log.write(" - Registration source is calibration image (brain extracted)\n")
                self.regfrom = self.calib_brain
            elif self.reg_from_pwi and self.pwi_brain:
                self.log.write(" - Registration source is perfusion weighted image (brain extracted)\n")
                self.regfrom = self.pwi_brain
            elif self.meanasl_brain:
                self.log.write(" - Registration source is mean ASL image (brain extracted)\n")
                self.regfrom = self.meanasl_brain

    def reg_asl2struc(self, **kwargs):
        """
        Registration of ASL images to structural image
        """
        import oxasl.reg as reg
        kwargs.update({"log" : self.log, "fsllog" : self.fsllog})
        if self.asl2struc is None:
            self.struc2asl = None
            self.preproc_struc()
            self.log.write("\nRegistering ASL data to structural data\n")
            if kwargs.get("do_flirt", True):
                self.regto, self.asl2struc = reg.reg_flirt(self.regfrom, self.struc_brain, **kwargs)
            if kwargs.get("do_bbr", False):
                self.regto, self.asl2struc = reg.reg_bbr(self.regfrom, self.struc_brain, initmat=self.asl2struc, **kwargs)
        
        if self.asl2struc is not None and self.struc2asl is None:
            self.struc2asl = np.linalg.inv(self.asl2struc)

        if self.asl2struc is not None:
            self.log.write(" - ASL->Structural transform\n")
            self.log.write(str(self.asl2struc) + "\n")
        if self.struc2asl is not None:
            self.log.write(" - Structural->ASL transform\n")
            self.log.write(str(self.struc2asl) + "\n")

    def segment_struc(self):
        """
        Segment the structural image
        """
        if None in (self.wm_seg, self.gm_seg, self.csf_seg):
            self.preproc_struc()
            if self.fslanat:
                self.log.write("\nGetting structural segmentation from FSL_ANAT output\n")
                self.csf_pv_struc = fsl.Image(os.path.join(self.fslanat, "T1_fast_pve_0"))
                self.gm_pv_struc = fsl.Image(os.path.join(self.fslanat, "T1_fast_pve_1"))
                self.wm_pv_struc = fsl.Image(os.path.join(self.fslanat, "T1_fast_pve_2"))
            
                try:
                    self.bias_struc = fsl.Image(os.path.join(self.fslanat, "T1_fast_bias"))
                    self.log.write(" - Bias field extracted sucessfully")
                except:
                    self.log.write(" - No bias field found")
            elif self.fastdir:
                raise NotImplementedError("Specifying FAST output directory")
            elif self.struc:
                self.log.write("\nRunning FAST to segment structural image\n")
                fast_result = fsl.fast(self.struc_brain, out=fsl.LOAD, log=self.fsllog)
                self.csf_pv_struc = fast_result["csf_pve"]
                self.gm_pv_struc = fast_result["gm_pve"]
                self.wm_pv_struc = fast_result["wm_pve"]
                self.bias_struc = fast_result["fast_bias"]
            else:
                raise ValueError("\nNo structural data provided - cannot segment\n")

            self.csf_seg_struc = self.csf_pv_struc.data > 0.5
            self.gm_seg_struc = self.gm_pv_struc.data > 0.5
            self.wm_seg_struc = self.wm_pv_struc.data > 0.5

    def generate_mask(self, **kwargs):
        """
        Generate mask for ASL data

        Apart from meanasl_img, all arguments are optional. The mask will be generated in
        different ways depending on what information is provided:

        - If a ready-made mask image is provided, this is returned
        - If a structural image and a transformation to ASL space is provided this will be used.
        Brain extraction will be performed if required
        - If a calibration image is provided, this is used. It is assumed to be in the same space
        as the ASL data
        - If none of the above are present, the ASL data itself is brain extracted to produce the mask

        :param meanasl_img: Averaged ASL image. Could be mean over timeseries or PWI
        :param mask_img: Ready-prepared mask in ASL space
        :param struc_img: Structural image (whole head)
        :param struc_brain_img: Structural image (brain extracted)
        :param struc2asl: Structural->ASL transformation
        :param calib_img: Calibration image registered to ASL data
        """
        self.get_regfrom()
        self.log.write("\nGenerating ASL data mask\n")
        if self.mask is not None:
            mask_source = "provided by user: %s" % self.mask.name
        elif self.struc is not None:
            # Preferred option is to use brain extracted structural
            self.preproc_struc()
            self.reg_asl2struc(do_flirt=True, do_bbr=False)
            brain_mask_asl = fsl.applyxfm(self.struc_brain_mask, self.regfrom, self.struc2asl, out=fsl.LOAD, interp="trilinear", log=self.fsllog)["out"]
            self.mask = fsl.fslmaths(brain_mask_asl).thr(0.25).bin().fillh().run()
            #fslcpgeom(regfrom_img, mask) FIXME
            mask_source = "generated from structural image: %s" % self.struc.name
        else:
            # Alternatively, use registration image (which will be BETed calibration or mean ASL image)
            self.mask = fsl.fslmaths(self.regfrom).bin().run()
            mask_source = "generated from registration ASL image"
        
        self.log.write(" - Mask %s\n" % mask_source)
        
        if self.report:
            mask_report = RstContent()
            mask_report.heading("Mask generation", level=0)
            mask_report.text("Mask source: %s" % mask_source)
            mask_report.heading("Masked brain image", level=1)
            mask_report.image("mask.png")
            self.report.add_rst("mask", mask_report)

            brain_img = np.copy(self.meanasl.data)
            brain_img[self.mask.data == 0] = 0
            self.report.add_lightbox_img("mask.png", Image(brain_img, header=self.meanasl.header))

    def motion_correct(self):
        """
        Motion Correction of ASL data
        
        Note motion correction of multi-volume calibration data is done in preprocessing
        """
        self.log.write("\nMotion Correction\n")
        if self.calib:
            # use supplied image as our reference for motion correction
            # Normally the calibration image since this will be most consistent if the data has a range of different TIs and background suppression etc
            # this also removes motion effects between asldata and calibration image
            self.log.write(" - Using calibration image as reference\n")
            mcflirt_result = fsl.mcflirt(self.asldata, ref=self.calib, mats=True, log=self.fsllog)
            self.asldata_mc = mcflirt_result["out"]
            self.moco_asl = mcflirt_result["mats"]

            # To reduce interpolation of the ASL data change the transformations so that we end up in the space of the central volume of asldata
            asl2calib = self.moco_asl[len(self.moco_asl)/2+1]
            calib2asl = np.linalg.inv(asl2calib)
            self.log.write("middle volume->calib:\n%s\n" % str(asl2calib))
            self.log.write("calib->middle volume:\n%s\n" % str(calib2asl))

            # Convert all the volumes to this space
            asl_moco_mats = [np.dot(mat, calib2asl) for mat in self.moco_asl]
            applyxfm_4d(self.asldata, self.asldata_mc, asl_moco_mats, output=fsl.LOAD)

            # Convert all calibration images to align with asldata
            applyxfm(self.calib, self.asldata, xfm=calib2asl)
            if self.cref:
                applyxfm(self.cref, self.asldata, xfm=calib2asl)
            if self.cblip:
                applyxfm(self.cblip, self.asldata, xfm=calib2asl)

        else:
            self.log.write(" - Using ASL data middle volume as reference\n")
            mcflirt_result = fsl.mcflirt(self.asldata, mats=True, log=self.fsllog)
            self.asldata_mc = mcflirt_result["out"]
            self.moco_asl = mcflirt_result["mats"]
            
        #cat $tempdir/asldata.mat/MAT* > $tempdir/asldata.cat # save the motion matrices for distortion correction if reqd

    def set_output_spaces(self):
        """
        Determine the output spaces we should be using and the transformations into them
        """
        if not self.output_spaces:
            self.output_spaces = {
                "native" : {},
            }

            if self.struc2std_warp:
                self.output_spaces["std"] = {"warp" : self.struc2std_warp}
            elif self.struc2std_mat:
                self.output_spaces["std"] = {"mat" : self.struc2std_mat}

            if "std" in self.output_spaces:
                if self.std_brain:
                    self.output_spaces["std"]["brain"] = self.std_brain
                else:
                    self.output_spaces["std"]["brain"] = os.path.join(os.environ["FSLDIR"], "data", "standard", "MNI152_T1_2mm")

                self.log.write("Standard brain is: %s" % self.output_spaces["std"]["brain"])
            else:
                self.log.write("No standard space transform found - output will be in native/structural space only\n")

    def basil(self):
        """
        Run BASIL modelling on ASL data
        """
        self.log.write("\nRunning BASIL Bayesian modelling on ASL data\n")
        # Single or Multi TI setup
        if self.asldata.ntis == 1:
            # Single TI data - don't try to infer arterial component of bolus duration, we don't have enough info
            self.log.write(" - Operating in Single TI mode - no arterial component, fixed bolus duration\n")
            self.artoff = True
            self.fixbolus = True
            self.singleti = True
            
        if self.wp:
            # White paper mode - this overrides defaults, but can be overwritten by command line 
            # specification of individual parameters
            self.log.write(" - Analysis in white paper mode: T1 default=1.65, BAT default=0, voxelwise calibration\n")
            t1_default = 1.65
            bat_default = 0.0
            self.calib_method = "voxel"
        else:
            t1_default = 1.3
            if self.casl:
                bat_default = 1.3
            else:
                bat_default = 0.7

        if not self.t1: self.t1 = t1_default
        if not self.t1b: self.t1b = 1.65
        if not self.bat: self.bat = bat_default
        if not self.bolusdur: self.bolusdur = 1.8
            
        # if we are doing CASL then fix the bolus duration, except where the user has 
        # explicitly told us otherwise
        if not self.fixbolus: self.fixbolus = self.casl
            
        from oxasl import basil
        if not self.basil_output:
            self.sub("basil_output")
        basil_options = dict(vars(self))

        # Initial BASIL run on mean data
        self.log.write(" - Doing initial fit on mean at each TI\n\n")
        init_wsp = self.basil_output.sub("init")
        basil_options["asldata"] = self.diffasl_mean
        basil.basil(wsp=init_wsp, **basil_options)
        final_step = 1
        while 1:
            if not getattr(init_wsp, "step%i" % (final_step+1), None):
                break
            final_step += 1
        basil_options["initmvn"] = getattr(init_wsp, "step%i" % final_step).finalMVN

        # Main run on full ASL data
        self.log.write("\n - Doing main fit on full ASL data\n\n")
        main_wsp = self.basil_output.sub("main")
        basil_options["asldata"] = self.asldata
        basil.basil(wsp=main_wsp, **basil_options)

def matrix_to_text(mat):
    """
    Convert matrix array to text using spaces/newlines as col/row delimiters
    """
    rows = []
    for row in mat:
        rows.append(" ".join([str(v) for v in row]))
    return "\n".join(rows)

def text_to_matrix(text):
    """
    Convert space or comma separated file to matrix
    """
    fvals = []
    ncols = -1
    lines = text.splitlines()
    for line in lines:
        # Discard comments
        line = line.split("#", 1)[0].strip()
        # Split by commas or spaces
        vals = line.replace(",", " ").split()
        # Ignore empty lines
        if not vals: continue
        # Check correct number of columns
        if ncols < 0: ncols = len(vals)
        elif len(vals) != ncols:
            raise ValueError("File must contain a matrix of numbers with fixed size (rows/columns)")
        # Check all data is numeric
        for val in vals:
            try:
                float(val)
            except:
                raise ValueError("Non-numeric value '%s' found in matrix text" % val)
        fvals.append([float(v) for v in vals])     
    return np.array(fvals)

def mkdir(dirname, fail_if_exists=False, warn_if_exists=True, log=sys.stdout):
    """
    Create a directory, including necessary subdirs
    """
    try:
        os.makedirs(dirname)
    except OSError as exc:
        if exc.errno == errno.EEXIST:
            if fail_if_exists: raise
            elif warn_if_exists: log.write("WARNING: mkdir - Directory %s already exists\n" % dirname)
    return os.path.abspath(dirname)

def tempdir(suffix, debug=False, log=sys.stdout):
    """
    Create a temporary directory

    :param debug: If True, creates directory in current working directory
    """
    if debug:
        tmpdir = os.path.join(os.getcwd(), "tmp_%s" % suffix)
        mkdir(tmpdir, log=log)
    else:
        tmpdir = tempfile.mkdtemp("_%s" % suffix)
    return tmpdir
