"""
oxasl.gui.analysis_tab.py

Copyright (c) 2019 University of Oxford
"""
import os

from oxasl.gui.widgets import TabPage, OptionError

class AslAnalysis(TabPage):
    """
    Tab page containing model fitting and analysis options
    """

    def __init__(self, parent, idx, n):
        TabPage.__init__(self, parent, "Analysis", idx, n)

        self.distcorr_choices = ["Fieldmap", "Calibration image"]

        self.section("Basic analysis options")

        self.outdir_picker = self.file_picker("Output Directory", pick_dir=True)
        self.mask_picker = self.file_picker("Brain Mask", optional=True)
        self.wp_cb = self.checkbox("Analysis which conforms to 'White Paper' (Alsop et al 2014)",
                                   handler=self.wp_changed)

        self.section("Initial parameter values")

        self.bat_num = self.number("Arterial Transit Time (s)", minval=0, maxval=2.5, initial=1.3)
        self.t1_num = self.number("T1 (s)", minval=0, maxval=3, initial=1.3)
        self.t1b_num = self.number("T1b (s)", minval=0, maxval=3, initial=1.65)
        self.ie_num = self.number("Inversion Efficiency", minval=0, maxval=1, initial=0.85)

        self.section("Analysis Options")

        self.spatial_cb = self.checkbox("Adaptive spatial regularization on perfusion", initial=True)
        self.infer_t1_cb = self.checkbox("Incorporate T1 value uncertainty")
        self.macro_cb = self.checkbox("Include macro vascular component")
        self.fixbolus_cb = self.checkbox("Fix label duration", initial=True)

        self.pv_cb = self.checkbox("Partial Volume Correction")
        self.mc_cb = self.checkbox("Motion Correction")

        self.sizer.AddGrowableCol(1, 1)
        #sizer.AddGrowableRow(5, 1)
        self.SetSizer(self.sizer)
        self.next_prev()

    def options(self):
        outdir = self.outdir_picker.GetPath()
        if not outdir:
            raise OptionError("Output directory must be specified")
        elif os.path.isfile(outdir):
            raise OptionError("Output directory already exists and is a file")

        options = {
            "savedir" : outdir,
            "wp" : self.white_paper(),
            "bat" : self.bat_num.GetValue(),
            "t1" : self.t1_num.GetValue(),
            "t1b" : self.t1b_num.GetValue(),
            "alpha" : self.ie_num.GetValue(),
            "spatial" : self.spatial_cb.IsChecked(),
            "infert1" : self.infer_t1_cb.IsChecked(),
            "inferart" : self.macro_cb.IsChecked(),
            "fixbolus" : self.fixbolus_cb.IsChecked(),
            "pvcorr" : self.pv_cb.IsChecked(),
            "mc" : self.mc_cb.IsChecked(),
        }

        if self.mask_picker.checkbox.IsChecked():
            options["mask"] = self.mask_picker.GetPath()

        return options

    def white_paper(self):
        """
        :return: True if white paper mode is enabled
        """
        return self.wp_cb.IsChecked()

    def update(self):
        self.mask_picker.Enable(self.mask_picker.checkbox.IsChecked())
        self.t1_num.Enable(not self.white_paper())
        self.bat_num.Enable(not self.white_paper())
        TabPage.update(self)

    def wp_changed(self, _):
        """
        Update values to white paper defaults when WP mode selected
        and our defaults when it is deselected
        """
        if self.white_paper():
            self.t1_num.SetValue(1.65)
            self.bat_num.SetValue(0)
        else:
            self.t1_num.SetValue(1.3)
            self.bat_num.SetValue(1.3)
        self.calibration.update()
        self.update()

    def labelling_changed(self, pasl):
        """
        Update values to PASL/pCASL defaults when labelling method changed
        """
        if pasl:
            self.bat_num.SetValue(0.7)
            self.ie_num.SetValue(0.98)
        else:
            self.bat_num.SetValue(1.3)
            self.ie_num.SetValue(0.85)
