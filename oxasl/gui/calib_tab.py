"""
oxasl.gui.calibration_tab.py

Copyright (c) 2019 University of Oxford
"""

from oxasl.gui.widgets import TabPage

class AslCalibration(TabPage):
    """
    Tab page containing options for calibration
    """

    def __init__(self, parent, idx, n):
        TabPage.__init__(self, parent, "Calibration", idx, n)

        self.calib_cb = self.checkbox("Enable Calibration", bold=True, handler=self.calib_changed)

        self.calib_image_picker = self.file_picker("Calibration Image")
        self.seq_tr_num = self.number("Sequence TR (s)", minval=0, maxval=10, initial=6)
        self.calib_gain_num = self.number("Calibration Gain", minval=0, maxval=5, initial=1)
        self.calib_mode_ch = self.choice("Calibration mode", choices=["Reference Region", "Voxelwise"])

        self.section("Reference tissue")

        self.ref_tissue_type_ch = self.choice("Type", choices=["CSF", "WM", "GM", "None"],
                                              handler=self.ref_tissue_type_changed)
        self.ref_tissue_mask_picker = self.file_picker("Mask", optional=True)
        self.ref_t1_num = self.number("Reference T1 (s)", minval=0, maxval=5, initial=4.3)
        self.seq_te_num = self.number("Sequence TE (ms)", minval=0, maxval=30, initial=0)
        self.ref_t2_num = self.number("Reference T2 (ms)", minval=0, maxval=1000, initial=750, step=10)
        self.blood_t2_num = self.number("Blood T2 (ms)", minval=0, maxval=1000, initial=150, step=10)
        self.coil_image_picker = self.file_picker("Coil Sensitivity Image", optional=True)

        self.sizer.AddGrowableCol(2, 1)
        self.SetSizer(self.sizer)
        self.next_prev()

    def options(self):
        options = {}
        if self.calib():
            options.update({
                "calib" : self.image("Calibration data", self.calib_image_picker.GetPath()),
                "calib_gain" : self.calib_gain_num.GetValue(),
                "tr" : self.seq_tr_num.GetValue(),
            })

            if self.refregion():
                options.update({
                    "calib_method" : "refregion",
                    "tissref" : self.ref_tissue_type_ch.GetString(self.ref_tissue_type()),
                    "te" : self.seq_te_num.GetValue(),
                    "t1r" : self.ref_t1_num.GetValue(),
                    "t2r" : self.ref_t2_num.GetValue(),
                    "t2b" : self.blood_t2_num.GetValue(),
                })
                if self.ref_tissue_mask_picker.checkbox.IsChecked():
                    options["refmask"] = self.ref_tissue_mask_picker.GetPath()
            else:
                options["calib_method"] = "voxelwise"

            if self.coil_image_picker.checkbox.IsChecked():
                options["cref"] = self.image("Calibration reference data", self.coil_image_picker.GetPath())

        return options

    def calib(self):
        """
        :return: True if calibration is enabled
        """
        return self.calib_cb.IsChecked()

    def refregion(self):
        """
        :return: True if reference region calibration is selected
        """
        return self.calib_mode_ch.GetSelection() == 0

    def ref_tissue_type(self):
        """
        :return reference tissue type index
        """
        return self.ref_tissue_type_ch.GetSelection()

    def ref_tissue_type_changed(self, _):
        """
        Update reference tissue parameters to currently selected reference tissue type
        """
        if self.ref_tissue_type() == 0: # CSF
            self.ref_t1_num.SetValue(4.3)
            self.ref_t2_num.SetValue(750)
        elif self.ref_tissue_type() == 1: # WM
            self.ref_t1_num.SetValue(1.0)
            self.ref_t2_num.SetValue(50)
        elif self.ref_tissue_type() == 2: # GM
            self.ref_t1_num.SetValue(1.3)
            self.ref_t2_num.SetValue(100)
        self.update()

    def calib_changed(self, _):
        """
        Update option visibility when calibration is enabled/disabled
        """
        self.distcorr.calib_changed(self.calib())
        self.update()

    def update(self):
        enable = self.calib()
        self.seq_tr_num.Enable(enable)
        self.calib_image_picker.Enable(enable)
        self.calib_gain_num.Enable(enable)
        self.coil_image_picker.checkbox.Enable(enable)
        if self.analysis.white_paper():
            self.calib_mode_ch.SetSelection(1)
        self.calib_mode_ch.Enable(enable and not self.analysis.white_paper())
        self.ref_tissue_type_ch.Enable(enable and self.refregion())

        if self.ref_tissue_type() == 3:
            # Ref tissue = None - enforce mask
            self.ref_tissue_mask_picker.checkbox.Enable(False)
            self.ref_tissue_mask_picker.checkbox.SetValue(enable and self.refregion())
            self.ref_tissue_mask_picker.Enable(enable and self.refregion())
        else:
            self.ref_tissue_mask_picker.checkbox.Enable(enable and self.refregion())
        self.ref_tissue_mask_picker.Enable(enable and self.refregion() and self.ref_tissue_mask_picker.checkbox.IsChecked())

        self.coil_image_picker.checkbox.Enable(enable and self.refregion())
        self.coil_image_picker.Enable(enable and self.refregion() and self.coil_image_picker.checkbox.IsChecked())
        self.seq_te_num.Enable(enable and self.refregion())
        self.blood_t2_num.Enable(enable and self.refregion())
        self.ref_t1_num.Enable(enable and self.refregion())
        self.ref_t2_num.Enable(enable and self.refregion())
        TabPage.update(self)
