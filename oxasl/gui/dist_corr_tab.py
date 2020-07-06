"""
oxasl.gui.dist_corr_tab.py

Copyright (c) 2019 University of Oxford
"""
import wx

from oxasl.gui.widgets import TabPage

class AslDistCorr(TabPage):
    """
    Tab page containing options for distortion correction
    """

    FIELDMAP = 0
    CALIB_IMAGE = 1

    def __init__(self, parent, idx, n):
        TabPage.__init__(self, parent, "Distortion Correction", idx, n, name="distcorr")

        self.distcorr_choices = ["Fieldmap", "Calibration image"]

        self.section("Distortion Correction")

        # Calibration image options
        self.distcorr_cb = wx.CheckBox(self, label="Apply distortion correction")
        self.distcorr_cb.Bind(wx.EVT_CHECKBOX, self._changed)

        self.distcorr_ch = wx.Choice(self, choices=self.distcorr_choices[:1])
        self.distcorr_ch.SetSelection(self.FIELDMAP)
        self.distcorr_ch.Bind(wx.EVT_CHOICE, self._changed)
        self.pack("", self.distcorr_cb, self.distcorr_ch, enable=False)

        # Calib image options
        self.section("Calibration Image Mode")
        self.calib_picker = self.file_picker("Phase-encode-reversed calibration image")

        # Fieldmap options
        self.section("Fieldmap Mode")
        self.fmap_picker = self.file_picker("Fieldmap image (in rad/s)")
        self.fmap_mag_picker = self.file_picker("Fieldmap magnitude image")
        self.fmap_be_picker = self.file_picker("Brain-extracted magnitude image", optional=True)

        # General options
        self.section("General")
        self.echosp_num = self.number("Effective EPI echo spacing", minval=0, maxval=10)
        self.pedir_ch = self.choice("Phase encoding direction", choices=["x", "y", "z", "-x", "-y", "-z"])

        self.sizer.AddGrowableCol(1, 1)
        #sizer.AddGrowableRow(5, 1)
        self.SetSizer(self.sizer)
        self.next_prev()

    def options(self):
        options = {}
        if self.distcorr():
            options.update({
                "echospacing" : self.echosp_num.GetValue(),
                "pedir" : self.pedir_ch.GetStringSelection(),
            })
            if self.distcorr_fmap():
                options.update({
                    "fmap" : self.image("Fieldmap", self.fmap_picker.GetPath()),
                    "fmapmag" : self.image("Fieldmap magnitude", self.fmap_mag_picker.GetPath()),
                    "fmapmagbrain" : self.image("Fieldmap brain", self.fmap_be_picker.GetPath()),
                })
            else:
                options.update({
                    "cblip" : self.image("Phase-reversed calibration data", self.calib_picker.GetPath()),
                })
        return options

    def distcorr(self):
        """
        :return: True if distortion correction is enabled
        """
        return self.distcorr_cb.IsChecked()

    def distcorr_fmap(self):
        """
        :return: True if field distortion correction is selected
        """
        return self.distcorr_ch.GetSelection() == self.FIELDMAP

    def update(self):
        self.distcorr_ch.Enable(self.distcorr())

        cal = self.distcorr() and not self.distcorr_fmap()
        self.calib_picker.Enable(cal)

        fmap = self.distcorr() and self.distcorr_fmap()
        self.fmap_picker.Enable(fmap)
        self.fmap_mag_picker.Enable(fmap)
        self.fmap_be_picker.Enable(fmap)

        self.pedir_ch.Enable(self.distcorr())
        self.echosp_num.Enable(self.distcorr())

        TabPage.update(self)

    def calib_changed(self, enabled):
        """
        If calibration enabled, add the calibration image option for distortion correction
        """
        sel = self.distcorr_ch.GetSelection()
        if enabled:
            choices = self.distcorr_choices
            sel = 1
        else:
            choices = self.distcorr_choices[:1]
            sel = 0
        self.distcorr_ch.Enable(False)
        self.distcorr_ch.Clear()
        self.distcorr_ch.AppendItems(choices)
        self.distcorr_ch.SetSelection(sel)
        self.update()
