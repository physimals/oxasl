"""
oxasl.gui.structure_tab.py

Copyright (c) 2019 University of Oxford
"""
import os

import wx

from oxasl.gui.widgets import TabPage, OptionError

class StructureTab(TabPage):
    """
    Tab page containing options for providing structural data
    """

    EXISTING_FSLANAT = 0
    NEW_FSLANAT = 1
    INDEP_STRUC = 2
    NONE = 3

    TRANS_MATRIX = 0
    TRANS_IMAGE = 1
    TRANS_FSLANAT = 2

    def __init__(self, parent, idx, n):
        TabPage.__init__(self, parent, "Structure", idx, n, name="structure")

        self.section("Structure")

        self.struc_ch = wx.Choice(self, choices=[
            "Existing FSL_ANAT output",
            "Run FSL_ANAT on structural image",
            "Independent structural data",
            "None"
        ])
        self.struc_ch.SetSelection(self.NONE)
        self.struc_ch.Bind(wx.EVT_CHOICE, self._changed)
        self.struc_ch.span = 2
        self.pack("Structural data from", self.struc_ch)

        self.fsl_anat_picker = self.file_picker("Existing FSL_ANAT directory", pick_dir=True)
        self.struc_image_picker = self.file_picker("Structural Image")
        self.brain_image_picker = self.file_picker("Brain image", optional=True)

        self.section("Registration")

        self.transform_choices = ["Use matrix", "Use warp image", "Use FSL_ANAT"]

        self.transform_cb = wx.CheckBox(self, label="Transform to standard space")
        self.transform_cb.Bind(wx.EVT_CHECKBOX, self._changed)
        self.pack("", self.transform_cb, enable=False)

        self.sizer.AddGrowableCol(2, 1)
        self.SetSizer(self.sizer)
        self.next_prev()

    def options(self):
        options = {}
        if self.struc_ch.GetSelection() == self.INDEP_STRUC:
            options["struc"] = self.image("Structural data", self.struc_image_picker.GetPath())
            if self.brain_image_picker.checkbox.IsChecked():
                options["struc_brain"] = self.image("Structural brain", self.brain_image_picker.GetPath())
        elif self.struc_ch.GetSelection() == self.EXISTING_FSLANAT:
            options["fslanat"] = self.fsl_anat_picker.GetPath()
            if not os.path.exists(options["fslanat"]):
                raise OptionError("FSL_ANAT directory does not exist")
        # FIXME new FSL_ANAT

        if self.transform_cb.IsChecked():
            options["output_std"] = True
        return options

    def update(self):
        mode = self.struc_ch.GetSelection()
        self.fsl_anat_picker.Enable(mode == self.EXISTING_FSLANAT)
        self.struc_image_picker.Enable(mode in (self.NEW_FSLANAT, self.INDEP_STRUC))

        self.brain_image_picker.checkbox.Enable(mode == self.INDEP_STRUC)
        self.brain_image_picker.Enable(mode == self.INDEP_STRUC and self.brain_image_picker.checkbox.IsChecked())

        TabPage.update(self)
        