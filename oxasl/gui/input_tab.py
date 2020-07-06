"""
oxasl.gui.structure_tab.py

Copyright (c) 2019 University of Oxford
"""
import os

import wx
import wx.grid

from oxasl.image import AslImage
from oxasl.gui.widgets import TabPage, NumberChooser, NumberList

class AslInputOptions(TabPage):
    """
    Tab page for specifying the ASL input data
    """

    def __init__(self, parent, idx, n):
        TabPage.__init__(self, parent, "Input Data", idx, n, name="input")

        self.groups = ["PLDs", "Repeats", "Label/Control pairs"]
        self.abbrevs = ["t", "r", "l"]

        self.section("Data contents")

        self.data_picker = self.file_picker("Input Image", handler=self.set_default_dir)
        self.ntis_int = self.integer("Number of PLDs", minval=1, maxval=100, initial=1)
        self.nrepeats_label = wx.StaticText(self, label="<Unknown>")
        self.pack("Number of repeats", self.nrepeats_label)

        self.section("Data order")

        self.choice1 = wx.Choice(self, choices=self.groups)
        self.choice1.SetSelection(2)
        self.choice1.Bind(wx.EVT_CHOICE, self._changed)
        self.choice2 = wx.Choice(self, choices=self.groups)
        self.choice2.SetSelection(0)
        self.choice2.Bind(wx.EVT_CHOICE, self._changed)
        self.pack("Grouping order", self.choice1, self.choice2)
        self.tc_ch = self.choice("Label/Control pairs", choices=["Label then control", "Control then label"],
                                 optional=True, initial_on=True)

        self.section("Acquisition parameters")

        self.labelling_ch = self.choice("Labelling", choices=["pASL", "cASL/pcASL"], initial=1,
                                        handler=self.labelling_changed)

        self.bolus_dur_ch = wx.Choice(self, choices=["Constant", "Variable"])
        self.bolus_dur_ch.SetSelection(0)
        self.bolus_dur_ch.Bind(wx.EVT_CHOICE, self._changed)
        self.bolus_dur_num = NumberChooser(self, minval=0, maxval=2.5, step=0.1, initial=1.8)
        self.bolus_dur_num.span = 2
        self.bolus_dur_num.spin.Bind(wx.EVT_SPINCTRL, self.bolus_dur_changed)
        self.bolus_dur_num.slider.Bind(wx.EVT_SLIDER, self.bolus_dur_changed)
        self.pack("Bolus duration (s)", self.bolus_dur_ch, self.bolus_dur_num)

        self.bolus_dur_list = NumberList(self, self.ntis())
        self.bolus_dur_list.span = 3
        self.bolus_dur_list.Bind(wx.grid.EVT_GRID_CELL_CHANGED, self._changed)
        self.pack("Bolus durations (s)", self.bolus_dur_list, enable=False)

        self.ti_list = NumberList(self, self.ntis())
        self.ti_list.span = 3
        self.ti_list.Bind(wx.grid.EVT_GRID_CELL_CHANGED, self._changed)
        self.pack("PLDs (s)", self.ti_list)

        self.readout_ch = wx.Choice(self, choices=["3D (eg GRASE)", "2D multi-slice (eg EPI)"])
        self.readout_ch.SetSelection(0)
        self.readout_ch.Bind(wx.EVT_CHOICE, self._changed)
        self.time_per_slice_num = NumberChooser(self, label="Time per slice (ms)",
                                                minval=0, maxval=50, step=1, initial=10)
        self.time_per_slice_num.span = 2
        self.pack("Readout", self.readout_ch, self.time_per_slice_num)
        self.time_per_slice_num.Enable(False)

        self.multiband_cb = wx.CheckBox(self, label="Multi-band")
        self.multiband_cb.Bind(wx.EVT_CHECKBOX, self._changed)
        self.slices_per_band_spin = wx.SpinCtrl(self, min=1, max=100, initial=5)
        self.slices_per_band_label = wx.StaticText(self, label="slices per band")
        self.pack("", self.multiband_cb, self.slices_per_band_spin, self.slices_per_band_label, enable=False)
        self.multiband_cb.Enable(False)

        self.sizer.AddGrowableCol(2, 1)
        self.SetSizer(self.sizer)
        self.next_prev()

    def options(self):
        options = {
            "order" : self.preview.order_preview.order,
            "casl" : self.labelling_ch.GetSelection() == 1,
        }

        if not self.tc_pairs():
            options["iaf"] = "diff"
        elif self.tc_ch.GetSelection() == 0:
            options["iaf"] = "tc"
        else:
            options["iaf"] = "ct"

        if options["casl"]:
            options["plds"] = self.ti_list.GetValues()
        else:
            options["tis"] = self.ti_list.GetValues()

        if self.bolus_dur_multi():
            options["taus"] = self.bolus_dur_list.GetValues()
        else:
            options["tau"] = self.bolus_dur_num.GetValue()

        if self.readout_2d():
            options["slicedt"] = self.time_per_slice_num.GetValue()

        if self.multiband():
            options["sliceband"] = self.slices_per_band_spin.GetValue()

        asldata = AslImage(self.image("Input data", self.data_picker.GetPath()).data, **options)

        # Update preview widget on successful AslImage construction
        nrpts = asldata.rpts[0]
        self.nrepeats_label.SetLabel("%i" % nrpts)
        self.preview.order_preview.n_tis = asldata.ntis
        self.preview.order_preview.n_repeats = nrpts
        self.preview.order_preview.tc_pairs = asldata.iaf in ("tc", "ct")
        self.preview.order_preview.tagfirst = asldata.iaf == "tc"
        self.preview.order_preview.Refresh()

        return {"asldata" : asldata}

    def ntis(self):
        """
        :return: Number of TIs/PLDs selected
        """
        return self.ntis_int.GetValue()

    def tc_pairs(self):
        """
        :return: True if data is tag-control pairs (in either order)
        """
        return self.tc_ch.checkbox.IsChecked()

    def bolus_dur_multi(self):
        """
        :return: True if multiple bolus durations are being used
        """
        return self.bolus_dur_ch.GetSelection() == 1

    def readout_2d(self):
        """
        :return: True if readout is 2D multi-slice
        """
        return self.readout_ch.GetSelection() == 1

    def multiband(self):
        """
        :return: True if multiband readout is selected
        """
        return self.multiband_cb.IsChecked()

    def set_default_dir(self, _):
        """
        Bit of a hack - set the default dir for other file pickers to the same dir
        as the main data
        """
        d = os.path.dirname(self.data_picker.GetPath())
        for w in [self.structure.fsl_anat_picker,
                  self.structure.struc_image_picker,
                  self.structure.brain_image_picker,
                  self.calibration.calib_image_picker,
                  self.calibration.ref_tissue_mask_picker,
                  self.calibration.coil_image_picker,
                  self.analysis.mask_picker,
                  self.distcorr.calib_picker,
                  self.distcorr.fmap_picker,
                  self.distcorr.fmap_mag_picker,
                  self.distcorr.fmap_be_picker]:
            try:
                # WX version dependent...
                w.SetInitialDirectory(d)
            except:
                w.SetPath(d)
        self.update()

    def update(self):
        self.ti_list.set_size(self.ntis())
        self.bolus_dur_list.set_size(self.ntis())

        self.time_per_slice_num.Enable(self.readout_2d())
        self.multiband_cb.Enable(self.readout_2d())
        self.slices_per_band_spin.Enable(self.multiband() and self.readout_2d())
        self.slices_per_band_label.Enable(self.multiband() and self.readout_2d())

        self.bolus_dur_num.Enable(not self.bolus_dur_multi())
        self.bolus_dur_list.Enable(self.bolus_dur_multi())

        self.tc_ch.Enable(self.tc_pairs())
        self.update_groups()

        TabPage.update(self)

    def labelling_changed(self, event):
        """
        Update naming of delay times depending on whether PASL or pCASL labelling is selected
        """
        if event.GetInt() == 0:
            self.bolus_dur_num.SetValue(0.7)
            self.ntis_int.label.SetLabel("Number of TIs")
            self.ti_list.label.SetLabel("TIs")
            self.preview.order_preview.tis_name = "TIs"
            self.groups[0] = "TIs"
        else:
            self.bolus_dur_num.SetValue(1.8)
            self.ntis_int.label.SetLabel("Number of PLDs")
            self.ti_list.label.SetLabel("PLDs")
            self.preview.order_preview.tis_name = "PLDs"
            self.groups[0] = "PLDs"
        self.analysis.labelling_changed(event.GetInt() == 0)
        self.update()

    def bolus_dur_changed(self, event):
        """ If constant bolus duration is changed, update the disabled list of
            bolus durations to match, to avoid any confusion """
        if self.bolus_dur_type() == 0:
            for c in range(self.bolus_dur_list.n):
                self.bolus_dur_list.SetCellValue(0, c, str(self.bolus_dur()[0]))
        event.Skip()

    def update_groups(self):
        """
        This hideous code is to update the choices available in the ordering menus,
        hide the second if it is not relevant (1 PLD) and derive the data ordering
        string to pass to the data order preview
        """
        g1 = self.choice1.GetString(self.choice1.GetSelection())
        if self.choice2.IsShown():
            g2 = self.choice2.GetString(self.choice2.GetSelection())
        else:
            g2 = self.groups[0]
        if g1 == g2:
            g2 = self.groups[(self.groups.index(g1) + 1) % 3]

        self.choice2.Show()
        choices_st, choices_end = 0, 3

        # If data is not pairs, don't offer TC pairs options
        # Change selections to something else, but make sure they aren't the same
        if not self.tc_pairs():
            self.choice2.Hide()
            if g1 == self.groups[2]:
                g1 = self.groups[0]
            if g2 in (g1, self.groups[2]):
                g2 = self.groups[1-self.groups.index(g1)]
            choices_end = 2

        # If only one TI/PLD, don't offer to group by TIs/PLDs
        # Change selections to something else. Second menu is set to group by PLDs
        # but is invisible now
        if self.ntis() == 1:
            self.choice2.Hide()
            if g1 == self.groups[0]:
                g1 = self.groups[1]
            g2 = self.groups[0]
            choices_st = 1

        group1_items = []
        group2_items = []
        for item in self.groups[choices_st:choices_end]:
            group1_items.append(item)
            if item != g1:
                group2_items.append(item)

        self.update_group_choice(self.choice1, group1_items, g1)
        self.update_group_choice(self.choice2, group2_items, g2)

        idx1 = self.groups.index(g1)
        idx2 = self.groups.index(g2)

        order = self.abbrevs[idx1]
        order += self.abbrevs[idx2]
        order += self.abbrevs[3-idx1-idx2]
        self.preview.order_preview.order = order

        # Need to do this as we may have unhidden the second menu
        self.GetSizer().Layout()

    def update_group_choice(self, w, items, sel):
        """
        Update a grouping menu to include a specified set of entries
        """
        w.Enable(False)
        w.Clear()
        w.AppendItems(items)
        w.SetSelection(w.FindString(sel))
        w.Enable(True)
