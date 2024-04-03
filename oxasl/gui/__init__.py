#!/usr/bin/env python
"""
Simple wxpython based GUI front-end to OXASL pipeline

Currently this does not use the GUI system from the FSL python libraries.
Possible improvements would include:

    - Use props library to hold run options and use the signalling mechanisms to communicate
    values. The built-in widget builder does not seem to be flexible enough however.
    - Use fsleyes embedded widget as the preview for a nicer, more interactive data preview

Requirements (beyond OXASL requirements):
    - wxpython
    - matplotlib

Copyright (c) 2008 University of Nottingham
"""

import os

import wx

from oxasl.gui.analysis_tab import AslAnalysis
from oxasl.gui.structure_tab import StructureTab
from oxasl.gui.calib_tab import AslCalibration
from oxasl.gui.input_tab import AslInputOptions
from oxasl.gui.dist_corr_tab import AslDistCorr
from oxasl.gui.widgets import PreviewPanel
from oxasl.gui.run_box import AslRun

class AslGui(wx.Frame):
    """
    Main GUI window
    """

    def __init__(self):
        wx.Frame.__init__(self, None, title="OXASL", size=(1200, 700), style=wx.DEFAULT_FRAME_STYLE)
        main_panel = wx.Panel(self)
        main_vsizer = wx.BoxSizer(wx.VERTICAL)

        local_dir = os.path.abspath(os.path.dirname(__file__))
        icon = wx.Icon()
        icon.CopyFromBitmap(wx.Bitmap(os.path.join(local_dir, "icon.png"), wx.BITMAP_TYPE_ANY))
        self.SetIcon(icon)

        banner = wx.Panel(main_panel, size=(-1, 80))
        banner.SetBackgroundColour((57, 71, 121))
        banner_fname = os.path.join(local_dir, "banner.png")
        wx.StaticBitmap(banner, -1, wx.Bitmap(banner_fname, wx.BITMAP_TYPE_ANY))
        main_vsizer.Add(banner, 0, wx.EXPAND)

        hpanel = wx.Panel(main_panel)
        hsizer = wx.BoxSizer(wx.HORIZONTAL)
        notebook = wx.Notebook(hpanel, id=wx.ID_ANY, style=wx.BK_DEFAULT)
        hsizer.Add(notebook, 1, wx.ALL|wx.EXPAND, 5)
        self.preview = PreviewPanel(hpanel)
        hsizer.Add(self.preview, 1, wx.EXPAND)
        hpanel.SetSizer(hsizer)
        main_vsizer.Add(hpanel, 2, wx.EXPAND)

        self.run_panel = wx.Panel(main_panel)
        runsizer = wx.BoxSizer(wx.HORIZONTAL)
        self.run_label = wx.StaticText(self.run_panel, label="Unchecked")
        self.run_label.SetFont(wx.Font(12, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        runsizer.Add(self.run_label, 1, wx.EXPAND)
        self.run_btn = wx.Button(self.run_panel, label="Run")
        runsizer.Add(self.run_btn, 0, wx.ALIGN_CENTER_VERTICAL)
        self.run_panel.SetSizer(runsizer)
        main_vsizer.Add(self.run_panel, 0, wx.EXPAND)
        self.run_panel.Layout()

        main_panel.SetSizer(main_vsizer)

        self.run = AslRun(self, self.run_btn, self.run_label)
        setattr(self.run, "preview", self.preview)
        tab_cls = [AslInputOptions, StructureTab, AslCalibration, AslDistCorr, AslAnalysis,]
        tabs = [cls(notebook, idx, len(tab_cls)) for idx, cls in enumerate(tab_cls)]

        for tab in tabs:
            notebook.AddPage(tab, tab.title)
            setattr(tab, "run", self.run)
            setattr(tab, "preview", self.preview)
            setattr(self.run, tab.name, tab)
            setattr(self.preview, tab.name, tab)
            for tab2 in tabs:
                if tab != tab2:
                    setattr(tab, tab2.name, tab2)
            tab.update()

        self.Layout()

def main():
    """
    Program entry point
    """
    app = wx.App(redirect=False)
    top = AslGui()
    top.Show()
    app.MainLoop()

if __name__ == '__main__':
    main()
