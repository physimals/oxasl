"""
oxasl.gui.run_box.py

Component of the OXASL GUI which runs the OXASL pipeline

Copyright (c) 2019 University of Oxford
"""
import sys
import traceback
from threading import Thread

import nibabel as nib

import wx
from wx.lib.pubsub import pub

from oxasl import Workspace
from oxasl.oxford_asl import oxasl
from oxasl.gui.widgets import OptionError

class LogWriter():
    """
    File-like interface to send OXASL log output to the GUI output window
    """
    def write(self, line):
        wx.CallAfter(pub.sendMessage, "run_stdout", line=line)

    def flush(self):
        pass

class OxaslRunner(Thread):
    """
    Runs OXASL pipeline in a background thread
    """
    def __init__(self, options):
        Thread.__init__(self)
        self.options = options

    def run(self):
        ret = -1
        try:
            wsp = Workspace(log=LogWriter(), **self.options)
            oxasl(wsp)
            ret = 0
        finally:
            wx.CallAfter(pub.sendMessage, "run_finished", retcode=ret)

class AslRun(wx.Frame):
    """
    Runs the OXASL pipeline, displaying output in a window
    """

    def __init__(self, parent, run_btn, run_label):
        wx.Frame.__init__(self, parent, title="Run", size=(600, 400), style=wx.DEFAULT_FRAME_STYLE)

        self.options = None
        self.run_btn = run_btn
        self.run_btn.Bind(wx.EVT_BUTTON, self.dorun)
        self.run_label = run_label
        self.preview_data = None

        self.sizer = wx.BoxSizer(wx.VERTICAL)
        self.output_text = wx.TextCtrl(self, style=wx.TE_READONLY | wx.TE_MULTILINE)
        font = wx.Font(8, wx.TELETYPE, wx.NORMAL, wx.NORMAL)
        self.output_text.SetFont(font)
        self.sizer.Add(self.output_text, 1, flag=wx.EXPAND)

        self.SetSizer(self.sizer)
        self.Bind(wx.EVT_CLOSE, self.close)
        pub.subscribe(self.write_output, "run_stdout")
        pub.subscribe(self.finished, "run_finished")

    def write_output(self, line):
        """
        Write a line onto the log output window
        """
        self.output_text.AppendText(line)

    def close(self, _):
        """
        Hide the log output window
        """
        self.Hide()

    def finished(self, retcode):
        """
        Called by background thread when pipeline has finished (successful or not)
        """
        if retcode != 0:
            self.write_output("\nWARNING: command failed\n")
        self.update()

    def dorun(self, _):
        """
        When 'Run' button clicked, run the pipeline displaying the output in the
        log window
        """
        if self.options:
            self.Show()
            self.Raise()
            self.output_text.Clear()
            self.run_btn.Enable(False)
            self.run_label.SetForegroundColour(wx.Colour(0, 0, 128))
            self.run_label.SetLabel("Running - Please Wait")
            OxaslRunner(self.options).start()

    def update(self):
        """
        Get the OXASL options, and check they are valid. If not valid,
        display the first error in the status label
        """
        self.options = None
        try:
            self.options = self.get_options()
            self.run_label.SetForegroundColour(wx.Colour(0, 128, 0))
            self.run_label.SetLabel("Ready to Go")
            self.run_btn.Enable(True)
        except (OptionError, ValueError, nib.filebasedimages.ImageFileError) as exc:
            self.run_btn.Enable(False)
            self.run_label.SetForegroundColour(wx.Colour(255, 0, 0))
            self.run_label.SetLabel(str(exc))
        except:
            # Any other exception is a program bug - report it to STDERR
            self.run_btn.Enable(False)
            self.run_label.SetForegroundColour(wx.Colour(255, 0, 0))
            self.run_label.SetLabel("Unexpected error - see console and report as a bug")
            traceback.print_exc(sys.exc_info()[1])

    def get_preview_data(self):
        """
        Get a perfusion weighted image for the preview
        """
        asldata = self.input.options()["asldata"]
        return asldata.perf_weighted().data

    def get_options(self):
        """
        Get OXASL options

        Exception text is reported by the GUI
        """
        options = {
            "output_native" : True,
            "output_struc" : True,
            "save_mask" : True,
        }
        options.update(self.input.options())
        options.update(self.analysis.options())
        options.update(self.structure.options())
        options.update(self.distcorr.options())
        options.update(self.calibration.options())
        #print(options)

        return options
