"""
Macroeco Desktop - A graphical interface for macroeco

Open file dialog
http://wiki.wxpython.org/Getting%20Started

Redirecting stdout and stderr
http://blog.pythonlibrary.org/2009/01/01/wxpython-redirecting-stdout-stderr/

Process and stdout to window (see Example at link below)
http://wxpython.org/Phoenix/docs/html/Process.html#process
"""

import wx
import os, sys
from threading import Thread
from multiprocessing import Process
from _main import main
import logging

def launch():
    app = wx.App(False)
    frame = MainWindow(None, 'Macroeco Desktop')
    app.MainLoop()

class RedirectText(object):
    def __init__(self,aWxTextCtrl):
        self.out=aWxTextCtrl

    def write(self,string):
        wx.CallAfter(self.out.WriteText, string)

# Class for window
class MainWindow(wx.Frame):

    def __init__(self, parent, title):
        wx.Frame.__init__(self, parent, title=title)
        self.t = None
        self.filename = ''
        self.dirname = ''
        self.InitUI()
        self.Show(True)


    def InitUI(self):

        # Header
        sizerhead =  wx.BoxSizer(wx.HORIZONTAL)
        head_font = wx.Font(18, wx.SWISS, wx.NORMAL, wx.BOLD)
        heading = wx.StaticText(self, label='Macroeco Desktop')
        sizerhead.Add(heading, 0, wx.EXPAND)
        heading.SetFont(head_font)

        # Step 1
        sizer1 = wx.BoxSizer(wx.HORIZONTAL)

        param_text = wx.StaticText(self,
            label=("1. Open or create a parameter file\n"
                   "     File can be edited below and saved"))
        self.open_button = wx.Button(self, label='Open')
        self.new_button = wx.Button(self, label='New')
        self.save_button = wx.Button(self, label='Save')
        self.save_button.Enable(False)

        sizer1.Add(param_text, 1, wx.EXPAND)
        sizer1.Add(self.open_button, 0, wx.EXPAND | wx.RIGHT, 6)
        sizer1.Add(self.new_button, 0, wx.EXPAND | wx.RIGHT, 6)
        sizer1.Add(self.save_button, 0, wx.EXPAND)

        # Bind open and new buttons
        self.Bind(wx.EVT_BUTTON, self.OnOpen, self.open_button)
        self.Bind(wx.EVT_BUTTON, self.OnNew, self.new_button)
        self.Bind(wx.EVT_BUTTON, self.OnSave, self.save_button)

        # Param window
        sizerpfile = wx.BoxSizer(wx.HORIZONTAL)
        self.pfile = wx.TextCtrl(self, wx.ID_ANY, size=(600,300),
                           style=wx.TE_MULTILINE|wx.HSCROLL)
        sizerpfile.Add(self.pfile, 1, wx.EXPAND)

        # Step 2
        sizer2 = wx.BoxSizer(wx.HORIZONTAL)
        run_text = wx.StaticText(self, label='2. Save parameter file and '
            + 'Run analysis')
        self.run_button = wx.Button(self, label='Save and Run')
        sizer2.Add(run_text, 1, wx.EXPAND)
        sizer2.Add(self.run_button, 0, wx.EXPAND)

        # Bind run button
        self.Bind(wx.EVT_BUTTON, self.OnRun, self.run_button)

        # Output window
        sizerlogbox = wx.BoxSizer(wx.HORIZONTAL)
        self.logbox = wx.TextCtrl(self, wx.ID_ANY, size=(600,150),
                           style=wx.TE_MULTILINE|wx.TE_READONLY|wx.HSCROLL)
        sizerlogbox.Add(self.logbox, 1, wx.EXPAND)

        # Redirect text here
        redir = RedirectText(self.logbox)
        sys.stdout = redir
        sys.stderr = redir

        # Restore run button
        self.Bind(wx.EVT_IDLE, self.OnIdle)

        # All items

        sizer_main = wx.BoxSizer(wx.VERTICAL)
        sizer_main.Add(sizerhead, 0, wx.EXPAND | wx.ALL, 12)

        sizer_main.Add(sizer1, 0, wx.EXPAND | wx.ALL, 12)
        sizer_main.Add(sizerpfile, 0, wx.EXPAND | wx.ALL, 12)
        sizer_main.Add(sizer2, 0, wx.EXPAND | wx.ALL, 12)
        sizer_main.Add(sizerlogbox, 0, wx.EXPAND | wx.ALL, 12)

        # Set up main layout
        self.SetSizer(sizer_main)
        self.SetAutoLayout(True)
        sizer_main.Fit(self)

    def defaultFileDialogOptions(self):
        ''' Return a dictionary with file dialog options that can be
            used in both the save file dialog as well as in the open
            file dialog. '''
        return dict(message='Choose a file', defaultDir=self.dirname,
                    wildcard='*.*')

    def askUserForFilename(self, **dialogOptions):
        dialog = wx.FileDialog(self, **dialogOptions)
        if dialog.ShowModal() == wx.ID_OK:
            userProvidedFilename = True
            self.filename = dialog.GetFilename()
            self.dirname = dialog.GetDirectory()
        else:
            userProvidedFilename = False
        dialog.Destroy()
        return userProvidedFilename

    def OnOpen(self,e):
        if self.askUserForFilename(style=wx.OPEN,
                                   **self.defaultFileDialogOptions()):
            parampath = os.path.join(self.dirname, self.filename)
            f = open(parampath, 'r')
            self.pfile.SetValue(f.read())
            f.close()

            self.save_button.Enable(True)

        self.logbox.SetValue('')
        print "File opened at " + os.path.join(self.dirname, self.filename)

    def OnNew(self,e):
        if self.askUserForFilename(style=wx.SAVE,
                                   **self.defaultFileDialogOptions()):
            self.OnSave(e, new_file=True)
            self.save_button.Enable(True)

    def OnSave(self, event, new_file=False):
        f = open(os.path.join(self.dirname, self.filename), 'w')
        f.write(self.pfile.GetValue())
        f.close()

        self.logbox.SetValue('')
        if new_file:
            print "File created at "+os.path.join(self.dirname, self.filename)
        else:
            print "File saved at " + os.path.join(self.dirname, self.filename)

    def OnRun(self,e):
        self.OnSave(e)
        self.logbox.SetValue('')
        self.RunMain()

    def RunMain(self):
        self.run_button.Enable(False)  # Turn the run button off
        parampath = os.path.join(self.dirname, self.filename)

        self.t = Thread(target=main, args=(parampath,))
        self.t.daemon = True  # Kills thread if app exits
        self.t.start()
        # self.t.join()

    def OnIdle(self, event):
        if self.t:  # If a thread has been started
            if not self.t.is_alive():  # And it's not alive
                logging.shutdown() # Release lock on log file
                self.run_button.Enable(True)  # Turn the run button on
