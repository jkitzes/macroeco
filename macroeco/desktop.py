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
import threading as thread

import main
from misc import get_log

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
        self.dirname = '.'
        self.parampath = 'parameters.txt'
        self.InitUI()
        self.Show(True)


    def InitUI(self):

        # Header
        sizerhead =  wx.BoxSizer(wx.HORIZONTAL)
        head_font = wx.Font(18, wx.SWISS, wx.NORMAL, wx.BOLD)
        heading = wx.StaticText(self, label='Macroeco Desktop')
        sizerhead.Add(heading, 0, wx.EXPAND)
        #heading.SetFont(head_font)

        # Step 1
        sizer1 = wx.BoxSizer(wx.VERTICAL)
        sizer1a = wx.BoxSizer(wx.HORIZONTAL)
        sizer1b = wx.BoxSizer(wx.HORIZONTAL)

        choose_text = wx.StaticText(self,
                                    label='1. Choose a parameters file'+' '*20)

        choose_button = wx.Button(self, label='Open')
        self.Bind(wx.EVT_BUTTON, self.OnOpen, choose_button)
        
        # Make attribute so easily modified by other methods
        self.choose_msg = wx.StaticText(self,
                                   label='')
        #self.choose_msg.SetFont(wx.Font(11, wx.SWISS, wx.SLANT, wx.NORMAL))

        sizer1a.Add(choose_text, 1, wx.EXPAND)
        sizer1a.Add(choose_button, 0, wx.EXPAND)
        sizer1b.Add(self.choose_msg, 1, wx.EXPAND)

        sizer1.Add(sizer1a, 0, wx.EXPAND)
        sizer1.Add(sizer1b, 0, wx.EXPAND)

        # Step 2
        sizer2 = wx.BoxSizer(wx.HORIZONTAL)
        run_text = wx.StaticText(self,
                                    label='2. Run analysis')
        self.run_button = wx.Button(self, label='Run')
        sizer2.Add(run_text, 1, wx.EXPAND)
        sizer2.Add(self.run_button, 0, wx.EXPAND)

        # Updating process
        self.process = None
        self.Bind(wx.EVT_BUTTON, self.OnRun, self.run_button)

        # Output window
        sizerlogbox = wx.BoxSizer(wx.HORIZONTAL)
        self.logbox = wx.TextCtrl(self, wx.ID_ANY, size=(400,400),
                           style = wx.TE_MULTILINE|wx.TE_READONLY|wx.HSCROLL)
        sizerlogbox.Add(self.logbox, 1, wx.EXPAND)

        # redirect text here
        redir=RedirectText(self.logbox)
        sys.stdout=redir
        sys.stderr=redir

        # Restore run button
        self.Bind(wx.EVT_IDLE, self.OnIdle)

        # All items
        sizer_main = wx.BoxSizer(wx.VERTICAL)
        sizer_main.Add(sizerhead, 0, wx.EXPAND | wx.ALL, 12)
        sizer_main.Add(sizer1, 0, wx.EXPAND | wx.ALL, 12)
        sizer_main.Add(sizer2, 0, wx.EXPAND | wx.ALL, 12)
        sizer_main.Add(sizerlogbox, 0, wx.EXPAND | wx.ALL, 12)

        # Set up main layout
        self.SetSizer(sizer_main)
        self.SetAutoLayout(True)
        sizer_main.Fit(self)

    def OnOpen(self,e):
        self.filename = ''
        self.dirname = ''
        dlg = wx.FileDialog(self, 'Choose a parameters file', self.dirname,
                            '', '*.*', wx.OPEN)
        if dlg.ShowModal() == wx.ID_OK:
            self.filename = dlg.GetFilename()
            self.dirname = dlg.GetDirectory()
            self.choose_msg.SetLabel('       Parameters file selected')
            self.parampath = os.path.join(self.dirname, self.filename)
        dlg.Destroy()

    def OnRun(self,e):
        self.logbox.SetValue('')
        self.RunMain()

    def RunMain(self):
        self.run_button.Enable(False)  # Turn the run button off
        self.t = thread.Thread(target=main.main, args=(self.parampath,))
        self.t.daemon = True  # Kills thread if app exits
        self.t.start()

    def OnIdle(self, event):
        if self.t:  # If a thread has been started
            if not self.t.is_alive():  # And it's not alive
                self.run_button.Enable(True)  # Turn the run button on

if __name__ == '__main__':
    # To execute, run `pythonw -m macroeco.desktop from root macroeco dir.
    app = wx.App(False)
    frame = MainWindow(None, 'Macroeco Desktop')
    app.MainLoop()
