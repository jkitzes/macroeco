#!/usr/bin/python

'''GUI to easily run MaxEnt macro-ecology analyses. 

Note: Depends on the file structure we ship with. 
'''



from Tkinter import *
from tkFileDialog   import askdirectory
import subprocess
import glob
import os



 
__author__ = "Chloe Lewis"
__copyright__ = "Copyright 2012, Regents of the University of California"
__license__ = None
__version__ = "0.5"
__maintainer__ = "Chloe Lewis"
__email__ = "chlewis@berkeley.edu"
__status__ = "Development"


class Chooser:
    '''GUI to choose the data, analysis script, and output directory.'''
    def __init__(self, master):
        self.structure={'data':"./data/*", 'code':"./code/*.py", 'output':os.path.abspath("./projects")}
        self.master = master
        Label(master,text="Datafile:").grid(row=0,column=0,padx=2)
        Label(master,text="Analysis:").grid(row=0,column=2,padx=2)
        Label(master,text="Results saved to:").grid(row=12,column=0,padx=2)
        self.output = StringVar()
        self.output.set(self.structure['output'])
        Label(master, textvariable=self.output).grid(row=13,column=0,columnspan=4,padx=2)
        self.dlist = Listbox(master,exportselection=0,selectmode=MULTIPLE)
        self.dlist.grid(row=1,column=0,rowspan=10,columnspan=2,padx=2)
        self.alist = Listbox(master,exportselection=0,selectmode=SINGLE)
        self.alist.grid(row=1,column=2,columnspan=2,rowspan=10,padx=2)

        self.fill_list('data',['./data/archival','./data/formatted'],self.dlist)
        self.fill_list('code',[],self.alist)

        self.projectbutton = Button(master, text="Change output directory", command=self.change_project).grid(row=14, column=0)
        self.quitbutton = Button(master, text="Exit", fg="red", command=master.quit).grid(row=14,column=2)
        self.runbutton = Button(master, text="Run", command=self.call_script).grid(row=14,column=3)

        self.dlist.selection_set(0)
        self.alist.selection_set(2) #TODO: also 0 when dir structure cleaned up
    def change_project(self):
        '''Change the directory in which output will be saved.'''

        path = askdirectory()
        self.structure['output'] = path
        self.output.set(path)

    def fill_list(self,target,hide,list):
        '''Fills the GUI. Assumes that data, code, projects are all siblings.'''
        files = glob.glob(self.structure[target])
        for h in hide:
            hi = files.index(h)
            files = files[:hi] + files[hi+1:]
        strip = self.structure[target].find('/',2)
        for item in files:
            list.insert(END, item[strip+1:].rstrip('/'))
        list.realcontent = files
        
        
    def call_script(self):
        '''Runs the chosen data,analysis,output triple.

        Note that the GUI stays open: user can choose another triple to run.'''
        METEbase = os.path.dirname(os.path.abspath(__file__))
        # print self.dlist.realcontent[int(self.dlist.curselection()[0])]
        # TODO: figure out why I can't clip off leading './'. Do I need them at all?
        data = os.path.join(METEbase,self.dlist.realcontent[int(self.dlist.curselection()[0])])
        script = os.path.join(METEbase,self.alist.realcontent[int(self.alist.curselection()[0])])
        print METEbase, data, script
        output = self.structure['output']

        
        subprocess.Popen(["python", script, data],cwd=output,shell=False,stdin=None,stdout=None,close_fds=True) 

root = Tk()

app = Chooser(root)
root.title("AutoMETE")
root.mainloop()

testing = '''Human test-script:

are lists correctly populated? Check against filesystem.

is default project directory correct?

run one d vs one a:
    log & output in default project directory?

change output directory:
    readable and correct in GUI?
    run multiple data vs multiple analysis:
        log and output correctly named, in new project directory?"

