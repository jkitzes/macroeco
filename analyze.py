from Tkinter import *
import subprocess
import glob
 

class Chooser:

    def __init__(self, master):
        self.structure={'data':"./data/*", 'code':"./code/*.py", 'projects':"./projects/*/"}
        Label(master,text="Datafile:").grid(row=0,column=0)
        Label(master,text="Analysis:").grid(row=0,column=1)
        Label(master,text="Output to:").grid(row=0,column=2)
        self.dlist = Listbox(master,exportselection=0,selectmode=SINGLE)
        self.dlist.grid(row=1,column=0,rowspan=10,padx=2)
        self.alist = Listbox(master,exportselection=0,selectmode=SINGLE)
        self.alist.grid(row=1,column=1,rowspan=10,padx=2)
        self.olist = Listbox(master,exportselection=0,selectmode=SINGLE)
        self.olist.grid(row=1,column=2,rowspan=5,padx=2)

        self.fill_list('data',self.dlist)
        self.fill_list('code',self.alist)
        self.fill_list('projects',self.olist)

        #TODO: add button for 'add project directory'
        
        self.quitbutton = Button(master, text="Exit", fg="red", command=master.quit).grid(row=12,column=0)
        

        self.runbutton = Button(master, text="Run", command=self.call_script).grid(row=12,column=2)


    def fill_list(self,target,list):
        '''Assumes a directory structure: data, code, projects are all siblings.'''
        files = glob.glob(self.structure[target])
        if "./data/archival" in files:
            hide_arch = files.index("./data/archival")
            files = files[:hide_arch] + files[hide_arch+1:]
        strip = self.structure[target].find('/',2)
        for item in files:
            list.insert(END, item[strip+1:].rstrip('/'))
        list.realcontent = files
        
        
    def call_script(self):
        data = '../.'+self.dlist.realcontent[int(self.dlist.curselection()[0])]
        script = '../.'+self.alist.realcontent[int(self.alist.curselection()[0])]
        output = self.olist.realcontent[int(self.olist.curselection()[0])]

        
        subprocess.Popen(["python", script, data],cwd=output,shell=False,stdin=None,stdout=None,close_fds=True) #Nb: this means everything had better be handled by individual scripts. TODO: design decision.

root = Tk()

app = Chooser(root)
root.title("AutoMETE")
root.mainloop()
