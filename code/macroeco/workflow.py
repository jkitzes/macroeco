#!/usr/bin/python

'''Workflow object manages the details of a reproducible workflow within the METE system.

   The parts are:
       Workflow: to create at the beginning of each script; it uses
       Parameters: Find, or ask for, and store the run names and parameters asked for by an analysis.
       Logging: using the standard Python library
       

'''
import xml.etree.ElementTree as etree
import sys, os, logging
from matplotlib.mlab import csv2rec


__author__ = "Chloe Lewis"
__copyright__ = "Copyright 2012, Regents of the University of California"
__credits__ = []
__license__ = None
__version__ = "0.5"
__maintainer__ = "Chloe Lewis"
__email__ = "chlewis@berkeley.edu"
__status__ = "Development"

paramfile = "parameters.xml"
logfile   = "macroeco.log"
loggername = "macroeco"



class Workflow:
    '''Manages the details of a reproducible workflow within the METE system.'''
    def __init__(self, asklist):
        ''' Set up logging, data-paths, output IDs, parameters.'''

        if len(sys.argv) < 2:
            print 'Error: need a path to data.'
            sys.exit()

        # 0th argument to script is the analysis name as running, tidy up:
        sPath, sExt = os.path.splitext(sys.argv[0])
        self.script = os.path.split(sPath)[-1]


        # The rest of the arguments are the data files to analyze.
        # Read them in to a data structure:
        self.data = {}
        for dfile in sys.argv[1:]:
            print dfile
            dname, dext = os.path.splitext(os.path.split(dfile)[-1])
            self.data[dname] = csv2rec(dfile,names=None)


        # Log everything tersely to console, INFO w/date to file.
        logging.basicConfig(level=logging.DEBUG,
                            format='%(asctime)s | %(levelname)s | %(filename)s:%(lineno)d | %(message)s', datefmt='%H:%M:%S')
        self.logger = logging.getLogger(loggername)
        fh = logging.FileHandler(logfile)
        fh.setLevel(logging.INFO)
        formatter = logging.Formatter('%(asctime)s | %(levelname)s | %(filename)s:%(lineno)d | %(message)s')
        fh.setFormatter(formatter)
        self.logger.addHandler(fh)

        #may need parameters, which may be in multiple runs
        self.runs = Parameters(self.script, asklist) #The asking stuff happens post-beta.
        
    def single_dataset_ID(self):
        #make an outputID
        #when do we know if we run a single dataset????
        pass
    
    def all_datasets_ID(self):
        pass     

        
class AllEntities:
    def __getitem__(self, key):
        return key        

class Parameters:
    '''Parameter values for any analysis:
        asked for as a dictionary of name:helpstring,
        available as self.params, a dictionary of name:value and "run_name":runname,
        written to the current working directory in %s.
    Parameters also tracks the run name, interactivity, and multiplicity.

    If %s is not sufficient, a dialog-box asks the user for values.'''%(paramfile, paramfile)
    
    def __init__(self, scriptname, asklist={}):
        '''Builds a dictionary of dictionaries to satisfy the asklist.

        Outer dictionary names are run names; inner are parameter name:value pairs.

        Finds parameters in this order:

            If there are interactive runs in %s for this scriptname:
                open populated dialog, return results with run_name.

            If there are only noninteractive runs in %s:
                if they satisfy the asklist for this scriptname:
                    return parameters with run_names added.
                    

        The argument asklist is a dictionary of 'name':'helpstring'.
        Helpstrings should be short and explain what kind of value is needed, e.g., 
            string,
            value in a range,
            value from a list.'''%(paramfile, paramfile)

        assert type(asklist) == type({}) #TODO: Integration test
        self.interactive = None
        self.script = scriptname
        self.params = {}
        self.logger = logging.getLogger(loggername)
        
        if len(asklist) == 0:
            return

        self.read_from_xml(asklist) 
        if self.is_asklist_fulfilled:
            return

        later='''        self.get_from_dialog(asklist)
        if self.is_asklist_fulfilled:
            self.write_to_xml()
            return
        
            else:
            raise Error #user wasn't helpful... TODO error type'''
         
    def read_from_xml(self, asklist):

        try:
            pf = open(paramfile,'r')
            pf.close()
        except IOError:
            self.logger.error('Could not open %s'%paramfile) #Note; can't write is also an IOError 

        parser = etree.XMLParser() # Without this, parsing works in iPython, console, not script. ??
        parser.parser.UseForeignDTD(True)
        parser.entity = AllEntities()
        try:
            pml = etree.parse(paramfile, parser=parser).getroot() #Depends on cwd - #TODO integration test
        except etree.ParseError:
            self.logger.error('ParseError trying to read %s'%paramfile)
        except:
            self.logger.error(sys.exc_info()[0])
        runcount = 0
        if len(pml) > 0:
            for child in pml:
                # check if any of these match scriptname TODO: and version (later)
                if child.get('scriptname') == self.script:
                    analysis = child
                    if 'interactive' in child.attrib:
                        ia = child.get('interactive')
                        if ia == 'F' or ia == 'False' or ia == 'f' or ia == 'false':
                            self.interactive = False
                        else:
                            self.interactive = True
                    else:
                        self.interactive = False #TODO: consider the default.
                    if len(analysis) > 0:
                        for run in analysis.getchildren():
                            if 'name' in run.attrib:
                                current_run = run.get('name')
                            else:
                                current_run = 'autoname'+str(runcount)
                                runcount += 1
                            self.params[current_run] = {}
                            for elt in run.getchildren():
                                if elt.tag == 'param':
                                    self.params[current_run][elt.get('name')] = elt.get('value')
        else:
            self.logger.error("Need run entries in %s"%paramfile)

        # if none, *or* if interactive, put up dialog asking for values for all the params

        # on dialog OK, write out  -- replace run with same name (tricky!)

    def is_asklist_fulfilled(self):
        return set(self.params.keys()).issubset(set.self.asklist.keys())
            

