#!/usr/bin/python

'''Find or ask for and store the run names and parameters asked for by an analysis. 

'''
import xml.etree.ElementTree as etree
import sys

__author__ = "Chloe Lewis"
__copyright__ = "Copyright 2012, Regents of the University of California"
__credits__ = []
__license__ = None
__version__ = "0.5"
__maintainer__ = "Chloe Lewis"
__email__ = "chlewis@berkeley.edu"
__status__ = "Development"

paramfile = "parameters.xml"
class AllEntities:
    def __getitem__(self, key):
        #key is your entity, you can do whatever you want with it here
        return key        

class Parameters:
    '''Parameter values for any analysis:
        asked for as a dictionary of name:helpstring,
        available as self.params, a dictionary of name:value and "run_name":runname,
        written to the current working directory in %s.
    Parameters also tracks the run name, interactivity, and multiplicity.

    If %s is not sufficient, a dialog-box asks the user for values.'''%(paramfile, paramfile)
    
    def __init__(self, scriptname, asklist):
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
            print 'IOError: Could not open %s'%paramfile #Note; can't write is also an IOError 

        parser = etree.XMLParser() # Without this, parsing works in iPython, console, not script. ??
        parser.parser.UseForeignDTD(True)
        parser.entity = AllEntities()
        try:
            pml = etree.parse(paramfile, parser=parser).getroot() #Depends on cwd - #TODO integration test
        except etree.ParseError:
            print 'ParseError trying to read %s'%paramfile
        except:
            print 'Error: unexpected error: ',sys.exc_info()[0]
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
            print "Error: need run entries for analysis %s in %s"%(self.script, paramfile)

        # if none, *or* if interactive, put up dialog asking for values for all the params

        # on dialog OK, write out  -- replace run with same name (tricky!)

    def is_asklist_fulfilled(self):
        return set(self.params.keys()).issubset(set.self.asklist.keys())
            

