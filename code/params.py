
#!/usr/bin/python

'''Manage the setup, commission, and teardown for a run of an analysis. 

'''
import xml.etree.ElementTree as etree

__author__ = "Chloe Lewis"
__copyright__ = "Copyright 2012, Regents of the University of California"
__credits__ = ["John Harte"]
__license__ = None
__version__ = "0.5"
__maintainer__ = "Chloe Lewis"
__email__ = "chlewis@berkeley.edu"
__status__ = "Development"


def get_parameters(scriptname, asklist):
    '''Returns parameter values as a dictionary.

    If there is no parameters file in the cwd, put up a dialog asking for them.
    If parameters given in dialog, write out to parameters file.

    Scripts are responsible for calling get_parameters().

    The argument asklist is a dictionary of 'name':'helpstring'.
    Helpstrings should be short and explain what kind of value is needed: string, value in a range, value from a list.

    NOTE: this function does *not* check the parameters against the help string.'''
    
    if len(asklist) == 0:
        return {}

    askparams = {}
    for p in asklist.keys():
        askparams[p] = asklist[p]
        
    # check for parameters.xml
    try:
        # if file, read, load, present in dialog
        pml = etree.parse('parameters.xml').getroot() #Depending on cwd
        for child in pml:
            # check if any of these match scriptname TODO: and version (later)
            print child, child.attrib
            if child.attrib['scriptname'] == scriptname:
                print('found reference to %s'%scriptname)
                for run in child:
                    print run, run.attrib
                    for p in run:
                        print p, p.value, p.text

    except:
        print 'could not get parameters for %s from file'%scriptname

    # if none, *or* if interactive, put up dialog asking for values for all the params

    # on dialog OK, write out  -- replace run with same name (tricky!)

    return askparams
    


# Cases for testing:
# No params file. Dialog, write, reload.
# Params file matches ask. Check that read matches expected.
# Params file doesn't match ask. Dialog, write, reload, check against ask.
## A proper subset P
## P proper subset A
## Neither
# Types of param: string, int, float, lists of any of those.

if __name__=="__main__":
    print get_parameters('test2.py', {'title':"OldTitle",'layer':"Between 0 and 2.5",'species':"List of integers, eg [1,4, 12]"})
