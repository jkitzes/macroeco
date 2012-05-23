#!/usr/bin/python

'''Annotated, skeleton script showing how to run a particular
analysis with a reproducible workflow.

This is also the setup expected by the GUI.'''

import sys, os
from macroeco.workflow import Workflow
from macroeco import test

 
__author__ = "Chloe Lewis"
__copyright__ = "Copyright 2012, Regents of the University of California"
__license__ = None
__version__ = "0.5"
__maintainer__ = "Chloe Lewis"
__email__ = "chlewis@berkeley.edu"
__status__ = "Development"


# The workflow object takes care of logging, parameters, multiple runs, interactivity
# Filenames are read from  the arguments used when this code was called (e.g., at the terminal)
# Parameters are found in parameters.xml in the current working directory of the script

wf = Workflow()
# If the analysis happens file-by-file:
for datafile, outputID, params in wf.single_datasets():
    test.single_patch(datafile, outputID, params, wf.runs.interactive)

# If the analysis combines all the files: 
for datafiles, outputID, params in wf.all_datasets():
    test.multiple_patch(datafiles, outputID, params)


#Examples of logging:
wf.logger.debug('Debug is logged to console only') #To console only
wf.logger.info('Info and error are logged to file in working directory')

        
            
    



