#!/usr/bin/python

'''Annotated, skeleton script showing how to run a particular
analysis with a reproducible workflow.

This is also the setup expected by the GUI.'''

import sys, os
from macroeco.workflow import Workflow
 
__author__ = "Chloe Lewis"
__copyright__ = "Copyright 2012, Regents of the University of California"
__license__ = None
__version__ = "0.5"
__maintainer__ = "Chloe Lewis"
__email__ = "chlewis@berkeley.edu"
__status__ = "Development"

# If your functions don't require any parameters, the Workflow object might
# still be useful to systematically track what script and data generated
# any output.

wf = Workflow()

wf.logger.debug('Something happening')
# If the analysis happens file-by-file:
for datafile, outputID in wf.single_datasets():
    wf.logger.debug('%s,%s'%(datafile, outputID))

# If the analysis combines all the files: 
for datafiles, outputID in wf.all_datasets():
    wf.logger.debug('%s,%s'%(datafiles, outputID))

        
            
    



