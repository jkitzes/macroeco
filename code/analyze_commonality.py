#!/usr/bin/python

'''Script that uses reproducible workflow to test commonality predictions
on data sets.
'''

import sys, os
from macroeco.workflow import Workflow
from macroeco import empirical
from macroeco import graph
from macroeco import cmn_analysis

__author__ = "Mark Wilber"
__copyright__ = "Copyright 2012, Regents of the University of California"
__license__ = None
__version__ = "0.5"
__maintainer__ = "Mark Wilber"
__email__ = "mqw@berkeley.edu"
__status__ = "Development"

wf = Workflow()
wf.logger.debug('Analyzing Commonality')

for datafile, outputID, params in wf.single_datasets():
    #common_arrays = cmn_analysis.get_common_arrays(empirical.Patch(datafile), eval(params['grid']))
    #graph.common_plot(common_arrays, outputID, params, interactive=wf.runs.interactive)
    a_d = cmn_analysis.merge_common_arrays(empirical.Patch(datafile), eval(params['grid']))
    graph.common_ADsquared_plot(a_d, outputID, params, interactive=wf.runs.interactive)
    
wf.logger.info('Commonality analysis complete')


