#!/usr/bin/python

'''Script that uses reproducible workflow to test commonality predictions
on data sets. Specifically, this script looks at Commonality as a function
of Distance ** 2 / Area

Parameters in parameters.xls
----------------------------
Required:

name='grids' value=list of tuples

    The list of tuples specifies the divisions of the plot
    example:
    name='grids' value='[(16,16), (4,4), (2,2)]'

Optional:

None
    

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

#TODO: Fail gracefully if there is no params['grid'], it might do this already
for datafile, outputID, params in wf.single_datasets():
    data = empirical.Patch(datafile)
    cmn = cmn_analysis.get_common_arrays(data, eval(params['grid']))
    graph.common_ADsquared_plot(cmn, outputID, params, interactive=wf.runs.interactive)
    
wf.logger.info('Commonality analysis complete')


