#!/usr/bin/python

'''Script that uses reproducible workflow to test commonality predictions
on data sets. 

Parameters in parameters.xls
----------------------------
Required:

name='grid' value=list of tuples

    The list of tuples specifies the divisions of the plot
    example:
    name='grid' value='[(16,16), (4,4), (2,2)]'

    NOTE: ONE MUST INCLUDE quotation marks around value list.

Optional:

None

Output
------
Per data set, per run, this analysis will output 2 plots and csv files corresponding
to the number of tuples in the parameter value 'grid'. Each tuple in 'grid' corresponds
to a different cell area.  All of these cell areas are plotted together with commonality
on the y and distance on the x.  For example, if 'grid' is "[(16,16), (4,4), (2,2)]", there
will be three csv files, each one containing commonality data for a single tuple.  There will
be one Commonality by Distance plot with three lines, each representing a different cell area.

Notes
-----
A division more than (50, 50) takes many hours compute.

    
'''

import sys, os
from macroeco.workflow import Workflow
from macroeco import empirical
from macroeco import output
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
    output.common_plot_save(cmn, outputID, params, areas=[], interactive=wf.runs.interactive)
    
wf.logger.info('Commonality analysis complete')


