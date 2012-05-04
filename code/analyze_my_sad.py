#!/usr/bin/python

'''Script that uses reproducible workflow to make and save graphs
of SADs of data sets given in the command line

This analysis does not grid the SAD, it just returns the full SAD
from the plot and tests it against the METE SAD'''

import sys, os
from macroeco.workflow import Workflow
from macroeco import predict_sad
from macroeco import empirical
from macroeco import graph
from macroeco import sad_analysis

__author__ = "Mark Wilber"
__copyright__ = "Copyright 2012, Regents of the University of California"
__license__ = None
__version__ = "0.5"
__maintainer__ = "Mark Wilber"
__email__ = "mqw@berkeley.edu"
__status__ = "Development"

wf = Workflow()
wf.logger.debug('Only analyzing the SAD from the entire plot. Grid is fixed at [(1,1)]')


for datafile, outputID, params in wf.single_datasets():
    sad = sad_analysis.get_gridded_sad_list(datafile, [(1,1)], clean=True)[0][0]
    graph.sad_cdf_obs_pred_plot(sad, outputID, params, interactive=wf.runs.interactive)

    

wf.logger.info("Analysis Complete")





    


