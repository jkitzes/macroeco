#!/usr/bin/python

'''Script that uses reproducible workflow to make and save graphs
of SADs 


Parameters in parameters.xls
----------------------------
Required:

name='grids' value=list of tuples

    The list of tuples specifies the divisions of the plot
    example:
    name='grids' value='[(16,16), (4,4), (2,2)]'

name='distr' value=distribution name
    
    The predicted distribution to be compared against the
    empirical SAD.  See predict_sad.py for options
    example:
    name='distr' value='mete' or name='distr' value='neg_binom'


Optional:

See graph.py


'''
#TODO: Need to descided how much flexiblity we want to give the user when graphing
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
wf.logger.debug('Analyzing SAD(s)')


for datafile, outputID, params in wf.single_datasets():
    sad = sad_analysis.get_gridded_sad_list(datafile, eval(params['grid']), clean=True)
    for i in xrange(len(sad)):
        for m in xrange(len(sad[i])):
            ID = outputID + '_grid_[' + str(i) + '][' + str(m) + ']'
            graph.sad_cdf_obs_pred_plot(sad[i][m], ID, params, interactive=wf.runs.interactive)
            graph.sad_rank_abund_plot(sad[i][m], ID, params, interactive=wf.runs.interactive)
            graph.obs_pred_rarity(sad[i][m], ID, params, interactive=wf.runs.interactive)
            graph.obs_pred_Nmax(sad[i][m], ID, params, interactive=wf.runs.interactive)

    

wf.logger.info("Analysis Complete")





    


