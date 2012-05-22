#!/usr/bin/python

'''Script that uses reproducible workflow to make and save graphs
of SADs 


Parameters in parameters.xls
----------------------------
Required:

name='grid' value=list of tuples

    The list of tuples specifies the divisions of the plot
    example:
    name='grid' value='[(16,16), (4,4), (2,2)]'

name='distr' value=distribution name
    
    The predicted distribution to be compared against the
    empirical SAD.  See predict_sad.py for options
    example:
    name='distr' value='mete' or name='distr' value='neg_binom'


Optional:

See graph.py

Output
------
Per data set per run, this analysis will output 4 plots and 2 csv's for each cell made by the 
parameter grids. If 'grids' is [(1,1)], the SAD for the entire plot is analyzed
and only 4 plots and 2 csv's are output.  If 'grid' is something like [(2,2), (4,4),
(6,6)], seperate plots are made for each tuple in the list.  Within each tuple (x,y) in
'grid' x*y sets of 4 plots and 2 csv's are made and output. In the output name, one
will see something that looks like *_grid_[i][j].  i represents the corresponding tuple
in 'grid' where j represents the cell number in the gridded plot, counting down columns
then across rows (where 0 is the first entry).  So for example, if 'grid' is 
[(2,2), (4,4)] the file *_grid_[1][13] corresponds to the tuple (4,4) and the 14th cell
in the grid counting down columns and across rows.  The 14th cell corresponds to cell (3,0).

The 4 plots consist of an observed vs. predicted cdf plot for the distribution specified in
parameters.xml, a rank-abundance plot, a bar graph of observed and predicted Nmax, and a
bar graph of observed and predicted rarity.  The 2 csv's are the data that led to the 
cdf plot and the rank-abundance plot.

'''

#TODO: Need to decide how much flexiblity we want to give the user when graphing
import sys, os
from macroeco.workflow import Workflow
from macroeco import predict_sad
from macroeco import empirical
from macroeco import output
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
            output.sad_cdf_obs_pred_plot(sad[i][m], ID, params, interactive=wf.runs.interactive)
            output.sad_rank_abund_plot(sad[i][m], ID, params, interactive=wf.runs.interactive)
            output.obs_pred_rarity(sad[i][m], ID, params, interactive=wf.runs.interactive)
            output.obs_pred_Nmax(sad[i][m], ID, params, interactive=wf.runs.interactive)

    

wf.logger.info("Analysis Complete")





    


