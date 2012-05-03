#!/usr/bin/python

'''Script that uses reproducible workflow to make and save graphs
of SADs of data sets given in the command line'''

import sys, os
from macroeco.workflow import Workflow
from macroeco import predict_sad
from macroeco import empirical
from macroeco import graph

__author__ = "Mark Wilber"
__copyright__ = "Copyright 2012, Regents of the University of California"
__license__ = None
__version__ = "0.5"
__maintainer__ = "Mark Wilber"
__email__ = "mqw@berkeley.edu"
__status__ = "Development"

wf = Workflow()

for datafile, outputID, params in wf.single_datasets():
    data = empirical.Patch(datafile)
    sad = data.sad_grid([(1,1)])[0][0]
    sad = sad[sad != 0]
    graph.make_obs_vs_pred_cum_plot(sad, outputID, params)


    


