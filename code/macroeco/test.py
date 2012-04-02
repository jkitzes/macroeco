#!/usr/bin/python

'''Annotated, skeleton script showing how to run a particular
analysis as a reproducible workflow.

This is also the setup expected by the GUI.'''

import sys, os
import matplotlib.pyplot as plt
import logging
from matplotlib.mlab import csv2rec
from macroeco.workflow import Workflow, loggername


 
__author__ = "Chloe Lewis"
__copyright__ = "Copyright 2012, Regents of the University of California"
__license__ = None
__version__ = "0.5"
__maintainer__ = "Chloe Lewis"
__email__ = "chlewis@berkeley.edu"
__status__ = "Development"

module_logger = logging.getLogger(loggername)

def single_patch(datapath, outputID, parameters, interactive=False):
    '''Dummied-up activity to use parameters and logging on a single dataset.'''
    # This is some dummied-up activity on dummied-up data
    module_logger.debug(str(parameters))
    fig= plt.figure()
    ax = fig.add_subplot(111)

    # Column headers are accessible along with the data
    data = csv2rec(datapath)
    x = data.dtype.names[1]
    y = data.dtype.names[2]
    phenom = data.dtype.names[0]


    # Here we use the names of our data columns, and a parameter
    ax.scatter(data[x], data[y], s=data[phenom]*10, c = parameters['color'])
    ax.set_xlabel(x)
    ax.set_ylabel(y)
    ax.set_title(phenom)

    # Save results so the user doesn't need to remember to do so
    fig.savefig(outputID + '.png')

    # need to suppress show() if this is expected to run overnight, on a server, etc.
    if interactive:
        plt.show()
    else:
        module_logger.debug('%s: Non-interactive run, suppressing plot.show()'%outputID)



def multiple_patch(datapaths, outputID, parameters):
    '''Dummied-up activity on multiple datasets.
    Just as a demo, *always* runs non-interactively, so *does not need* an interactive keyword.'''
    fig = plt.figure()
    rangeX = fig.add_subplot(121)
    rangeY = fig.add_subplot(122)
    N = range(len(datapaths))
    xmins = []
    xmaxs = []
    ymins = []
    ymaxs = []
    for dname in datapaths:
        data = csv2rec(dname)
        x = data.dtype.names[1]
        y = data.dtype.names[2]
        xmins.append(min(data[x]))
        xmaxs.append(max(data[x]))
        ymins.append(min(data[y]))
        ymaxs.append(max(data[y]))

    rangeX.hlines(N, xmins, xmaxs, linewidth = 20, colors=parameters['color'])
    rangeX.set_yticks([])
    
    rangeY.vlines(N, ymins, ymaxs, linewidth = 20, colors=parameters['color'])
    rangeY.set_xticks([])

    rangeX.set_title('X-extents of datasets')
    rangeY.set_title('Y-extents of datasets')
    
    # Save results so the user doesn't need to remember to do so
    fig.savefig(outputID + '.png')
    module_logger.debug('%s: Non-interactive function.'%outputID)
