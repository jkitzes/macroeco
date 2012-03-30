#!/usr/bin/python

'''Annotated, skeleton script showing how to run a particular
analysis as a reproducible workflow.

This is also the setup expected by the GUI.'''

import sys, os
from matplotlib.mlab import csv2rec
import matplotlib.pyplot as plt
from macroeco.workflow import Workflow

 
__author__ = "Chloe Lewis"
__copyright__ = "Copyright 2012, Regents of the University of California"
__license__ = None
__version__ = "0.5"
__maintainer__ = "Chloe Lewis"
__email__ = "chlewis@berkeley.edu"
__status__ = "Development"

def single_patch(data, outputID):
    '''Dummied-up activity on a single dataset'''
    # This is some dummied-up activity on dummied-up data
    fig= plt.figure()
    ax = fig.add_subplot(111)

    # Column headers are accessible along with the data
    x = data.dtype.names[1]
    y = data.dtype.names[2]
    phenom = data.dtype.names[0]


    # Here we use the names of our data columns, and a parameter
    ax.scatter(data[x], data[y], s=data[phenom]*10, c = thisrun['color'])
    ax.set_xlabel(x)
    ax.set_ylabel(y)
    ax.set_title(phenom)

    # Save results so the user doesn't need to remember to do so
    fig.savefig(outputID + '_' + runID + '.png')

    # Parameters keeps track if we're running interactively or automatically;
    #    don't use show() if this is expected to run overnight, etc.
    if wf.runs.interactive:
        plt.show()
    else:
        wf.logger.debug('%s: Non-interactive run, suppressing plot.show()'%outputID)



def multiple_patch(data, outputID):
    '''Dummied-up activity on multiple datasets'''
    fig = plt.figure()
    rangeX = fig.add_subplot(121)
    rangeY = fig.add_subplot(122)
    N = range(len(data.keys()))
    xmins = []
    xmaxs = []
    ymins = []
    ymaxs = []
    for dname in data.keys():
        x = data[dname].dtype.names[1]
        y = data[dname].dtype.names[2]
        xmins.append(min(data[dname][x]))
        xmaxs.append(max(data[dname][x]))
        ymins.append(min(data[dname][y]))
        ymaxs.append(max(data[dname][y]))
    #    don't use show() if this is expected to run overnight, etc.

    rangeX.hlines(N, xmins, xmaxs, linewidth = 20, colors=thisrun['color'])
    rangeX.set_yticks([])
    
    rangeY.vlines(N, ymins, ymaxs, linewidth = 20, colors=thisrun['color'])
    rangeY.set_xticks([])

    rangeX.set_title('X-extents of datasets')
    rangeY.set_title('Y-extents of datasets')
    
    # Save results so the user doesn't need to remember to do so
    fig.savefig(outputID + '_' + runID + '.png')

    # Parameters keeps track if we're running interactively or automatically;
    #    don't use show() if this is expected to run overnight, etc.
    if wf.runs.interactive:
        plt.show()
    else:
        wf.logger.debug('%s: Non-interactive run, suppressing plot.show()'%outputID)

# The rest of this is the workflow organization of the script. Should be packaged into a class
# (github issue 76) TODO

wf = Workflow({'color':'Plotting color'})


# If there is more than one set of parameters stored, run them all: 
for runID in wf.runs.params.keys():
    thisrun = wf.runs.params[runID]
    wf.logger.info('Running %s as a single-dataset analysis with parameters %s'%(wf.script, runID))

    for dataset in wf.data.keys():
        # Concatenate the arguments into a unique-enough identifier for this analysis.
        outputID = '_'.join([wf.script, dataset])
        single_patch(wf.data[dataset], outputID)
                     

    if len(sys.argv) > 2:
        wf.logger.info('Running %s as a multiple-dataset analysis with parameters %s'%(wf.script, runID))

        outputID = '_'.join([wf.script] + wf.data.keys())
        multiple_patch(wf.data, outputID)

        
            
    



