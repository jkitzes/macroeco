#!/usr/bin/python

'''This module provides functions for graphing results'''
from __future__ import division
import matplotlib.pyplot as plt
import numpy as np
import logging
from macroeco.workflow import Workflow, loggername
from macroeco import predict_sad

__author__ = "Mark Wilber"
__copyright__ = "Copyright 2012, Regents of the University of California"
__license__ = None
__version__ = "0.5"
__maintainer__ = "Mark Wilber"
__email__ = "mqw@berkeley.edu"
__status__ = "Development"

module_logger = logging.getLogger(loggername)

def make_obs_vs_pred_cum_plot(sad, outputID, params):

    cdf = get_obs_vs_pred_cdf(sad)
    jitt_x, jitt_y = jitter(cdf['pred'], cdf['obs'], jitter_scale=0.007)
    plt.plot([0] + list(cdf['obs']),[0] +  list(cdf['obs']),\
                                            color='gray', linestyle='--')
        
    plt.scatter(jitt_x, jitt_y, color='grey') 
    plt.scatter(cdf['pred'], cdf['obs'], color='black')
    for s in xrange(len(cdf)):
        plt.text(cdf['pred'][s], cdf['obs'][s] + 0.01, "n=" + str(cdf['n'][s]),\
                                fontsize=7)
    plt.ylabel('Observed cdf')
    plt.xlabel('Predicted (METE) cdf')
    plt.title('Observed vs. predicted values of the cdf')
    plt.xlim(-0.1,1.1)
    plt.ylim(-0.1,1.1)
    plt.legend(('Slope 1', 'All species cum.\nprob. (jittered)',\
                    'Cum. prob.\nfor given n'),loc='best', prop={'size':12})
    plt.savefig(outputID)
    plt.clf()

def jitter(x, y, jitter_scale=0.001):
    '''Jitter the points'''
    assert len(x) == len(y), "x and y must be the same length"
    jitt_x = np.copy(x)
    jitt_y = np.copy(y)
    xy_tuple = []
    for i in xrange(len(x)):
        xy_tuple.append((x[i], y[i]))

    xy_tuple = np.array(xy_tuple)
    unq_xy = np.unique(xy_tuple)
    for xy in unq_xy:
        indices = np.where(xy_tuple == xy)[0]
        if len(indices > 1):
            for ind in xrange(len(indices)):
                if ind > 0:
                    jitt_x[indices[ind]] += (jitter_scale * ind)
    return jitt_x, jitt_y

def get_obs_cdf_values(sad):
    '''Generates predicted cdf values from an observed
    SAD.
    
    Parameters
    ----------
    sad : 1D np.array
        array containing an SAD

    Returns
    -------
    : 1D np.array
        an array comtaining the cdf for the observed sad
    '''
    unq_sad = np.unique(sad)
    S = len(sad)
    cdf = []
    count = 0
    for i in unq_sad:
        tot_in = sum((i == sad))
        count += tot_in
        #Removing or adding (1/(2*S)) can change the graphs
        for s in xrange(tot_in):
            cdf.append((count / S))# - (1/(2*S)))
    assert len(cdf) == len(sad), "Lengths don't match"
    return np.array(cdf)

def get_obs_vs_pred_cdf(sad):
    '''Generates a structured array with n, observed cdf, and predicted cdf
    from the observed sad

     Parameters
    ----------
    sad : 1D np.array
        an array containing an SAD

    Returns
    -------
    : Structured np.array, dtype=[('n', np.int), ('obs', np.float),
    ('pred', np.float)]
        Length of the returned structured array is the same as sad

    '''


    obs_cdf = get_obs_cdf_values(sad)
    sad_sorted = np.sort(sad)
    pred_cdf = predict_sad.mete_lgsr_cdf(len(sad), np.sum(sad))
    cpred_cdf = []
    for n in sad_sorted:
        cpred_cdf.append(pred_cdf['cdf'][n - 1])
    cpred_cdf = np.array(cpred_cdf)
    obs_vs_pred = np.empty(len(sad_sorted), dtype=[('n', np.int), ('obs', np.float),\
                                               ('pred', np.float)])
    obs_vs_pred['n'] = sad_sorted
    obs_vs_pred['obs'] = obs_cdf
    obs_vs_pred['pred'] = cpred_cdf
    return obs_vs_pred
