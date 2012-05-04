#!/usr/bin/python

'''
Contains functions to help analyze a given sad

MORE COMMENTING LATER
'''

from __future__ import division
from macroeco import predict_sad
from macroeco import empirical
import numpy as np


__author__ = "Mark Wilber"
__copyright__ = "Copyright 2012, Regents of University of California"
__credits__ = "John Harte"
__license__ = None
__version__ = "0.1"
__maintainer__ = "Mark Wilber"
__email__ = "mqw@berkeley.edu"
__status__ = "Development"

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
        an array containing the cdf for the observed sad
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
    if len(cdf) != len(sad): 
            raise ValueError("Lengths don't match")
    return np.array(cdf)

def get_obs_vs_pred_cdf(sad, distr):
    '''Generates a structured array with n, observed cdf, and predicted cdf
    from the observed sad

     Parameters
    ----------
    sad : 1D np.array
        an array containing an SAD
    distr : string
        The predicted distribution:
        'mete' - METE
        'mete_approx' - METE with approximation
        'plognorm' - Poisson lognormal
        'trun_plognorm' - Truncated poisson lognormal
        'neg_binom' - Negative binomial
        'lgsr' - Fisher's log series

    Returns
    -------
    : Structured np.array, dtype=[('n', np.int), ('obs', np.float),
    ('pred', np.float)]
        Length of the returned structured array is the same as sad

    '''


    obs_cdf = get_obs_cdf_values(sad)
    sad_sorted = np.sort(sad)
    pred_cdf = predict_sad.get_sad_cdf(len(sad), np.sum(sad), distr, sad=sad)
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

def get_gridded_sad_list(datapath, grid, clean=False):
    '''Function takes in a datapath and returns a gridded sad
    
    Parameters
    ----------
    datapath : string
        Path to a properly formatted data file
        
    grid : list
        A list of tuples containing the desired gridding scale i.e [(1,1), (2,2)]
        See Patch.sad_grid() in empirical.py for more information
    
    clean : bool
        If False return SAD as is.  If True, remove zeros from SAD.

    Returns
    -------
    : list
        Length of list is equal to length of grid.  Each element in the list is a
        list of arrays.

    '''
    data = empirical.Patch(datapath)
    sad = data.sad_grid(grid)

    if clean:
        sad_clean = []
        cuts = []
        for cut in sad:
            for cell in cut:
                cell = cell[cell != 0]
                cuts.append(cell)
            sad_clean.append(cuts)
        return sad_clean

    else:
        return sad

def aic(neg_L, k, n, corrected=True):
    '''
    Calculates the AIC of a given model

    Parameters
    ----------
    neg_L : float
        The negative log likelihood of a model
    k : float
        The number of parameters of a model
    n : int
        Number of observations
    corrected : bool
        If True, returns the corrected AIC. If False, returns
        uncorrected AIC
    
    
   Returns
   -------
   : float
        AICc or AIC for a given model
    '''

    aic = (2 * neg_L) + (2 * k)
    if corrected:
        aicc = aic + ((2 * k * (k + 1)) / (n - k - 1))
        return aicc
    else:
        return aic


