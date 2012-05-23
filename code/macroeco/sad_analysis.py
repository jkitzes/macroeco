#!/usr/bin/python

'''
Contains functions to help analyze a given sad

This module is the interface between the theoretical functions
found in predict_sad and the empirical sad's generated by 
empirical

'''

from __future__ import division
from macroeco import predict_sad
from macroeco import empirical
import numpy as np
import scipy.stats as stats


__author__ = "Mark Wilber"
__copyright__ = "Copyright 2012, Regents of University of California"
__credits__ = "John Harte"
__license__ = None
__version__ = "0.1"
__maintainer__ = "Mark Wilber"
__email__ = "mqw@berkeley.edu"
__status__ = "Development"

def get_obs_cdf_values(sad):
    '''Generates cdf values from an observed
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
    sad = np.array(sad)
    if len(np.where(sad == 0)[0]) != 0:
        raise ValueError("SAD cannot contain zeros")

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
        'geo' - Geometric
        'lgsr' - Fisher's log series

    Returns
    -------
    : Structured np.array, dtype=[('n', np.int), ('obs', np.float),
    ('pred', np.float)]
        Length of the returned structured array is the same as sad. 'n' is the number of
        individuals within a species.

    '''

    sad = np.array(sad)
    sad_sorted = np.sort(sad)
    obs_cdf = get_obs_cdf_values(sad_sorted)
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
        Length of list is equal to length of grid.  Each element in the list is an
        array of arrays.

    '''
    data = empirical.Patch(datapath)
    sad = data.sad_grid(grid)

    if clean:
        sad_clean = []
        for cut in sad:
            cuts = []
            for cell in cut:
                clean_cell = cell[cell != 0]
                cuts.append(clean_cell)
            sad_clean.append(cuts)
        return sad_clean

    else:
        return sad

def get_obs_and_pred_rarity(sad, distr, n=10):
    '''Generates the number of observed and predicted rare species.
    Rarity is defined as the number of species with abundance less than n.

    Parameters
    ----------
    sad : 1D np.array
        an array containing an SAD

    distr : string
        The predicted SAD distribution:
        'mete' - METE
        'mete_approx' - METE with approximation
        'plognorm' - Poisson lognormal
        'trun_plognorm' - Truncated poisson lognormal
        'neg_binom' - Negative binomial
        'geo' - Geometric
        'lgsr' - Fisher's log series

    n : int
        species are considered rare if they have abundances less than n

    Returns
    -------
    : tuple
    Tuple of (observed, predicted) rare species

    '''
    sad = np.array(sad)
    if len(np.where(sad == 0)[0]) != 0:
        raise ValueError("SAD cannot contain zeros")
    obs_rare = np.sum(sad < n)
    pred_abund = predict_sad.make_rank_abund(predict_sad.macroeco_pmf(len(sad),\
                                        np.sum(sad), distr, sad=sad), len(sad))
    pred_rare = np.sum(pred_abund < n)
    return (obs_rare, pred_rare)

def get_obs_and_pred_Nmax(sad, distr):
    '''Gets the Nmax for observed and predicted sads. 

    Parameters
    ----------
    sad : 1D np.array
        an array containing an SAD
    distr : string
        The predicted SAD distribution:
        'mete' - METE
        'mete_approx' - METE with approximation
        'plognorm' - Poisson lognormal
        'trun_plognorm' - Truncated poisson lognormal
        'neg_binom' - Negative binomial
        'geo' = Geometric
        'lgsr' - Fisher's log series

    Returns
    -------
    : tuple
    Tuple of (observed, predicted) Nmax

    '''
    sad = np.array(sad)
    if len(np.where(sad == 0)[0]) != 0:
        raise ValueError("SAD cannot contain zeros")
    obs_Nmax = np.max(sad)
    pred_abund = predict_sad.make_rank_abund(predict_sad.macroeco_pmf(len(sad),\
                                        np.sum(sad), distr, sad=sad), len(sad))
    pred_Nmax = np.max(pred_abund)
    return (obs_Nmax, pred_Nmax)

def get_obs_pred_abund(sad, distr):
    '''Gets the observed and predicted abundances for the sad

     Parameters
    ----------
    sad : 1D np.array
        an array containing an SAD
    distr : string
        The predicted SAD distribution:
        'mete' - METE
        'mete_approx' - METE with approximation
        'plognorm' - Poisson lognormal
        'trun_plognorm' - Truncated poisson lognormal
        'neg_binom' - Negative binomial
        'geo' - Geometric
        'lgsr' - Fisher's log series

    Returns
    -------
    :structured array
        Returns a structured array with dtype=[('observed', np.int), ('predicted', np.int)]


    '''
    sad = np.array(sad)
    if len(np.where(sad == 0)[0]) != 0:
        raise ValueError("SAD cannot contain zeros")

    obs_abund = np.copy(sad)
    obs_abund.sort()
    pred_pmf = predict_sad.macroeco_pmf(len(sad), sum(sad), distr,  sad=sad)
    pred_abund = predict_sad.make_rank_abund(pred_pmf, len(sad))
    obs_pred_abund = np.empty(len(sad), dtype=[('observed', np.int), ('predicted', np.int)])
    obs_pred_abund['observed'] = obs_abund
    obs_pred_abund['predicted'] = pred_abund
    return obs_pred_abund




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

def lognorm_KS(sad):
    '''Function uses the Kolomogrov-Smirnov test to determine if an
    empirical sad is lognormal

    Parameters
    ----------
    sad : array like object
        The empirical SAD

    Returns
    -------
    : tuple
        (D-statistic, two-sided p-value)

    '''
    
    sad = np.array(sad)
    logsad = np.log(sad)
    norm_logsad = (logsad - logsad.mean()) / logsad.std()
    return stats.kstest(norm_logsad, 'norm')

def get_values_for_sad(sad, distr):
    '''Function takes in sad and returns state variables
    and other information based on specified distribution

    Parameters
    ----------
    sad : array like object
        The empirical SAD
    
    distr : string
        The distribution to which to compare the empirical sad:
        'mete' - METE
        'mete_approx' - METE with approximation
        'plognorm' - Poisson lognormal
        'trun_plognorm' - Truncated poisson lognormal
        'neg_binom' - Negative binomial
        'geo' - Geometric
        'lgsr' - Fisher's log series
    
    Returns
    -------
    : dict
        A dictionary containing state variables and other requested values
        'S' - Total species in SAD
        'N' - Total individuals in SAD
        'Nmax_obs' - Number of individuals in most abundant species
        'rarity_obs' - Number of species with n < 10
        'distr' - dictionary of distribution specific values
            'name' - name of distribution (see above)
            'nll' - negative-log likelihood
            'AICc' - Corrected AIC
            'Nmax_pred' - Predicted Nmax
            'rarity_pred' = Predicted rarity


    '''
    assert type(sad) == tuple or type(sad) == list or type(sad) == np.ndarray,\
                'SAD must be tuple, list, or ndarray'

    sad = np.array(sad)
    if len(np.where(sad == 0)[0]) != 0:
        raise ValueError("SAD cannot contain zeros")
    value_dict = {}
    value_dict['S'] = len(sad)
    value_dict['N'] = sum(sad)
    Nmax = get_obs_and_pred_Nmax(sad, distr)
    rarity = get_obs_and_pred_rarity(sad, distr)
    value_dict['Nmax_obs'] = Nmax[0]
    value_dict['rarity_obs'] = rarity[0]
    value_dict['distr'] = {}
    value_dict['distr']['name'] = distr
    value_dict['distr']['nll'] = predict_sad.nll(sad, distr)
    value_dict['distr']['params'] = predict_sad.distr_parameters(len(sad), sum(sad),\
                                        distr, sad=sad)
    value_dict['distr']['AICc'] = aic(value_dict['distr']['nll'],\
                                    len(value_dict['distr']['params']), len(sad))
    value_dict['distr']['Nmax_pred'] = Nmax[1]
    value_dict['distr']['rarity_pred'] = rarity[1]

    return value_dict











