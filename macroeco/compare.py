"""
===========================
Compare (:mod:`macroeco.compare`)
===========================

This module contains functions that compare the goodness of fit of a
distribution/curve to data or the fit of two distributions/curves to each 
other.

Comparison Functions
====================

.. autosummary::
   :toctree: generated/

   get_AIC
   get_AICC
   get_AIC_weights
   get_nll
   get_empirical_cdf
   get_sum_of_squares


"""

from __future__ import division

import numpy as np
import scipy.stats as stats
import pandas as pd

from distributions import *

# NOTE: get_* functions usually refer to a method within a class.  I would
# suggest dropping all of the get prefixes


def get_AIC(values, params):
    """
    Calculate AIC given values of a pdf/pmf and a set of model parameters.
    """
    k = len(params)  # Num parameters
    L = get_nll(values)
    return 2*k + 2*L

def get_AICC(values, params):
    """
    Calculate AICC given values of a pdf/pmf and a set of model parameters.

    Notes
    -----
    Should be used when the number of observations is < 40.

    References
    ----------
    .. [#]
        Burnham, K and Anderson, D. (2002) Model Selection and Multimodel
        Inference: A Practical and Information-Theoretic Approach (p. 66). New
        York City, USA: Springer.

    """

    k = len(params)  # Num parameters
    n = len(values)  # Num observations
    return get_AIC(values, params) + (2*k * (k + 1)) / (n - k - 1)

def get_AIC_weights(aic_values):
    """
    Calculates the aic_weights for a given set of models

    Parameters
    ----------
    aic_values : array-like object
        Array-like object containing AIC values from different models
    
    Returns
    -------
    (weights, delta) : tuple
        First element contains the relative AIC weights, second element
        contains the delta AIC values.

    Notes
    -----
    AIC weights can be interpreted as the probability that a given model is the
    best model in comparison to the other models
    """

    aic_values = _to_arrays(aic_values)[0]
    minimum = np.min(aic_values) 
    delta = aic_values - minimum
    values = np.exp(-delta / 2)
    weights = values / np.sum(values)

    return weights, delta

def get_nll(values):
    """
    Calculate negative log likelihood from an array of pdf/pmf values.
    """
    return -np.sum(np.log(values))

def get_empirical_cdf(data):
    """
    Generates an empirical cdf from empirical data

    Parameters
    ----------
    data : array-like object
        Empirical data 

    Returns
    --------
    : array
        The empirical cdf corresponding to the inputted data

    """

    vals = pd.Series(data).value_counts()
    ecdf = pd.DataFrame(data).set_index(keys=0)
    probs = pd.DataFrame(vals.sort_index().cumsum() / np.float(len(data)))
    ecdf = ecdf.join(probs)

    return np.array(ecdf[0])

class _gen_loss_function(object):
    """
    Generic class for loss function between observed and predicted data

    """

    def __init__(self, loss_fxn_str):
        """
        Parameters
        ----------
        loss_fxn_str : string
            A Python string representing the loss function between observed
            (obs) and predicted (pred).

        Ex. 'np.abs(obs - pred)' or (obs - pred)**2
        """
        self.loss_fxn = loss_fxn_str

    def total_loss(self, obs, pred):
        """
        Total loss for observed and predicted 

        Parameters
        ----------
        obs, pred : array-like objects
            observed and predicted data

        Returns
        -------
        : float
            The sum of the loss function
        """

        obs, pred = _to_arrays(obs, pred)
        return np.sum(eval(self.loss_fxn))

get_sum_of_squares = _gen_loss_function('(obs - pred)**2').total_loss
        
def _to_arrays(*args):
    '''
    Converts all args to np.arrays
    '''
    return tuple([np.array(ta) for ta in args])    
