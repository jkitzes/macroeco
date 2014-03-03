"""
===========================
Main (:mod:`macroeco.main`)
===========================

This module contains functions that compare the goodness of fit of a
distribution/curve to data or the fit of two distributions/curves to each 
other.

.. autosummary::
   :toctree: generated/

   main

"""

from __future__ import division

import numpy as np
import scipy.stats as stats

from distributions import *


def get_AIC(values, params):
    """
    Calculate AIC given values of a pdf/pmf and a set of model parameters.
    """
    k = len(params)  # Num parameters
    L = get_nll(values)
    return 2*k + 2*L

def get_nll(values):
    """
    Calculate negative log likelihood from an array of pdf/pmf values.
    """
    return -np.sum(np.log(values))

def get_empirical_cdf(data):

    min, max = 0, np.ceil(np.max(data))
    x = np.arange(min, max+2) # x max is 1 above emp_result max
    counts, _ = np.histogram(data, bins=x, normed=True)
    emp_cdf = np.cumsum(counts)

    return x[:-1], emp_cdf
