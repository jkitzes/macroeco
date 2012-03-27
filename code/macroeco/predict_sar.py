#!/usr/bin/python

'''
Predict SAR from a base scale SAD and a spatial abundance distribution.

Functions
---------
- `predict_sar` -- Predict SAR
'''

from __future__ import division
import numpy as np


def predict_sar(sad, S, a_list, ssad, k_eq):
    '''
    Predict the SAR for a given SAD, list of area fractions, and ssad.

    Parameters
    ----------
    sad : ndarray
        Species abundance distribution, should sum to 1 or nearly so. Support 
        must be >= 1 (ie, no P(0) at start)
    S : int or float
        Number of species in landscape
    a_list : list
        List of area fractions at which to calculate SAD
    ssad : function
        Spatial abundance distribution function from ssad module. May not be 
        tgeo, because cannot accept vectors for N and a.
    k_eq : function
        Function that calculates k as a function of n (mean abund per cell). If 
        ssad does not require k_array, k_eq should be something that evaluates 
        to False.

    Returns
    -------
    sar : ndarray
        Array of mean prediction of S found at each a in a_list
    '''
    # TODO: Expand function to do upscaling when a > 1

    sar = []
    size = sad.shape[0]
    N_range = np.arange(1, size + 1)

    # Loop through each area fraction
    for i, a in enumerate(a_list):
        assert a < 1, "a must be < 1"

        if not k_eq:
            p_pres = 1 - ssad(0, N_range, np.repeat(a, size), summary = False)
        else:
            k_array = k_eq(N_range * a)
            p_pres = 1 - ssad(0, N_range, np.repeat(a, size), k_array, summary 
                              = False)

        sar.append(sum(S * sad * p_pres))
    
    return np.array(sar)
