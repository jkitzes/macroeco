#!/usr/bin/python

'''
Predict SAR from a base scale SAD and a spatial abundance distribution.

Functions
---------
- `predict_sar` -- Predict SAR
- `universal_sar` -- Upscale and downscale Universal SAR
- `universal_sar_curve` -- Universal

References
----------

Harte, J., Smith, A. B., and Storch, D. 2009. Biodiversity scales from plots to biomes 
with a universal species-area curve. Ecology Letters, 12:789-797.

Harte, J. 2011. Maximum Entropy and Ecology: A Theory of Abundance,
Distribution, and Energetics. Oxford University Press.


'''

from __future__ import division
import numpy as np
import scipy.optimize
import sys

class DownscaleError(Exception):
    '''Catch downscale errors'''
    def __init__(self, value=None):
        Exception.__init__(self)
        self.value = value
    def __str__(self):
        return '%s' % self.value


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

def universal_sar(S, N, anchor_area, upscale=0, downscale=0):
    '''
    Predict the universal SAR curve for the given S and N found at 
    the given anchor scale

    Parameters
    ----------
    S : int
        Total number of species at the given anchor area
    N : int
        Total number of individuals at the given anchor area
    anchor_area : float
        The area from which the SAR will be upscaled or downscaled.
    upscale : int
        Number of iterations up from the anchor scale.  Each iteration doubles the
        previous area.
    downscale : int
        Number of iterations down from the anchor scale. Each iteration halves the 
        previous area.

    Returns
    -------
    : 1D structured np.array
        The structured array has fields 'species' and 'area'

    Notes
    -----
    This function uses equations 3, 7, 8, and 9 found in Harte et al. (2009).
    When possible, the approximation sum(x**n / n) ~= log(1 / log( 1/x)) was
    used to decrease runtime.  
    
    
    '''

    assert S < N, "S must be less than N"
    assert S > 1, "S must be greater than 1"
    assert N > 0, "S must be greater than 0"
    if upscale == 0 and downscale == 0:
        return np.array((S, anchor_area), dtype=[('species', np.float),\
                                                ('area', np.float)])
    areas = _generate_areas_(anchor_area, upscale, downscale)
    sar = np.empty(len(areas), dtype=[('species', np.float), ('area', np.float)])
    sar['area'] = areas
    if upscale != 0: 
        sar['species'][downscale:] = _upscale_sar_(areas[downscale:], N, S)
    if downscale != 0:
        sar['species'][:downscale + 1] = _downscale_sar_(areas[:downscale + 1], N, S)

    return sar

def _upscale_sar_(up_areas, N, S):
    '''
    This function is used to upscale from the anchor area.

    up_areas -- A 1D area of length upscale + 1.  The areas to be upscaled to.

    N -- Number of individuals at anchor scale (int)

    S -- Number of species at anchor scale (int)

    returns:
        1D array of species at a given upscaled area
    '''

    
    num_ind = np.empty(len(up_areas))
    spp = np.empty(len(up_areas))
    for i in xrange(len(up_areas)):
        if i == 0:
            num_ind[i] = N
            spp[i] = S
        else:
            num_ind[i] = 2 * num_ind[i - 1]
            N2A = num_ind[i]
            n = np.linspace(1, N2A, num=N2A)
            eq1 = lambda x: (sum((x**n)/n) * ((N2A) / ((x - x**(N2A + 1)) / \
                            (1 - x)) * (1 / x))) - ((N2A) * ((1 - x) / \
                            (x - x**(N2A + 1))) * (1 - (x**N2A / (N2A + 1)))) - spp[i - 1]
            x = scipy.optimize.brentq(eq1, 1e-10, 1 - 1e-10, disp=True)
            eq2 = lambda S2A: np.log(1 / np.log(1 / x)) * ((N2A) / ((x - x**(N2A + 1)) / (1 - x))) - S2A
            S2A = scipy.optimize.brentq(eq2, spp[i - 1], 2 * spp[i - 1], disp=True)
            spp[i] = S2A
    return spp

def _downscale_sar_(down_areas, N, S):
    '''
    This function is used to downscale from the anchor area.

    up_areas -- A 1D area of length downscale + 1.  The areas to be downscaled to.

    N -- Number of individuals at anchor scale (int)

    S -- Number of species at anchor scale (int)

    returns:
        1D array of species at a given downscaled areas
    '''

    num_ind = np.empty(len(down_areas))
    spp = np.empty(len(down_areas))
    for i in xrange(len(down_areas)):
        if i == 0:
            num_ind[i] = N
            spp[i] = S
        else:
            num_ind[i] = 0.5 * num_ind[i - 1]
            if num_ind[i] <= 1:
                raise DownscaleError('Cannot downscale %i iterations from anchor scale.\
                                     One or less individuals per cell.' % (len(down_areas) - 1))
            N_A = num_ind[i - 1]
            S_A = spp[i - 1]
            n = np.linspace(1, N_A, num=N_A)
            eq1 = lambda x: sum((x**n) / n) - ((S_A / N_A) * sum(x**n))
            x = scipy.optimize.brentq(eq1, 1e-10, min((sys.float_info[0] / S)**(1/float(N)),\
                              2), disp=True)
            ShalfA = (S_A * (1 / x)) - ((N_A) * ((1 - x) / (x - x**(N_A + 1))) * \
                     (1 - ((x**N_A) / (N_A + 1))))
            spp[i] = ShalfA
    return spp[::-1]



def _generate_areas_(anchor_area, upscale, downscale):
    '''
    Utility function that makes the area list
    
    '''

    areas = np.empty(upscale + downscale + 1)
    areas[downscale] = anchor_area
    for i in range(downscale)[::-1]:
        areas[i] = areas[i + 1] / 2
    for i in xrange(downscale + 1, len(areas)):
        areas[i] = areas[i - 1] * 2

    return areas

def universal_sar_curve(S=20, N=40, num_iter=10):
    '''
    universal_sar_curve(S=20, N=40, num_iter=10)

    Parameters
    ----------
    S : int
        Total number of species at the given anchor area
    N : int
        Total number of individuals at the given anchor area
    num_iter : int
        Number of upscaling iterations.
        WARNING: Running more than 10 ten begins to take along time
    
    Returns
    -------
    :1D structured np.array
        The structured array has fields 'z' and 'N/S' for plotting 
        convenience

    Notes
    -----
    This function calculates the universal SAR curve.  Different ratios
    of N/S will only put you at different places along the same curve. 
    Iterations of more than 15 take a long time. The equations used in this
    function were taken from Harte et al. (2009) and Harte (2011)

    '''

    num_ind = np.empty(num_iter + 1)
    spp = np.empty(num_iter + 1)
    univ_SAR = np.empty(num_iter + 1, dtype=[('z', np.float), ('N/S', np.float)])
    z = lambda beta: 1 / (np.log(2) * np.log(1 / beta)) #From Harte et al. (2009)

    for i in xrange(num_iter + 1):
        if i == 0:
            num_ind[i] = N
            spp[i] = S
            n = np.linspace(1, N, num=N)
            eq = lambda x: (((x - x ** (N + 1)) / ( 1 - x)) / sum((x ** n) / n)) - (N / S)
            x = scipy.optimize.brentq(eq, 1e-10, min((sys.float_info[0] / S)**(1/float(N)),\
                              2), disp=True)
            univ_SAR['z'][i] = z(-np.log(x))
            univ_SAR['N/S'][i] = N / S
        else:
            num_ind[i] = 2 * num_ind[i - 1]
            N2A = num_ind[i]
            n = np.linspace(1, N2A, num=N2A)
            eq1 = lambda x: (sum((x**n)/n) * ((N2A) / ((x - x**(N2A + 1)) / \
                            (1 - x)) * (1 / x))) - ((N2A) * ((1 - x) / \
                            (x - x**(N2A + 1))) * (1 - (x**N2A / (N2A + 1)))) - spp[i - 1]
            x = scipy.optimize.brentq(eq1, 1e-10, 1 - 1e-10, disp=True)
            univ_SAR['z'][i] = z(-np.log(x))
            #Using approximations for summations. See module notes
            eq2 = lambda S2A: np.log(1 / np.log(1 / x)) * ((N2A) / ((x - x**(N2A + 1)) / (1 - x))) - S2A
            S2A = scipy.optimize.brentq(eq2, spp[i - 1], 2 * spp[i - 1], disp=True)
            spp[i] = S2A
            univ_SAR['N/S'][i] = N2A / S2A
    return univ_SAR








