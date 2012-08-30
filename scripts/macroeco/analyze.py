#!/usr/bin/python

'''
This module contains functions for analyzing empirical and theoretical
distributions.  


Functions
---------
-`empirical_cdf` -- Empirical cdf for given data
-`aic` -- Calculate AIC value
-`aicc` -- Calculate corectted AIC value
-`aic_wieghts` -- Calculate AIC weights for models
-`ks_two_sample_test` -- Kolmogrov-Smirnov two sample test
-`likelihood_ratio_test` -- Calculated likelihood ratio for nested models


'''

from __future__ import division
import numpy as np
import scipy.stats as stats
from distributions import *

import scipy.optimize 
import scipy.special
import math as m
import scipy.integrate as integrate
import sys


__author__ = "Mark Wilber"
__copyright__ = "Copyright 2012, Regents of University of California"
__credits__ = "John Harte"
__license__ = None
__version__ = "0.1"
__maintainer__ = "Mark Wilber"
__email__ = "mqw@berkeley.edu"
__status__ = "Development"

class Compare_Distr(object):
    '''
    Comparison object allows a list of data to any number of distributions


    Examples
    --------

    '''

    def __init__(self, data_list, dist_list):
        '''
        Parameters
        ----------
        data_list : list
            List of np.arrays containing data
        dist_list : list
            List of distribution objects that are already

        '''
        self.data_list = [np.array(data) for data in data_list]
        self.dist_list = dist_list

    def compare_aic(self, crt=False):
        '''
        Get the aic or aicc values for every data set and for every
        distribution

        Parameters
        ----------
        crt : bool
            If True, calculates the corected AIC for the given data. If False,
            calculates AIC.

        Returns
        -------
        : list
            A list of arrays.  The list has length = to number of data sets in
            self.data_list.  Each array within list has the length of
            self.dist_list.

        '''
        aic_vals = []
        for dist in self.dist_list:
            #NOTE: Is distribution object already fit? Let's say not
            #NOTE: What about 'a' parameter?  Does the object already have
            #it?
            dist.fit(self.data_list)
            nlls = dist.nll(self.data_list)
            #NOTE: dist.par_num is the number of parameters that
            #distribution has. Need to make this an attribute
            k = np.repeat(dist.par_num, len(nlls))
            if crt:
                obs = np.array([len(data) for data in self.data_list])
                aic_vals.append(aicc(nlls, k, obs))
            else:
                aic_vals.append(aic(nlls, k))
        return list(np.array(aic_vals).T)

    def compare_aic_weights(self, crt=False):
        '''
        Compare AIC weights across the different models. Output is a list of
        arrays with each array having length equal to the number of models
        proposed and the length of the list is the lenth of self.data_lists
        
        Parameters
        ----------
        crt : bool
            If True, calculates the corected AIC weights for the given data. 
            If False, calculates AIC weights.

        Returns
        -------
        : list
            a list ofarrays with each array having length equal to the number
            of models proposed and the length of the list is the lenth of 
            self.data_lists

        '''
        aic_vals = self.compare_aic(crt=crt)
        return [aic_weights(mods_aic) for mods_aic in aic_vals]
    
    #Maybe return a dict instead?
    def compare_rads(self):
        '''
        Compares rank abundance distributions for all data in data_list and to
        the given distributions

        Returns
        -------
        : dict
            Has len(self.dist_list) + 1.  All the distribution class names
            passed to the constructor are key words as well 'obs' which
            references the observed data, self.data_list. Each keyword looks up
            a list of arrays.  Each list is len(self.data_list) long and
            contains the predicted rads for the empirical data sets for the
            given distribution. 

        '''
        rads_dict = {}
        rads_dict['obs'] = self.data_list
        for i, dist in enumerate(self.dist_list):
            dist.fit(self.data_list)
            #NOTE: Need to make sure this is implemented
            #Different Identifier?
            rads_dict[dist.__str__().split('.')[1].split(' ')[0] + str(i)] = dist.rad()
        return rads_dict

    def compare_cdfs(self):
        '''
        Compares cdfs for all data in data_lists and to the empirical cdfs

        Returns
        -------
        :dict
            Has len(self.dist_list) + 1.  All the distribution class names
            passed to the constructor are key words as well 'obs' which
            references the observed data, self.data_list. Each keyword looks up
            a list of arrays.  Each list is len(self.data_list) long and
            contains the predicted cdfs for the empirical data sets for the
            given distribution. 


        '''

        cdfs_dict = {}
        cdfs_dict['obs'] = [empirical_cdf(data) for data in self.data_list]
        for i, dist in enumerate(self.dist_list):
            dist.fit(self.data_list)
            #Might need to reference dist differently...
            cdfs_dict[dist.__str__().split('.')[1].split(' ')[0] + str(i)] =\
                        dist.cdf(self.data_list) 
        return cdfs_dict


def empirical_cdf(emp_data):
    '''
    Generates an empirical cdf from a empirical data

    Parameters
    ----------
    emp_data : array-like object
        Empirical data 

    Returns
    --------
    :ndarray
        An empirical cdf
    '''

    try:
        len(emp_data); emp_data = np.array(emp_data)
    except:
        emp_data = np.array([emp_data])
    unq_vals = np.unique(emp_data)
    leng = len(emp_data)
    cdf = np.empty(len(emp_data))
    count = 0
    for i in unq_vals:
        loc = np.where((i == emp_data))[0]
        count += len(loc)
        cdf[loc] = count / leng
    return cdf

def aic(neg_L, k):
    '''
    Calculates the AIC of a given model

    Parameters
    ----------
    neg_L : float
        The negative log likelihood of a model
    k : float
        The number of parameters of a model
    
    
   Returns
   -------
   : float
        AIC for a given model
    '''
    try:
        len(neg_L); neg_L = np.array(neg_L)
    except:
        neg_L = np.array([neg_L])
    try:
        len(k); k = np.array(k)
    except:
        k = np.array([k])
    assert len(k) == len(neg_L)
    aic = (2 * neg_L) + (2 * k)
    return aic

def aicc(neg_L, k, n):
    '''
    Calculates the corrected AIC of a given model

    Parameters
    ----------
    neg_L : array-like object
        The negative log likelihood of models
    k : array-like object
        The number of parameters of models
    n : array-like object
        Number of observations for each model

    Returns
    -------
    : nparray
        AICc for a given models

    '''
    try:
        len(n); n = np.array(n)
    except:
        n = np.array([n])
    try:
        len(k); k = np.array(k)
    except:
        k = np.array([k])
    try:
        len(neg_L); neg_L = np.array(neg_L)
    except:
        neg_L = np.array([neg_L])
    assert len(neg_L) == len(k) and len(neg_L) == len(n) and len(k) == len(n),\
            "neg_L, k, and n must all have the same length"
    aic_value = aic(neg_L, k)
    return aic_value + ((2 * k * (k + 1)) / (n - k - 1))

def aic_weights(aic_values):
    '''
    Calculates the aic_weights for a given set of models

    Parameters
    ----------
    aic_values : array-like object
        Array-like object containing AIC values from different models
    
    Returns
    -------
    : np.ndarray
        Array containing the relative AIC weights

    Notes
    -----
    AIC weights can be interpreted as the probability that a given model is the
    best model in comparison to the other models

    '''

    #NOTE: Check length of aic_values
    if type(aic_values) == float or type(aic_values) == int:
        raise ValueError("Parameter must be array-like object")
    aic_values = np.array(aic_values)
    minimum = np.min(aic_values) 
    delta = np.array([x - minimum for x in aic_values])
    values = np.exp(-delta / 2)
    weights = np.array([x / sum(values) for x in values])
    return weights

def ks_two_sample(data1, data2):
    '''Function uses the Kolomogrov-Smirnov two-sample test to determine if the
    two samples come from the same distribution.  Note that the KS-test is only
    valid for continuous distributions

    Parameters
    ----------
    data1 : array-like object
        Array-like object which contains a set of data to compare
    data2 : array-like object
        Array-like object which contains a set of data to compare

    Returns
    -------
    : tuple
        (D-statistic, two-sided p-value)

    '''
    
    data1 = np.array(data1)
    data2 = np.array(data2)
    return stats.ks_2samp(data1, data2)

def likelihood_ratio_test(neg_LL_null, neg_LL_alt, df):
    '''
    This functions compares of two nested models using the likelihood ratio
    test.

    Parameters
    ----------
    neg_LL_null : float
        The negative log-likelihood of the null model
    neg_LL_alt : float
        The negative log-likelihood of the alternative model
    df : int
        the degrees of freedom calulated as (number of free parameters in
        alternative model) - (number of free parameters in null model)
    
    Returns
    -------
    : tuple
        (test_statistic, p-value)

    Notes
    -----
    The LRT only applies to nested models.  
    '''

    test_stat = 2 * neg_LL_null - (2 * neg_LL_alt)
    return (test_stat, stats.chisqprob(test_stat, df))


##TEST CLASSES, REMOVE BEFORE RELEASE##
class RootError(Exception):
    '''Error if no root or multiple roots exist for the equation generated
    for specified values of S and N'''

    def __init__(self, value=None):
        Exception.__init__(self)
        self.value = value
    def __str__(self):
        return '%s' % self.value

class DownscaleError(Exception):
    '''Catch downscale errors'''
    def __init__(self, value=None):
        Exception.__init__(self)
        self.value = value
    def __str__(self):
        return '%s' % self.value

class Distribution(object):

    def __init__(self, **kwargs):
        '''
        Generic constructor

        **kwargs : keyword parameters for distribution

        '''
        self.params = kwargs
    
    #Better option than lower bound?
    def cdf(self, n_list, lower_bound=1):
        '''
        Cumulative distribution function.  Determined by summing pmf

        Parameters
        ----------
        n : int, float or array-like object
            Values at which to calculate the cdf

        Returns
        -------
        : ndarray (1D)
            Returns array with cdf values for the given values of n.

        '''
        S_list = self.params.get('S', None)
        N_list = self.params.get('N', None) 
        assert S_list != None, "S parameter not given"
        assert N_list != None, "N parameter not given"
        cdf_list = []
        for n, S, N in zip(n_list, S_list, N_list):
            try:
                len(n); n = np.array(n)
            except:
                n = np.array([n])
            max_n = np.max(n)
            self.params['S'] = [S]
            self.params['N'] = [N]
            cdf = np.cumsum(self.pmf([np.arange(lower_bound, max_n + 1)]))
            cdf_list.append(np.array([cdf[x - lower_bound] for x in n]))
        self.params['S'] = S_list
        self.params['N'] = N_list
        return cdf_list

    def rad(self):
        '''
        Rank abundance distribution. Calculated using pmf

        Parameters
        ----------
        None

        Returns
        -------
        : ndarray (1D)
            Returns and array of length S with the expected abundances given
            pmf

        '''
        S_list = self.params.get('S', None)
        N_list = self.params.get('N', None) 
        assert S_list != None, "S parameter not given"
        assert N_list != None, "N parameter not given"
        rad_list = []
        for S, N in zip(S_list, N_list):
            self.params['S'] = [S]
            self.params['N'] = [N]
            full_pmf = self.pmf([np.arange(1, N + 1)])[0]
            rad_list.append(make_rank_abund(full_pmf, S))
        self.params['S'] = S_list
        self.params['N'] = N_list
        return rad_list

    def fit(self, data_list, sad=True):
        '''
        Generic fit to sad or ssads

        Parameters
        ----------
        data : list of np.arrays objects
            Data used to calculate fit
        sad : bool
            If True, calculate S and N parameter, else only calculate N.

        Return
        ------
        None
        '''
        #Might not be the best check here
        if sad:
            self.params['N'] = [np.sum(data) for data in data_list]
            self.params['S'] = [len(data) for data in data_list]
            return self
        else:
            self.params['N'] = [len(data) for data in data_list]
            return self

    def nll(self, n_list):
        '''
        Calcuates the negative log-likelihood for a given distribution

        Parameters
        ----------
        n : array-like object
            Values for which to calculate nll

        Returns
        : float
            Negative log-likelihood for given values

        '''
        pmfs = self.pmf(n_list)
        return [-sum(np.log(pmf)) for pmf in pmfs]

class lgsr(Distribution):
    '''
    Fisher's log series distribution (Fisher et al. 1943, Hubbel 2001).

    Methods
    -------
    pmf(n, param_out); S and N parameters passed into __init__
        Probability mass function
    cdf(n); S and N parameters passed into __init__
        Cumulative distribution function
    rad(); S and N parameters passed into __init__
        Rank abundance distribution

    Parameters
    ----------
    S : int
        Total number of species in landscape
    N : int
        Total number of indviduals in landscape

    Notes
    -----
    S and N are passed into the constructor (__init__) as keyword arguments.
    Example: lgsr(S=34, N=345)

    '''
    def __init__(self, **kwargs):
        self.params = kwargs
        self.par_num = 1
        
    def pmf(self, n_list, param_out=False):
        '''
        Probability mass function of Fisher log series generated with S and N.

        Parameters
        ----------
        n : int, float or array-like object
            Abundances at which to calculate the pmf
        
        Passed into __init__
        ---------------------
        S : int
            Total number of species in landscape
        N : int
            Total number of indviduals in landscape

        Returns
        -------
        : ndarray (1D)
            Returns array with pmf for values for the given values n. If 
            param_out = True, returns the array as well as the parameter 
            estimates.

        Notes
        -----
        Multiplying the pmf by S yields the predicted number of species
        with a given abundance.

        '''
        S_list = self.params.get('S', None)
        N_list = self.params.get('N', None) 
        assert S_list != None, "S parameter not given"
        assert N_list != None, "N parameter not given"
        assert len(n_list) == len(S_list), "n_list must be the same length as S"
        pmf_list = []
        for n, S, N in zip(n_list, S_list, N_list):
            assert S < N, "S must be less than N"
            assert S > 1, "S must be greater than 1"
            assert N > 0, "N must be greater than 0"
            try:
                len(n); n = np.array(n)
            except:
                n = np.array([n])
            assert np.max(n) <= N, "Maximum n cannot be greater than N"

            start = -2
            stop = 1 - 1e-10
    
            eq = lambda x, S, N: (((N/x) - N) * (-(np.log(1 - x)))) - S
    
            x = scipy.optimize.brentq(eq, start, stop, args=(S,N), disp=True)
            pmf = stats.logser.pmf(n, x)

    
            if param_out == True:
                pmf_list.append((pmf, {'x' : x}))
            else:
                pmf_list.append(pmf)
        return pmf_list



