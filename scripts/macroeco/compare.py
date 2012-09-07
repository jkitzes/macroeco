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

__author__ = "Mark Wilber"
__copyright__ = "Copyright 2012, Regents of University of California"
__credits__ = "John Harte"
__license__ = None
__version__ = "0.1"
__maintainer__ = "Mark Wilber"
__email__ = "mqw@berkeley.edu"
__status__ = "Development"

#NOTE: This needs to be updated when more distributions/sars are ready!
dist_dict  = {'logser' : logser(), 'plognorm' : plognorm(),
              'logser_ut' : logser_ut(), 'logser_ut_appx' : logser_ut_appx(),
              'plognorm_lt' : plognorm_lt(), 'sugihara': sugihara(),
              'broken_stick' : broken_stick(), 'geo_ser' : geo_ser(),
              'lognorm' : lognorm(), 'binm' : binm(), 'pois' : pois(), 'nbd' :
              nbd(), 'fnbd' : fnbd(), 'geo' : geo(), 'fgeo' : fgeo(), 'tgeo' :
              tgeo()} 

sar_dict = {'METE_sar' : METE_sar(), 'powerlaw' : powerlaw(), 'logser_ut_tgeo':
            logser_ut_tgeo(), 'logser_ut_fgeo' : logser_ut_fgeo(),
            'logser_ut_binm' : logser_ut_binm(), 'plognorm_lt_binm' :
            plognorm_lt_binm(), 'plognorm_lt_tgeo' : plognorm_lt_tgeo(),
            'plognorm_lt_fgeo' : plognorm_lt_fgeo()}


class CompareDistribution(object):
    '''
    Comparison object compares a list of data to any number of distributions

    '''
    
    #TODO: Error Checking
    def __init__(self, data_list, dist_list, clean=False, patch=False):
        '''
        Parameters
        ----------
        data_list : list or tuple output from Patch object
            List of np.arrays containing data
        dist_list : list
            List of distribution objects or list of distribution names
        patch : True
            If True, expects data_list to be input from the sad or ssad method
            in Patch.  If false, expects a list of iterables
        clean : bool
            If true, removes zeros from data_list.  Necessary for SAD
            comparisons.

        '''
        if patch:
            self.items = data_list[0]
            self.criteria = []
            self.data_list = []
            for obj in data_list[1]:
                self.criteria.append(obj[0])
                self.data_list.append(obj[1])
        else:
            self.data_list = [np.array(data) for data in data_list]
            self.criteria = None
        if clean:
            dlist = []
            for data in self.data_list:
                dlist.append(data[np.where((data != 0))[0]])
            self.data_list = dlist

        if np.all([type(dist) == str for dist in dist_list]):
            self.dist_list = np.empty(len(dist_list), dtype=object)
            for kw in list(dist_dict.viewkeys()):
                self.dist_list[np.where((kw == np.array(dist_list)))[0]] = \
                                dist_dict[kw]
            self.dist_list = list(self.dist_list)
        else:
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
            dist.fit(self.data_list)

            nlls = nll(dist.pmf(self.data_list)[0])
            #NOTE: dist.par_num is the number of parameters that
            k = np.repeat(dist.par_num, len(nlls))
            if crt:
                obs = np.array([len(data) for data in self.data_list])
                aic_vals.append(aicc(nlls, k, obs))
            else:
                aic_vals.append(aic(nlls, k))
        return list(np.array(aic_vals).T)

    def compare_aic_measures(self, crt=False):
        '''
        Compare AIC weights, delta_AIC, and AIC values across the different 
        models. Output is a three item tuple where each item is a list of 
        arrays with each array having length equal to the number of models 
        proposed and the length of the list is the lenth of self.data_lists.
        See Returns for tuple description.
        
        Parameters
        ----------
        crt : bool
            If True, calculates the corected AIC weights for the given data. 
            If False, calculates AIC weights.

        Returns
        -------
        : tuple
            first element is a list of arrays with each array having length 
            equal to the number of models proposed and the length of the list 
            is the lenth of self.data_lists.  Second element is are the delta
            AIC values in the same format as the first tuple object. The third
            object are the AIC values in the same the output of the 
            compare_aic method. 

        '''
        aic_vals = self.compare_aic(crt=crt)
        aic_wghts = []; delta_aic = []
        for mods_aic in aic_vals:
            taic_wghts, tdelta_aic = aic_weights(mods_aic)
            aic_wghts.append(taic_wghts)
            delta_aic.append(tdelta_aic)
        return aic_wghts, delta_aic, aic_vals
    
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
            #Different Identifier?
            rads_dict[get_name(dist)] = dist.rad()
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
            cdfs_dict[get_name(dist)] = dist.cdf(self.data_list)[0] 
        return cdfs_dict
    
#
    def compare_LRT(self, null_mdl):
        '''
        Performs a likelihood ratio test (LRT) on the distributions with in
        self.dist_list with the parameter nll_mdl as the null model. While this
        function will generate output on non-nested models, the models must be
        nested for the output to be meaningful.

        Parameters
        ----------
        null_mdl : distribution object
            The null distribution object to use in the LRT.

        Returns
        -------
        : dict
            A dictionary with keywords 'null_model, alternative model.' Each
            keyword references a list of length len(self.data_list) which
            contains tuples that contain the output of the function
            likelihood_ratio_test.  The LRT is performed on each data set in
            self.data_list for each given model pair.

        '''
        LRT_list = {}
        null_mdl.fit(self.data_list)
        null_nlls = nll(null_mdl.pmf(self.data_list)[0])
        for i, dist in enumerate(self.dist_list):
            dist.fit(self.data_list)
            alt_nlls = nll(dist.pmf(self.data_list)[0])
            k = dist.par_num - null_mdl.par_num
            df = np.repeat(k, len(alt_nlls))
            lrt = likelihood_ratio_test(null_nlls, alt_nlls, df)
            comp_kw = get_name(null_mdl) + ", " + get_name(dist)
            LRT_list[comp_kw] = lrt
        return LRT_list

    def summary(self, mins=10, crt=False):
        '''
        Summarizes the given datasets and the predicted rads. Looks at
        total balls sampled ('balls'), number of urns ('urns'), the max balls
        in a given urn ('max'), number of urns with less than MIN balls ('tot
        <= MIN'), and the fit of the the distributions in self.dist_list to the
        data in self.data_list

        Parameters
        ----------
        mins : int
            Bins with balls less than or equal to 10
        crt : bool
            If True, corrected AIC, if False, not

        Returns
        -------
        : dict
            Dictionary of dictionaries of length self.dist_list + 1.  Each
            sub-dictionary other than 'obs' contains the keywords balls, urns,
            max, tot <= MIN, and aic

        '''
        summary = {}
        rads = self.compare_rads()
        for kw in rads.iterkeys():
            summary[kw] = {}
            summary[kw]['balls'] = [np.sum(data) for data in rads[kw]]
            summary[kw]['urns'] = [len(data) for data in rads[kw]]
            summary[kw]['max'] = [np.max(data) for data in rads[kw]]
            summary[kw]['tot_min'] = [sum(data <= mins) for data in 
                                                                rads[kw]]
        aic_vals = self.compare_aic_measures(crt=crt)
        names = [get_name(dist) for dist in self.dist_list]
        for i, nm in enumerate(names):
            summary[nm]['aic'] = list(np.array(aic_vals[2]).T)[i]
            summary[nm]['aic_d'] = list(np.array(aic_vals[1]).T)[i]
            summary[nm]['aic_w'] = list(np.array(aic_vals[0]).T)[i]
        return summary

class CompareSARCurve(object):
    '''
    Object allows comparison between sar curve objects

    '''
    
    def __init__(self, sar_list, curve_list, full_sad):
        '''
        NOTE: Requiring the full_sad makes this function a lot less flexible.
        Maybe we shouldn't do that.  Maybe just give N and S.  Should we get
        the anchor area from the a_list?

        Parameters
        ----------
        sar_list : list of tuples or list of outputs from Patch().sar
            A list of tuples where each tuple contains two array-like objects
            of the same length.  The first element in the tuple is the
            area list and the second element is the species count for the sar.
            The maximum area in the area list should be the anchor area from
            which the full_sad was generated.
        curve_list : list
            A list of SARCurve objects or list of SARCurve object names (str)
        full_sad : list of array-like objects
            List of complete sads.  Each sad corresponds to an element in
            sar_list. 
        patch : bool
            If True, sar_list should be list of outputs from Patch().sar method
        '''
        assert len(sar_list) == len(full_sad), "sar_list and full_sad must " \
                                              + " be the same length"
        self.full_sad = [np.array(sad) for sad in full_sad]
        self.a_list = []
        self.sar_list = []
        for sar in sar_list:
            self.a_list.append(np.array(sar[0]))
            self.sar_list.append(np.array(sar[1]))

        if np.all([type(cur) == str for cur in curve_list]):
            self.curve_list = np.empty(len(curve_list), dtype=object)
            for kw in list(sar_dict.viewkeys()):
                self.curve_list[np.where((kw == np.array(curve_list)))[0]] = \
                                sar_dict[kw]
            self.curve_list = list(self.curve_list)
        else:
            self.curve_list = curve_list

        self.full_sad = [np.array(sad) for sad in full_sad]

    def compare_curves(self):
        '''
        Method generates predicted SAR curves from the given observed data and
        curve objects for comparison

        Returns
        -------
        : list of dicts
            The list is the same length self.sar_list and each dictionary is
            the length of self.curve_list + 1.  Each keyword in a dictionary
            references either the observed SAR ('obs') or the SAR generate by
            one of the curve objects. 
        '''
        pred_sar = []
        for sar, a, sad in zip(self.sar_list, self.a_list, self.full_sad):
            psar = {}
            psar['obs'] = np.array(zip(sar, a), dtype=[('species', np.float),
                                        ('area', np.float)])
            for cur in self.curve_list:
                cur.fit((a, sar), sad)
                psar[get_name(cur)] = cur.vals(a, np.max(a))
            for kw in psar.iterkeys():
                psar[kw].sort(order='area')
            pred_sar.append(psar)
        return pred_sar

def nll(pdist):
    '''
    Parameters
    ----------
    pdist : list of arrays
        List of pmf values on which to compute the negative log-likelihood

    Returns
    -------
    :list
        List of nll values

    '''
    return [-sum(np.log(dist)) for dist in pdist]

    

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

    emp_data = cnvrt_to_arrays(emp_data)[0]
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
    neg_L : array-like object
        The negative log likelihood of a model
    k : array-like object
        The number of parameters of a model
    
    
   Returns
   -------
   : float
        AIC for a given model
    '''
    neg_L, k = cnvrt_to_arrays(neg_L, k)
    assert len(k) == len(neg_L), "neg_L and k must have the same length"
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
    neg_L, k, n = cnvrt_to_arrays(neg_L, k, n)
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
    : tuple
        First element contains the relative AIC weights, second element
        contains the delta AIC values.

    Notes
    -----
    AIC weights can be interpreted as the probability that a given model is the
    best model in comparison to the other models

    '''
    aic_values = cnvrt_to_arrays(aic_values)[0]
    aic_values = np.array(aic_values)
    minimum = np.min(aic_values) 
    delta = np.array([x - minimum for x in aic_values])
    values = np.exp(-delta / 2)
    weights = np.array([x / sum(values) for x in values])
    return weights, delta

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
    data1, data2 = cnvrt_to_arrays(data1, data2) 
    data1 = np.array(data1)
    data2 = np.array(data2)
    return stats.ks_2samp(data1, data2)

def likelihood_ratio_test(nll_null, nll_alt, df_list):
    '''
    This functions compares of two nested models using the likelihood ratio
    test.

    Parameters
    ----------
    nll_null : array-like object
        The negative log-likelihood of the null model
    nll_alt : array-like object
        The negative log-likelihood of the alternative model
    df_list : array-like object
        the degrees of freedom calulated as (number of free parameters in
        alternative model) - (number of free parameters in null model)
    
    Returns
    -------
    : list of tuples
        (test_statistic, p-value)

    Notes
    -----
    The LRT only applies to nested models.  
    '''
    
    nll_null, nll_alt, df_list = cnvrt_to_arrays(nll_null, nll_alt, df_list)
    assert len(nll_null) == len(nll_alt) and len(nll_null) == len(df_list) and\
           len(nll_alt) == len(df_list), "nll_null, nll_alt, and df_list " + \
                                          "must have the same length"
    test_stat = 2 * nll_null - (2 * nll_alt)
    return [(ts, stats.chisqprob(ts, df)) for ts, df in zip(test_stat, df_list)]

def cnvrt_to_arrays(*args):
    '''
    Converts all args to np.arrays
    '''
    arg_list = []
    for arg in args:
        try:
            len(arg); arg = np.array(arg)
        except:
            arg = np.array([arg])
        arg_list.append(arg)
    return tuple(arg_list)

def get_name(obj):
    '''
    Return the name of the object
    '''
    return obj.__class__.__name__
