#!/usr/bin/python

''' This module contains classes and functions for comparing empirical and
predicted macroecological metrics. 

Classes
-------
CompareDistribution : Base class to for CompareSAD, CompareSSAD, CompareIED,
CompareSED, CompareASED

CompareSAD : Compares predicted species abundance distributions (SAD) with
empirical SADs

CompareSSAD : Compares predicted species-level spatial abundance distributions
(SSAD) with empirical SSADS

CompareSAR : Compares predicted species-area relationship (SAR) curves with
empirical SAR curves

CompareIED : Compares predicted individual energy distributions (IED) with
empirical IEDs  

CompareSED : Compares predicted species-level energy distributions (SED) with
empirical SEDs 

CompareASED : Compares predicted average species-level energy distributions
(ASED) with empirical ASEDs.

Functions
---------
-`empirical_cdf` -- Empirical cdf for given data
-`aic` -- Calculate AIC value
-`aicc` -- Calculate corectted AIC value
-`aic_wieghts` -- Calculate AIC weights for models
-`ks_two_sample_test` -- Kolmogrov-Smirnov two sample test
-`likelihood_ratio` -- Calculated likelihood ratio for nested models
-`variance` -- Calculates the variance for given datasets
-`skew` -- Calculates the skew for given datasets
-`kurtosis` -- Calculates the kurtosis for given data sets
-`bootstrap` -- Get bootstrapped samples from a dataset
-`'mean_squared_error` -- Calculates the MSE between an obs and pred data set


'''

from __future__ import division
import numpy as np
import scipy.stats as stats
from distributions import *
import copy
import random
import time
import logging


class CompareDistribution(object):
    '''
    Comparison object compares a list of data to any number of distributions

    '''
    
    #TODO: Error Checking
    def __init__(self, data_list, dist_list, observed_index):
        '''
        Parameters
        ----------
        data_list : list of iterables or list of tuples of iterables 
            data_list is any list of iterables or list of tuples of iterables
            that will be passed to the fit functions of the distribution
            objects in dist_list. data_list will be passed to fit functions for
            each distribution.  data_list undergoes no validation in __init__
        dist_list : list
            List of distribution objects or strings that have the same name as 
            a distribution object. If they are strings, they will be evaled
        observed_index : int
            The index of the desired observed metric in the tuples within
            data_list.  If 0, data_list can be a list of data
            rather than a list of tuples of data.  The index specified by
            object_ind will be considered the observed data.
        
        Notes
        -----
        All distribution objects are fit in the __init__ method.

        '''

        self.dist_list = make_dist_list(dist_list)

        # Fit the distributions objects 
        [dist.fit(data_list) for dist in self.dist_list]
        
        # Set the observed data
        if observed_index == 0 and np.all([type(dt) != type((1,)) for dt in
                                                            data_list]):
            self.observed_data = [np.array(dt) for dt in data_list]
        elif np.all([type(dt) == type((1,)) for dt in data_list]):
            self.observed_data = [np.array(dt[observed_index]) for dt in
                                                                     data_list]
        else:
            self.observed_data = [np.array(dt) for dt in data_list]

        # Set this in __init__ so other methods can check if compare_rads() has
        # been called
        self.rads = None
        self.cdfs = None
       
        # If attributes have not been instantiated, set to None
        try:
            self.sad_spp_list
        except:
            self.sad_spp_list = None
        try:
            self.criteria
        except:
            self.criteria = None

    def compare_mse(self, mse_base='cdf'):
        '''
        This function compares the mean squared error (mse) for each distribution
        against the observed data, self.observed_data.  Perfect predicted data
        would yield a mse of 0.  The lower the mse the better the predicted
        values fit the data. If mse_base='cdf' the mse is calculated from the
        cdf. If mse_base='rad', the mse is calculated from the rank_abundance
        distribution.

        Parameters
        -----------
        mse_base : str
           Either 'cdf' or 'rad'.  If 'cdf' the mse values are computed
           from the cumulative density function.  It 'rad' the mse values are
           computed from the rank abundance distribution. Default is 'cdf'

        Returns
        -------
        : dict
            A dictionary of length self.dist_list with keywords being the
            distribution names.  Each keyword looks up a list of length
            self.observed_data in which are the mse values comparing that
            distribution's predicted values (cdf or rad) to the corresponding
            observed values.

        Notes
        -----
        Calculating the mse from the cdf is the least bias approximater

        '''
        if mse_base == 'cdf':
            if self.cdfs == None:
                vals = self.compare_cdfs()
            else:
                vals = self.cdfs
        elif mse_base == 'rad':
            if self.rads == None:
                vals = self.compare_rads()
            else:
                vals = self.rads
        else:
            raise NameError('%s value for mse_base not recognized' % mse_base)


        mse = {}
        for kw in vals.iterkeys():
            if kw != 'observed':
                if not np.all([len(j) == 0 for j in vals[kw]]):
                    mse[kw] = [mean_squared_error(vals['observed'][i], 
                                   vals[kw][i]) for i in xrange(len(vals[kw]))]
                else:
                    logging.warning('MSE values for %s set to NaN' % kw)
                    mse[kw] = [np.NaN for i in xrange(len(self.observed_data))]
        return mse


    def compare_aic(self, crt=False):
        '''
        Get the aic or aicc values for every data set and for every
        distribution

        Parameters
        ----------
        crt : bool
            If True, calculates the corrected AIC for the given data. If False,
            calculates AIC.

        Returns
        -------
        : list
            A list of arrays.  The list has length = to number of data sets in
            self.observed_data.  Each array within list has the length of
            self.dist_list. The first element of the array corresponds to the
            first distribution in dist_list, the second corresponds to the
            second distribution, etc.

        '''
        aic_vals = []
        for dist in self.dist_list:
            
            try:
                nlls = nll(dist.pmf(self.observed_data))
            except NotImplementedError:
                try:
                    nlls = nll(dist.pdf(self.observed_data))
                except NotImplementedError:
                    logging.warning('%s has neither a PMF nor a PDF. AIC set'
                                            % get_name(dist) + ' to infinity')
                    nlls = np.repeat(np.inf, len(self.observed_data)) 
                    
            #NOTE: dist.par_num is the number of parameters of distribution
            k = np.repeat(dist.par_num, len(nlls))
            if crt:
                obs = np.array([len(data) for data in self.observed_data])
                aic_vals.append(aicc(nlls, k, obs))
            else:
                aic_vals.append(aic(nlls, k))
        return list(np.array(aic_vals).T)

    def compare_aic_measures(self, crt=False):
        '''
        Compare AIC weights, delta_AIC, and AIC values across the different 
        models. Output is a three item tuple where each item is a list of 
        arrays with each array having length equal to the number of models 
        proposed and the length of the list is the length of self.observed_data.
        See Returns for tuple description.
        
        Parameters
        ----------
        crt : bool
            If True, calculates the corrected AIC weights for the given data. 
            If False, calculates AIC weights.

        Returns
        -------
        : tuple
            The first element is a list of arrays with each array having length
            equal to the number of models proposed and the length of the list
            is the length of self.observed_data. The first element contains
            the AIC weights. The second element is the delta AIC values in the
            same format as the first tuple object. The third object are the AIC
            values in the same format as the output of the compare_aic method. 

        Notes
        -----
        The given AIC values in each array correspond to the distributions in
        self.dist_list. 

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
            passed to the constructor are key words as well as 'observed' which
            references the observed data, self.observed_data. Each keyword looks up
            a list of arrays.  Each list is len(self.observed_data) long and
            contains the predicted rads for the empirical data sets for the
            given distribution.

        Note
        ----
        If self.rads has already been set in another method (i.e. is not None).
        This method will not overwrite it.  To reset self.rads, set self.rads
        = None and then run self.compare_rads().

        '''
        if self.rads == None:
            rads_dict = {}
            rads_dict['observed'] = copy.deepcopy(self.observed_data)
            for i, dist in enumerate(self.dist_list):
                #Different Identifier?
                rads_dict[get_name(dist)] = dist.rad()

            self.rads = rads_dict
        return self.rads

    def compare_cdfs(self):
        '''
        Compares cdfs for all data in data_lists and to the empirical cdfs

        Returns
        -------
        :dict
            Has len(self.dist_list) + 1.  All the distribution class names
            passed to the constructor are key words as well 'observed' which
            references the observed data, self.observed_data. Each keyword looks up
            a list of arrays.  Each list is len(self.observed_data) long and
            contains the predicted cdfs for the empirical data sets for the
            given distribution. 


        '''
        if self.cdfs == None:

            cdfs_dict = {}
            cdfs_dict['observed'] = [empirical_cdf(data) for data in
                                                            self.observed_data]
            for i, dist in enumerate(self.dist_list):
                try:
                    cdfs_dict[get_name(dist)] = dist.cdf(self.observed_data)
                except NotImplementedError:
                    logging.warning('CDF method not implemented for %s' %
                                                                get_name(dist))
                    cdfs_dict[get_name(dist)] = [np.array([]) for i in
                                               xrange(len(self.observed_data))]
                    
            self.cdfs = cdfs_dict
        return self.cdfs
    

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
            keyword references a list of length len(self.observed_data) which
            contains tuples that contain the output of the function
            likelihood_ratio (chisquared, p-value).  The LRT is performed on
            each data set in self.observed_data for each given model pair.

        '''
        LRT_list = {}
        null_mdl.fit(self.observed_data)

        try:
            null_nlls = nll(null_mdl.pmf(self.observed_data))
        except:
            null_nlls = nll(null_mdl.pdf(self.observed_data))
        for i, dist in enumerate(self.dist_list):
            
            try:
                alt_nlls = nll(dist.pmf(self.observed_data))
            except:
                alt_nlls = nll(dist.pdf(self.observed_data))

            k = dist.par_num - null_mdl.par_num
            df = np.repeat(k, len(alt_nlls))
            lrt = likelihood_ratio(null_nlls, alt_nlls, df)
            comp_kw = get_name(null_mdl) + ", " + get_name(dist)
            LRT_list[comp_kw] = lrt
        return LRT_list

    def compare_rarity(self, mins_list):
        '''
        This method takes in the output from self.compare_rads and a list of
        minimum values against which to compare the observed and predicted
        rads.  and outputs a dictionary with length self.dist_list + 1 (all
        distributions + observed).  Each keyword in this dict looks up a dict
        of len(mins_list) where the keywords are the values against which the
        rads will be <=.  Each one of these sub-dictionaries looks up a list
        with len(self.observed_data).

        Parameters
        ----------
        mins_list : array-like object
            A list of numbers.  Each number number will be used in the
            following function: rad <= mins_list[i].

        Returns
        -------
        : dict
            Returns a dictionary with length self.dist_list + 1 (all
        distributions + observed).  Each keyword in this dict looks up a dict
        of len(mins_list) where the keywords are the values against which the
        rads will be <=.  Each one of these sub-dictionaries looks up a list
        with len(self.observed_data).


        '''

        # Don't remake rads if they have already been made
        if self.rads == None:
            rads = self.compare_rads()
        else:
            rads = self.rads

        mins_list = make_array(mins_list)

        rarity = {}
        keys = list(rads.viewkeys())
        for kw in keys:
            rarity[kw] = {}
            for mins in mins_list:
                rarity[kw][mins] = [sum(data <= mins) for data in rads[kw]]
        return rarity

    def compare_moments(self):
        '''
        Compare the higher order moments (variance, skew, kurtosis) for the
        given distributions and observed data. 

        Returns
        -------
        : dict
            A dictionary with keywords variance, skew, and kurtosis.  Each
            keyword looks up a dictionary len(dist_list) + 1 keywords.  The
            keywords are 'observed' and the distribution object names. Each of
            these keywords looks up a list of floats with the same length as
            data_list.

        '''

        if self.rads == None:
            rads = self.compare_rads()
        else:
            rads = self.rads

        var = {}
        skw = {}
        kurt = {}

        for kw in rads.iterkeys():
            var[kw] = variance(rads[kw])
            skw[kw] = skew(rads[kw])
            kurt[kw] = kurtosis(rads[kw])
        moments = {}
        moments['variance'] = var
        moments['skew'] = skw
        moments['kurtosis'] = kurt

        return moments

    def summary(self, mins_list=[10], crt=False):
        '''
        Summarizes the given datasets and the predicted rads. Looks at
        total balls sampled ('balls'), number of urns ('urns'), the max balls
        in a given urn ('max'), number of urns with less than MIN balls ('tot
        <= MIN'), and the fit of the distributions in self.dist_list to the
        data in self.observed_data

        'balls' is the sum of the observed data.  For a Species Abundance
        Distribution 'balls' would represent individuals.  For an Individual
        Energy Distribution 'balls' would represent energy.
        
        'urns' is the length of the observed data. For a Species Abundance
        Distribution 'urns' would represent species and for a Individual Energy
        Distribution 'urns' would represent individuals.

        Parameters
        ----------
        mins_list : list
            Bins with balls less than or equal to 10
        crt : bool
            If True, corrected AIC, if False, not

        Returns
        -------
        : dict
            Dictionary of dictionaries of length self.dist_list + 1.  Each
            sub-dictionary other than 'observed' contains the keywords balls,
            urns, max, tot_min, aic, aic_d, aic_w, and par_num.  Each of these
            keywords contains a list that is the same length as the number of
            sads under consideration.


            urns = total number of items in self.observed_data.  Could be
                   species (SAD, ASED), cells (SSAD), or individuals (IED, SED)
            balls = Items that are placed in urns. Could be individuals (SAD,
                    SSAD), energy (ASED, IED, SED).
            max = Maximum number of balls in an urn
            tot_min = Total number of urns with with <= a given number of balls
            aic = AIC
            aic_d = Delta AIC
            aic_w = AIC weights
            par_num = Parameter number of the given distribution
            tot_min = total counts less than or equal numbers in min_list
            vars = Additional variables computed for the given distribution


        '''
        summary = {}

        # Check that rads is already set, if not set it
        if self.rads == None:
            rads = self.compare_rads()
            if type(rads) == type((1,)):
                rads = rads[0]
        else:
            rads = self.rads
        
        rarity = self.compare_rarity(mins_list=mins_list)
        for kw in rads.iterkeys():
            summary[kw] = {}
            summary[kw]['balls'] = [np.sum(data) for data in rads[kw]]
            summary[kw]['urns'] = [len(data) for data in rads[kw]]
            summary[kw]['max'] = [np.max(data) for data in rads[kw]]
            summary[kw]['tot_min'] = rarity[kw]

        aic_vals = self.compare_aic_measures(crt=crt)
        names = [get_name(dist) for dist in self.dist_list]
        for i, nm in enumerate(names):
            summary[nm]['aic'] = list(np.array(aic_vals[2]).T)[i]
            summary[nm]['aic_d'] = list(np.array(aic_vals[1]).T)[i]
            summary[nm]['aic_w'] = list(np.array(aic_vals[0]).T)[i]
            summary[nm]['par_num'] = np.repeat(self.dist_list[i].par_num,
                                        len(list(np.array(aic_vals[2]).T)[i]))
            summary[nm]['vars'] = self.dist_list[i].var

        return summary

class CompareSAD(CompareDistribution):
    '''
    Object inherits CompareDistribution and uses it to compare species
    abundance distributions (SAD)

    Attributes
    ----------
    self.observed_data : A list of arrays
        Each array in this list is an SAD.  Each of these SADs will be compared
        to the distributions in self.dist_list
    self.dist_list : a list of distribution objects
        Each object is a distribution object to which the SADs in
        self.observed_data will be compared.  
    self.criteria : a list of dictionaries or None
        If not None, each dictionary specifies the divisions made on the plot
        that generated each SAD in self.observed_data.  self.criteria should be
        the same length as self.observed_data
    self.sad_spp_list : list of arrays or None
        If not None, each array contains the species strings for the
        corresponding SAD in self.observed_data.  The length of
        self.sad_spp_list should be the same length as self.observed_data and
        the length of any array within self.sad_spp_list should be the same
        length the corresponding array in self.observed_data. The index of any
        species name within any array within self.sad_spp_list references the
        species count with the same index in self.observed_data.

    '''
    
    def __init__(self, data_list, dist_list, patch=False):
        '''
        Parameters
        ----------
        data_list : list of iterables or output from Patch.sad
            List of np.arrays containing data
        dist_list : list
            List of distribution objects or strings that have the same name as 
            a distribution object. If they are strings, they will be evaled 
        patch : bool
            If True, expects the output from the Patch.sad method and if False, 
            expects a list of iterables. Presumably, each iterable is an SAD.

        Notes
        -----
        If data_list is a list of tuples containing iterables, the 1st entry
        (0th element) in each tuple is considered the observed SADs
        '''
        if patch == True:
            self.criteria, sad_data, self.sad_spp_list = unpack(data_list)
            super(CompareSAD, self).__init__(sad_data, dist_list, 0) 
        else:
            super(CompareSAD, self).__init__(data_list, dist_list, 0)

class CompareSSAD(CompareDistribution):
    '''
    Object inherits CompareDistribution and uses it to compare species-level
    spatial abundance distributions (SSAD)

    Attributes
    ----------
    self.observed_data : A list of arrays
        Each array in this list is an SSAD.  Each of these SSADs will be
        compared to the distributions in dist_list
    self.dist_list : a list of distribution objects
        Each object is a distribution object to which the SSADs in
        self.observed_data will be compared.  
    self.criteria : a list of dictionaries or None
        If not None, each dictionary specifies the divisions made on the plot
        that generated each SAD in self.observed_data.  self.criteria should be
        the same length as self.observed_data
    self.sad_spp_list : List of strings or None
        If not None, self.sad_spp_list is a list of strings where each string
        refers to a species.  The length of self.sad_spp_list should be the same
        length as self.observed_data.  Each species string has the same index
        within the list as its corresponding SSAD in self.observed_data. 

    '''
    
    def __init__(self, data_list, dist_list, patch=False):
        '''
        Parameters
        ----------
        data_list : list of iterables or output from Patch.ssad
            List of np.arrays containing data
        dist_list : list
            List of distribution objects or strings that have the same name as 
            a distribution object. If they are strings, they will be evaled 
        patch : bool
            If True, expects the output from the Patch.sad method and if False, 
            expects a list of iterables. Presumably, each iterable is an SSAD.


        Notes
        -----
        If data_list is a list of tuples containing iterables, the 1st entry
        (0th element) in each tuple is considered the observed SSADs
        '''
        if patch == True:

            self.sad_spp_list = list(data_list[1].viewkeys())
            ssad_data = [np.array(data_list[1][nm]) for nm in
                                                            self.sad_spp_list]
            self.criteria = data_list[0]

            super(CompareSSAD, self).__init__(ssad_data, dist_list, 0) 
        else:
            super(CompareSSAD, self).__init__(data_list, dist_list, 0)



class CompareIED(CompareDistribution):
    '''
    Class compares predicted individual energy distributions (IED) for the
    entire community to observed IEDs

    Attributes
    ----------
    self.observed_data : list of arrays
        Observed individual energy distributions (IED)
    self.ied_spp_lists : list of arrays
        Each array contains species strings which pair to the values
        contained in the corresponding array in self.ied_list. The length of
        self.ied_spp_lists should be the same length as self.ied_list.
    self.sad_spp_list : list of arrays
        If not None, each array contains the species strings for the
        corresponding SAD in self.sad_list.  The length of self.sad_spp_list
        should be the same length as self.sad_list and the length of any array
        within self.sad_spp_list should be the same length the corresponding
        array in self.sad_list. The index of any species name within any array
        within self.sad_spp_list references the species count with the same
        index in self.sad_list.
    self.criteria : a list of dictionaries or None
        If not None, each dictionary specifies the divisions made on the plot
        that generated each SAD and IED in self.sad_list and self.ied_list.
        self.criteria should be the same length as self.sad_list and
        self.ied_list. 
    self.dist_list : a list of distribution objects
        Each object is a distribution to which the IEDs in self.ied_list will
        be compared.

    '''

    def __init__(self, data_list, dist_list, patch=False):
        '''
        Parameters
        ----------
        data_list : list of tuples or output from Patch object
            A list containing tuples of length two.  The first object in a
            tuple an iterable containing the community individual energy
            distribution.  The second object in a tuple is an iterable
            containing the empirical species abundance distribution.
            See patch argument for more information.
        dist_list : list of strings or objects
            Each string corresponds to a name of a psi distribution to which to
            compare to the observed data. 
        patch: bool
            If True, expects a tuple of length 2 with the first object being
            the output from Patch.ied and the second element being the
            output from Patch.sad. If False expects what argument data_list
            describes. sads and energy should be made with the same criteria.

        Notes
        -----
        The __init__ method always removes zeros from the SADs
        
        If data_list is a list of tuples containing iterables, the 1st entry
        (0th element) in each tuple is considered the observed IEDs 
        '''

        if patch == True:
            # Unpack sad. Store spp_lists in items
            sad_criteria, sad_list, self.sad_spp_list = \
                                                    unpack(data_list[1])

            # Unpack ied
            ied_criteria, ied_list, self.ied_spp_lists = \
                                                        unpack(data_list[0])
            self.criteria = sad_criteria

            super(CompareIED, self).__init__(zip(ied_list, sad_list),
                                                                dist_list, 0)
            
        else:
            super(CompareIED, self).__init__(data_list, dist_list, 0)
            self.ied_spp_lists = None
    


class CompareSED(CompareDistribution):
    '''
    Class compares predicted species-level energy distribution(s) with the
    observed species-level energy distribution(s)

    Attributes
    ----------
    self.observed_data : list of iterables
        Observed species energy distributions (SED)
    self.criteria : a list of dictionaries or None
        If not None, each dictionary specifies the divisions made on the plot
        that generated each SED, IED, and SAD and IED in self.sed_list,
        self.ied_list, and self.sad_list.  All self.criteria should have the
        same length.
    self.dist_list : a list of distribution objects
        Each object is a distribution to which the IEDs in self.ied_list will
        be compared.
    self.sad_spp_list : list of strings or None
        If not None, each string in self.spp_names is a species ID which
        corresponds to an array in self.sed_list.

    '''

    def __init__(self, data_list, dist_list, patch=False):
        '''
        Parameters
        ----------
        data_list : list of tuples or output from Patch object
            A list of tuple where each tuple has length 3.  The first object in
            a tuple is an iterable containing the empirical species energy
            distribution.  The second object is a tuple is a community
            individual energy distribution.  The third object in a tuple is an
            empirical species abundance distribution.        
        dist_list : list of strings or objects
            Each string corresponds to a name of a psi distribution to which to
            compare to the observed data. 
        patch : bool
            If True, expects a tuple of length 3 with the first object being
            the complete output from Patch.sed, the second object being the
            output from Patch.ied and the third element being the output from
            Patch.sad. If False expects what argument data_list describes.
            Empirical sads and energy distributions should be made with the
            same criteria (See Patch class for criteria explanation).

        Notes
        -----
        If data_list is a list of tuples containing iterables, the 1st entry
        (0th element) in each tuple is considered the observed SEDs. 
        '''

        if patch:
            # TODO: Check length of input objects!

            if not ((len(data_list[0]) == len(data_list[1])) and\
                    (len(data_list[1]) == len(data_list[2]))):
                raise IndexError('SED, IED, and SAD patch returns' +\
                    ' must have the same length. Use the same criteria for' +\
                    ' each.')

            #Sort species energy
            sed_criteria = []
            sed_list = []
            spp_names = []
            for obj in data_list[0]:
                spp = list(obj[1].viewkeys()); spp.sort()
                spp_names.append(spp)
                for kw in spp_names[-1]:
                    sed_list.append(obj[1][kw])
                    sed_criteria.append(obj[0])

            #Sort community energy
            ied_criteria = []
            ied_list = []
            for i, obj in enumerate(data_list[1]):

                # For consistency I am copying the ied data for each species 
                num = len(spp_names[i])
                tcri = [obj[0] for i in xrange(num)]
                ied_criteria += tcri
                teng = [obj[1] for i in xrange(num)]
                ied_list += teng

            #Sort sad
            sad_criteria = []
            sad_list = []
            for i, obj in enumerate(data_list[2]):

                # Copy sad data for each species
                num = len(spp_names[i])
                tcri = [obj[0] for i in xrange(num)]
                sad_criteria += tcri
                tsad = [obj[1] for i in xrange(num)]
                sad_list += tsad
            
            self.sad_spp_list = []
            for i in xrange(len(spp_names)):
                self.sad_spp_list += spp_names[i]
            self.criteria = sad_criteria

            super(CompareSED, self).__init__(zip(sed_list, ied_list, sad_list),
                                                                  dist_list, 0)

        else: 
            
            super(CompareSED, self).__init__(data_list, dist_list, 0)
        
    def compare_rads(self, return_spp=False):
        '''
        Comparison of species level energy distributions rank abundance
        distributions.

        Parameters
        ----------
        return_spp : bool
            If True, the returns a tuple with a species list as the second
            element.

        Returns
        -------
        : dict
            Has len(self.dist_list) + 1.  All the distribution class names
            passed to the constructor are key words as well as 'observed' which
            references the observed data. Each keyword looks up
            a list of arrays.  Each list is len(self.ied_list) long and
            contains the predicted reds for the empirical data sets for the
            given distribution. 
        : list or None
            Returns self.sad_spp_list which could be a list of lists or None.
            These names are the species names that correspond numerically with
            the arrays in within each distribution. Only returned if
            return_spp=True.

        '''
        if return_spp:
            return super(CompareSED, self).compare_rads(), self.sad_spp_list
        else:
            return super(CompareSED, self).compare_rads()


    def compare_cdfs(self, return_spp=False):
        '''
        Comparison of species level energy distributions cdfs

        Parameters
        ----------
        return_spp : bool
            If True, the returns a tuple with a species list as the second
            element.

        Returns
        -------
        : dict
            Has len(self.dist_list) + 1.  All the distribution class names
            passed to the constructor are key words as well as 'observed' which
            references the observed data. Each keyword looks up
            a list of arrays.  Each list is len(self.ied_list) long and
            contains the predicted reds for the empirical data sets for the
            given distribution. 
        : list or None
            Returns self.sad_spp_list which could be a list of lists or None.
            These names are the species names that correspond numerically with
            the arrays within each distribution. Only returned if
            return_spp=True.

        '''
        if return_spp: 
            return super(CompareSED, self).compare_cdfs(), self.sad_spp_list
        else:
            return super(CompareSED, self).compare_cdfs()

class CompareASED(CompareDistribution):
    '''
    Compares theoretical and observed ased's

    Attributes
    ----------
    self.observed_data : list of arrays
        Observed average species energy distributions (ASED)
    self.sad_spp_list : list of arrays
        If not None, each array contains the species strings for the
        corresponding SAD in self.sad_list.  The length of self.sad_spp_list
        should be the same length as self.sad_list and the length of any array
        within self.sad_spp_list should be the same length the corresponding array
        in self.sad_list. The index of any species name within any array
        within self.sad_spp_list references the species count with the same
        index in self.sad_list.
    self.criteria : a list of dictionaries or None
        If not None, each dictionary specifies the divisions made on the plot
        that generated each SAD and IED in self.sad_list and self.ied_list.
        self.criteria should be the same length as self.sad_list and
        self.ied_list. 
    self.dist_list : a list of distribution objects
        Each object is a distribution to which the IEDs in self.ied_list will
        be compared.

    '''

    def __init__(self, data_list, dist_list, patch=False):
        '''
        Parameters
        ----------
        data_list : list of tuples or output from Patch object
            A list containing tuples of length three. The first object in the
            tuple is an iterable containing the average energy distribution.
            The second object in a tuple an iterable containing the community
            individual energy distribution.  The third object in a tuple is an
            iterable containing the empirical species abundance
            distribution.See patch argument in this method for information
            about Patch object output.
        dist_list : list of strings or objects
            Each string corresponds to a name of a ased distribution to which to
            compare to the observed data. 
        patch : bool
            If True, expects a tuple of length 3 with the first object being
            the complete output from Patch.ased, the second object being
            the output from Patch.ied and the third element being the
            output from Patch.sad. If False expects what argument data_list
            describes. Empirical sads and energy distributions should be made 
            with the same criteria.

        Notes
        -----
        If data_list is a list of tuples containing iterables, the 1st entry
        (0th element) in each tuple is considered the observed ASEDs. 
        '''

        if patch:

            # Unpack sad. Store spp_lists in items
            sad_criteria, sad_list, sad_spp_lists = \
                                                unpack(data_list[2])

            # Unpack ased 
            ased_criteria, ased_list, ased_species = \
                                                        unpack(data_list[0])

            # Unpack ied
            ied_criteria, ied_list, ied_spp = unpack(data_list[1])

            self.criteria = sad_criteria
            self.sad_spp_list = ased_species

            super(CompareASED, self).__init__(zip(ased_list, ied_list,
                                                       sad_list), dist_list, 0)


        else:
            super(CompareASED, self).__init__(data_list, dist_list, 0)

class CompareSAR(object):
    '''
    Object allows comparison between species-area relationships

    Attributes
    ----------
    self.sar_list : list of arrays
        A list of arrays in which each array is the number of species a
        given areas.  The areas are specified in self.a_list and correspond
        exactly self.sar_list.
    self.a_list : list of arrays
        A list of arrays in which each array is the area (or area fraction) at
        which the number of species specified in self.sar_list are found.
        Indices correspond exactly with self.sar_list.
    self.full_sad : list of arrays
        A list of species abundance distributions (SAD) computed at the anchor
        scale for each given SAR.  The length of self.full_sad should equal the
        length of self.sar_list and self.a_list.
    self.curve_list : list of objects
        A list of SAR curve objects to which the empirical SARs in
        self.sar_list will be compared. 

    '''
    
    def __init__(self, sar_list, curve_list, full_sad, max_a=True, 
                                                                  patch=False):
        '''
        Parameters
        ----------
        sar_list : list of tuples or list of outputs from Patch().sar
            A list of tuples where each tuple contains two array-like objects
            of the same length.  The first element in the tuple is the
            area list and the second element is the species count for the sar.
            The maximum area in the area list should be the anchor area from
            which the full_sad was generated.  If patch=True, accepts the
            output from Patch.sar
        curve_list : list
            A list of SARCurve objects or list of SARCurve object names (str)
        full_sad : list of array-like objects
            List of complete sads.  Each sad corresponds to an element in
            sar_list. 
        max_a : bool
            If max_a is True, compare sets all areas to fractions in area_list.
        patch : bool
            If True, sar_list should be a list of outputs from Patch().sar
        '''

        assert len(sar_list) == len(full_sad), "sar_list and full_sad must " \
                                              + " be the same length"
        self.sar_list = []
        self.a_list = []
        if patch:
             for sar_obj in sar_list:
                 unzipped_sar = unpack(sar_obj[0])
                 self.sar_list.append(np.array(unzipped_sar[0]))
                 self.a_list.append(np.array(unzipped_sar[1]))
        else:
            unzipped_sar = unpack(sar_list)
            self.a_list = [np.array(areas) for areas in unzipped_sar[0]]
            self.sar_list = [np.array(sar) for sar in unzipped_sar[1]]

        # Set to area fractions if max_a is true
        if max_a:
            self.a_list = [ars / np.max(ars) for ars in self.a_list]

        self.full_sad = [np.array(sad) for sad in full_sad]

        self.curve_list = make_dist_list(curve_list)


    def compare_curves(self, iter_vals=False, use_rad=False, form='sar'):
        '''
        Method generates predicted SAR curves from the given observed data and
        curve objects for comparison

        Parameters
        ----------
        use_rad : bool
            If False, uses the sad pmf to calculate the SAR.  If True, uses the
            sad rank abundance distribution to calculate the SAR.
        iter_val : bool
            If True, uses the iterative method to calculate SAR. If False uses
            the one shot method.
        form : string
            Default value is 'sar' which calculates the SAR given the
            parameters. You can also use 'ear' which calculates the EAR with
            the given parameters.

        Returns
        -------
        : list of dicts
            The list is the same length self.sar_list and each dictionary is
            the length of self.curve_list + 1.  Each keyword in a dictionary
            references either the observed SAR ('observed') or the SAR generate by
            one of the curve objects.

        Notes
        -----
        If possible, the SARs are computed using an iterative method.
        Otherwise, they are calculated with a one-shot method.
        '''
        pred_sar = []
        for sar, a, sad in zip(self.sar_list, self.a_list, self.full_sad):
            psar = {}
            psar['observed'] = np.array(zip(sar, a), dtype=[('items', np.float),
                                        ('area', np.float)])
            for cur in self.curve_list:
                cur.fit(sad, (a, sar))

                if iter_vals:
                    try:
                        psar[cur.get_name()] = cur.iter_vals(a,
                                                    use_rad=use_rad, form=form)
                    except AttributeError:
                        psar[cur.get_name()] = cur.iter_vals(a, use_rad=True,
                                                                     form=form)
                else:
                    try:
                        psar[cur.get_name()] = cur.vals(a, use_rad=use_rad,
                                                                     form=form)
                    except AttributeError:
                        psar[cur.get_name()] = cur.vals(a, use_rad=True,
                                                                     form=form)
                    
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
    Generates an empirical cdf from empirical data

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

def aic(neg_L, k, loglik=True):
    '''
    Calculates the AIC of a given model

    Parameters
    ----------
    neg_L : array-like object
        The negative log likelihood of the models or a list of pdfs/pmfs,
        depending on nll
    k : array-like object
        The number of parameters of the model
    loglik : bool
        If True, assumes neg_L is an array-like object of negative log
        likelihood.  If False, assumes neg_L is a list of pdfs/pmfs.
    
   Returns
   -------
   : float
        AIC for a given model
    '''
    if loglik:
        neg_L, k = cnvrt_to_arrays(neg_L, k)
    else:
        neg_L = nll(neg_L)
        neg_L, k = cnvrt_to_arrays(neg_L, k)
        
    assert len(k) == len(neg_L), "neg_L and k must have the same length"
    aic = (2 * neg_L) + (2 * k)
    return aic

def aicc(neg_L, k, n=None, loglik=True):
    '''
    Calculates the corrected AIC of a given model

    Parameters
    ----------
    neg_L : array-like object
        The negative log likelihood of models or list of pdfs/pmfs
    k : array-like object
        The number of parameters of models
    n : array-like object
        Number of observations for each model. Can be left as None if neg_L is
        list of pdfs/pmfs and loglik = True
    loglik : bool
        If True, assumes neg_L is a array-like object of negative log
        likelihood.  If False, assumes neg_L is a list of pdfs/pmfs.

    Returns
    -------
    : np.array
        AICc for a given models

    '''
    if loglik:
        assert n != None, 'n argument must be given if loglik is True'
        neg_L, k, n = cnvrt_to_arrays(neg_L, k, n)
    else:
        n = np.array([len(tneg_L) for tneg_L in neg_L])
        neg_L = nll(neg_L)
        neg_L, k = cnvrt_to_arrays(neg_L, k)

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

def likelihood_ratio(nll_null, nll_alt, df_list):
    '''
    This functions compares of two nested models using the likelihood ratio
    test.

    Parameters
    ----------
    nll_null :  array-like object
        The negative log-likelihood of the null model
    nll_alt : array-like object
        The negative log-likelihood of the alternative model
    df_list : array-like object
        the degrees of freedom calculated as (number of free parameters in
        alternative model) - (number of free parameters in null model)
    
    Returns
    -------
    : list of tuples
        (test_statistic, p-value)

    Notes
    -----
    The LRT only applies to nested models. The variable test_stat is known as
    the G^2 statistic.
    '''
    
    nll_null, nll_alt, df_list = cnvrt_to_arrays(nll_null, nll_alt, df_list)
    assert len(nll_null) == len(nll_alt) and len(nll_null) == len(df_list) and\
           len(nll_alt) == len(df_list), "nll_null, nll_alt, and df_list " + \
                                          "must have the same length"
    # Calculate G^2 statistic
    ll_null = nll_null * -1; ll_alt = nll_alt * -1
    test_stat = 2 * (ll_null - ll_alt) 
    return [(ts, stats.chisqprob(ts, df)) for ts, df in zip(test_stat, df_list)]

def variance(data_sets):
    '''Calculates the variance of the given data_sets
    
    Parameters
    ----------
    data_sets : list
        A list of np.arrays on which the kurtosis will be calculated

    '''

    variance_list = []
    for data in data_sets:
        variance_list.append(np.var(data, ddof=1))

    return variance_list

def skew(data_sets):
    '''Calculates the skew of some given data

    Parameters
    ----------
    data_sets : list
        A list of np.arrays on which the kurtosis will be calculated

    Returns
    -------
    : list
        A list of kurtosis values with the same length as data_sets

    '''

    skewness_list = []
    for data in data_sets:
        skewness_list.append(stats.skew(data))

    return skewness_list 

def kurtosis(data_sets):
    '''Calculates the kurtosis using an online algorithm for the given list of
    datasets

    Parameters
    ----------
    data_sets : list
        A list of np.arrays on which the kurtosis will be calculated

    Returns
    -------
    : list
        A list of kurtosis values with the same length as data_sets

    '''
    kurtosis_list = []
    for data in data_sets:
        kurtosis_list.append(stats.kurtosis(data))

    return kurtosis_list

def bootstrap(data_sets, num_samp=1000):
    '''Bootstrap a data_set within data_sets num_samp times. With replacement

    Parameters
    ----------
    data_sets : list
        A list of np.arrays on which the kurtosis will be calculated
    num_samp : int
        Number of bootstrap samples to take

    Returns
    -------
    : a list
        A list of lists of arrays.  Each list contains num_samp bootstrapped
        arrays
    '''
    
    random.seed(time.time())

    bootstraps = []
    for data in data_sets:
        bt_data = []
        n = len(data)
        for j in xrange(num_samp):
            bt_data.append(np.array([random.choice(data) for j in xrange(n)]))
        bootstraps.append(bt_data)
    
    return bootstraps

def bootstrap_moment(data1, data2, moment, CI=.95, num_samp=1000):
    '''
    A bootstrap two-sample test of a moment. Returns the test_statistic 
    distribution and the confidence interval as specified by parameter CI.  The
    confidence interval is the difference of the moment from data1 minus the
    moment from data2.

    Parameters
    ----------
    data1 : array-like object
        An array like object containing data
    data2 : array-like object
        An array-like object containing data
    moment : list 
        List of strings (mean, skew, kurtosis, and/or variance).
        Will calculate the bootstrap CI's for all the moments in the list
    CI : float
        The desired confidence interval
    num_samp : int
        Number of bootstrap samples

    Returns
    -------
    res : dict
        A dictionary with key words equivalent to the strings found in moment.
        Each keyword looks up tuple with two elements.  The first element is
        the observed difference between the moment of data1 and the moment of
        data2.  The second element is a tuple containing the confidence
        interval (lower_bound, upper_bound) on the difference between the
        specified moment of data1 and data2.
    
    Notes
    -----
    From the returned confidence interval, one is CI confident that the
    returned confidence interval contains the true difference between the
    moment of data1 and data2.  Therefore, if the confidence interval does not
    contain 0 you can be CI confident that the moments are different.

    Bootstrapping in typically only appropriate for sample sizes >= 25. 


    '''

    data1 = np.array(data1)
    data2 = np.array(data2)
    # Bootstrap the data
    data1_boot = bootstrap([data1], num_samp=num_samp)[0]
    data2_boot = bootstrap([data2], num_samp=num_samp)[0]

    def calc_ci(stat1, stat2):
        """ Calculate CI """
        
        diff = stat1 - stat2
        lci = (1 - CI) / 2.
        uci = 1 - lci
        ci = (stats.scoreatpercentile(diff, 100 * lci),\
          stats.scoreatpercentile(diff, 100 * uci))
        return ci

    
    res = {}
    # Set the higher order moment
    if 'skew' in moment:

        stat_1 = np.array(skew(data1_boot))
        stat_2 = np.array(skew(data2_boot))
        
        stat_dist = skew([data1])[0] - skew([data2])[0]
        ci = calc_ci(stat_1, stat_2)

        res['skew'] = (stat_dist, ci)

    if 'variance' in moment:
        stat_1 = np.array(variance(data1_boot))
        stat_2 = np.array(variance(data2_boot))
        
        stat_dist = variance([data1])[0] - variance([data2])[0]
        ci = calc_ci(stat_1, stat_2)

        res['variance'] = (stat_dist, ci)

    if 'kurtosis' in moment:
        stat_1 = np.array(kurtosis(data1_boot))
        stat_2 = np.array(kurtosis(data2_boot))
        
        stat_dist = kurtosis([data1])[0] - kurtosis([data2])[0]
        ci = calc_ci(stat_1, stat_2)

        res['kurtosis'] = (stat_dist, ci)

    if "mean" in moment:
        stat_1 = np.array([np.mean(bs) for bs in data1_boot])
        stat_2 = np.array([np.mean(bs) for bs in data2_boot])
        
        stat_dist = np.mean(data1) - np.mean(data2) 
        ci = calc_ci(stat_1, stat_2)

        res['mean'] = (stat_dist, ci)

    return res 

def mean_squared_error(obs, pred):
    '''
    Calculates the mean squared error between observed and predicted data sets.
    The data sets must be of the same length
    
    Parameters
    ----------
    obs : array-like object
        The observed data
    pred : array-like object
        The predicted data

    Returns
    -------
    : float
        The mean squared error
    '''

    if len(obs) != len(pred):
        raise ValueError('obs and pred parameters must have the same length')

    obs, pred = cnvrt_to_arrays(obs, pred)

    return sum((pred - obs)**2) / len(obs)


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

def make_dist_list(dist_list):
    '''
    If the dist_list is all strings, eval them.  Else return as is
    '''

    if np.all([type(dist) == str for dist in dist_list]):

        ret_dist_list = np.empty(len(dist_list), dtype=object)

        for i, dist_obj in enumerate(dist_list):

            # Clean strings
            dist_obj = dist_obj.strip()
            try:
                ret_dist_list[i] = eval(dist_obj + '()')
            except:
                # Do this if passing in a gen_sar sad and ssad
                # Assumes the sad and ssad are separated by '-'
                try:
                    sad, ssad = tuple(dist_obj.split('-'))
                    if sad.find('(') != 1 and sad.find(')') != -1:
                        sad_obj = eval(sad.strip())
                    else:
                        sad_obj = eval(sad.strip() + '()')
                    if ssad.find('(') != 1 and ssad.find(')') != -1:
                        ssad_obj = eval(ssad.strip())
                    else:
                        ssad_obj = eval(ssad.strip() + '()')
                    ret_dist_list[i] = gen_sar(sad_obj, ssad_obj)
                except:
                    raise NameError("Could not evaluate '%s' as an object name"
                    % dist_obj + '. It may not exist or may be improperly' + 
                    ' formatted. Please check your distribution list in ' 
                    + 'your parameters.xml file or in the dist_list' + 
                    " argument '%s'" % str(dist_list))

        ret_dist_list = list(ret_dist_list)
    else:
        ret_dist_list = dist_list

    return ret_dist_list

def unpack(zipped_data):
    '''
    Unpacks zipped data

    '''

    unzipped_data = zip(*zipped_data)
    unzipped_data = [list(tup) for tup in unzipped_data]
    return tuple(unzipped_data)
 
