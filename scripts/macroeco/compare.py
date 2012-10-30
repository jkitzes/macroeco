#!/usr/bin/python

''' This module contains classes and functions for comparing empirical and
predicted macroecological metrics. 

Classes
-------
CompareDistribution : Base class to for CompareSAD and CompareSSAD

CompareSAD : Compares predicted species abundance distributions with empirical
SADs

CompareSSAD : Compares predicted species-level spatial abundance distributions
with empirical SSADS

CompareSAR : Compares predicted species-area relationship (SAR) curves 
with empirical SAR curves

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


class CompareDistribution(object):
    '''
    Comparison object compares a list of data to any number of distributions

    '''
    
    #TODO: Error Checking
    def __init__(self, data_list, dist_list, patch=None):
        '''
        Parameters
        ----------
        data_list : list of iterables or output from Patch object
            List of np.arrays containing data
        dist_list : list
            List of distribution objects or strings that have the same name as 
            a distribution object. If they are strings, they will be evaled 
        patch : str
            If 'sad', expects the output from the Patch.sad method and if
            'ssad' expects the output from the Patch.ssad method. If None, 
            expects a list of iterables.

        '''

        if patch == 'sad':
            self.criteria, self.data_list, self.spp_list = unpack(data_list)
        elif patch == 'ssad':
            self.spp_list = list(data_list[1].viewkeys())
            self.data_list = [np.array(data_list[1][nm]) for nm in
                                                            self.spp_list]
            self.criteria = data_list[0]
        else:
            self.data_list = [np.array(data) for data in data_list]
            self.criteria = None
            self.spp_list = None
       
        self.dist_list = make_dist_list(dist_list)

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
            #NOTE: dist.par_num is the number of parameters of distribution
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
            If True, calculates the corrected AIC weights for the given data. 
            If False, calculates AIC weights.

        Returns
        -------
        : tuple
            first element is a list of arrays with each array having length 
            equal to the number of models proposed and the length of the list 
            is the lenth of self.data_lists. The second element is the delta
            AIC values in the same format as the first tuple object. The third
            object are the AIC values in the same format as the output of the 
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
            passed to the constructor are key words as well as 'observed' which
            references the observed data, self.data_list. Each keyword looks up
            a list of arrays.  Each list is len(self.data_list) long and
            contains the predicted rads for the empirical data sets for the
            given distribution. 

        '''
        rads_dict = {}
        rads_dict['observed'] = self.data_list
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
            passed to the constructor are key words as well 'observed' which
            references the observed data, self.data_list. Each keyword looks up
            a list of arrays.  Each list is len(self.data_list) long and
            contains the predicted cdfs for the empirical data sets for the
            given distribution. 


        '''

        cdfs_dict = {}
        cdfs_dict['observed'] = [empirical_cdf(data) for data in self.data_list]
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

    def compare_rarity(self, rads, mins_list):
        '''
        This script takes in the output from self.compare_rads and a list of
        minimum values against which to compare the observed and predicted
        rads.
        and outputs a dictionary with length self.dist_list + 1 (all
        distributions + observed).  Each keyword in this dict looks up a dict
        of len(mins_list) where the keywords are the values against which the
        rads will be <=.  Each one of these subdictionaries looks up a list
        with len(self.data_list).

        Parameters
        ----------
        rads : dict
            The output from self.compare_rads(). 
        mins_list : array-like object
            A list of numbers.  Each number number will be used in the
            following function: rad <= mins_list[i].

        Returns
        -------
        : dict
            Returns a dictionary with length self.dist_list + 1 (all
        distributions + observed).  Each keyword in this dict looks up a dict
        of len(mins_list) where the keywords are the values against which the
        rads will be <=.  Each one of these subdictionaries looks up a list
        with len(self.data_list).


        '''
        mins_list = make_array(mins_list)

        rarity = {}
        keys = list(rads.viewkeys())
        for kw in keys:
            rarity[kw] = {}
            for mins in mins_list:
                rarity[kw][mins] = [sum(data <= mins) for data in rads[kw]]
        return rarity


    def summary(self, rads, mins_list=[10], crt=False):
        '''
        Summarizes the given datasets and the predicted rads. Looks at
        total balls sampled ('balls'), number of urns ('urns'), the max balls
        in a given urn ('max'), number of urns with less than MIN balls ('tot
        <= MIN'), and the fit of the the distributions in self.dist_list to the
        data in self.data_list

        Parameters
        ----------
        rads : dict
            The dictionary returned by self.compare_rads()
        mins_list : list
            Bins with balls less than or equal to 10
        crt : bool
            If True, corrected AIC, if False, not

        Returns
        -------
        : dict
            Dictionary of dictionaries of length self.dist_list + 1.  Each
            sub-dictionary other than 'observed' contains the keywords balls, urns,
            max, tot_min, aic, aic_d, aic_w, and par_num.  Each of these
            keywords contains a list that is the same length as the number of
            sads under consideration.

            aic = AIC
            aic_d = Delta AIC
            aic_w = AIC weights
            par_num = Parameer number of the given distribution


        '''
        summary = {}
        rarity = self.compare_rarity(rads, mins_list=mins_list)
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

        return summary

class CompareSAD(CompareDistribution):
    '''
    Object inherits CompareDistribution and uses it to compare species
    abundance distributions (SAD)

    Attributes
    ----------
    self.data_list : A list of arrays
        Each array in this list is an SAD.  Each of these SADs will be compared
        to the distributions in self.dist_list
    self.dist_list : a list of distribution objects
        Each object is a distribution object to which the SADs in
        self.data_list will be compared.  
    self.criteria : a list of dictionaries or None
        If not None, each dictionary specifies the divisions made on the plot
        that generated each SAD in self.data_list.  self.criteria should be the
        same length as self.data_list
    self.spp_list : list of arrays or None
        If not None, each array contains the species strings for the
        corresponding SAD in self.data_list.  The length of self.spp_list
        should be the same length as self.data_list and the length of any array
        within self.spp_list should be the same length the corresponding array
        in self.data_list. The indice of any species name within any array
        within self.spp_list references the species count with the same
        indice in self.data_list.

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
            expects a list of iterables. Persumably, each iterable is an SAD.
        '''
        if patch == True:
            super(CompareSAD, self).__init__(data_list, dist_list, patch='sad')
        else:
            super(CompareSAD, self).__init__(data_list, dist_list, patch=None)

class CompareSSAD(CompareDistribution):
    '''
    Object inherits CompareDistribution and uses it to compare species-level
    spatial abundance distributions (SSAD)

    Attributes
    ----------
    self.data_list : A list of arrays
        Each array in this list is an SSAD.  Each of these SSADs will be
        compared to the distributions in dist_list
    self.dist_list : a list of distribution objects
        Each object is a distribution object to which the SSADs in
        self.data_list will be compared.  
    self.criteria : a list of dictionaries or None
        If not None, each dictionary specifies the divisions made on the plot
        that generated each SAD in self.data_list.  self.criteria should be the
        same length as self.data_list
    self.spp_list : List of strings or None
        If not None, self.spp_list is a list of strings where each string
        refers to a species.  The length of self.spp_list should be the same
        length as self.data_list.  Each species string has the same index
        within the list as its corresponding SSAD in self.data_list. 

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
            expects a list of iterables. Persumably, each iterable is an SAD.
        '''
        if patch == True:
            super(CompareSSAD, self).__init__(data_list, dist_list, \
                                                                  patch='ssad')
        else:
            super(CompareSSAD, self).__init__(data_list, dist_list, patch=None)

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


    def compare_curves(self):
        '''
        Method generates predicted SAR curves from the given observed data and
        curve objects for comparison

        Returns
        -------
        : list of dicts
            The list is the same length self.sar_list and each dictionary is
            the length of self.curve_list + 1.  Each keyword in a dictionary
            references either the observed SAR ('observed') or the SAR generate by
            one of the curve objects. 
        '''
        pred_sar = []
        for sar, a, sad in zip(self.sar_list, self.a_list, self.full_sad):
            psar = {}
            psar['observed'] = np.array(zip(sar, a), dtype=[('items', np.float),
                                        ('area', np.float)])
            for cur in self.curve_list:
                cur.fit(sad, (a, sar))
                psar[get_name(cur)] = cur.vals(a)
            for kw in psar.iterkeys():
                psar[kw].sort(order='area')
            pred_sar.append(psar)
        return pred_sar

class CompareIED(object):
    '''
    Class compares predicted individual energy distributions (IED) for the entire
    community to observed IEDs

    Attributes
    ----------
    self.ied_list : list of arrays
        Each array in this list is an IED.  Each of these IEDs will be compared
        to the distributions in self.dist_list
    self.ied_items : list of arrays
        Each array contains species strings which pair to the values
        contained in the corresponding array in self.ied_list. The length of
        self.ied_items should be the same length as self.ied_list.
    self.sad_list : list of arrays
        Each array in this list is an species abundance distribution (SAD).  
        self.sad_list should have the same length as self.ied_list.  Each SAD
        corresponds the IED with the same index and should have been generated
        with the same criteria. 
    self.sad_items : list of arrays
        If not None, each array contains the species strings for the
        corresponding SAD in self.sad_list.  The length of self.sad_items
        should be the same length as self.sad_list and the length of any array
        within self.sad_items should be the same length the corresponding array
        in self.sad_list. The indice of any species name within any array
        within self.sad_items references the species count with the same
        indice in self.sad_list.
    self.*_criteria : a list of dictionaries or None
        If not None, each dictionary specifies the divisions made on the plot
        that generated each SAD and IED in self.sad_list and self.ied_list.
        self.*_criteria should be the same length as self.sad_list and
        self.ied_list. 
    self.dist_list : a list of distribution objects
        Each object is a distribution to which the IEDs in self.ied_list will
        be compared.

    '''

    def __init__(self, data_list, dist_list, patch=False):
        '''
        Parameters
        ----------
        data_list : tuple of lists or output from Patch object
            Tuple of length two. The first object in the tuple is a list of
            np.arrays that are the empirical community energy distributions.
            The second object is a list of np.arrays that are the empirical
            sads. See patch argument for more information.
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
        '''

        if patch:
            # Unpack sad. Store spp_lists in items
            self.sad_criteria, self.sad_list, self.sad_items = \
                                                    unpack(data_list[1])

            # Unpack ied
            self.ied_criteria, self.ied_list, self.ied_items = \
                                                        unpack(data_list[0])
            
        else:
            self.sad_list = [np.array(data) for data in data_list[1]]
            self.sad_criteria = None
            self.sad_items = None
            self.ied_list = [np.array(data) for data in data_list[0]]
            self.ied_criteria = None

        # Clean SADs
        self.sad_list = [data[data != 0] for data in self.sad_list]

        # Check if dist_list is a list of strings or list of objects
        self.dist_list = make_dist_list(dist_list)

        assert len(self.sad_list) == len(self.ied_list), 'SAD list and' \
                                        + 'IED list must be the same length'

    def compare_reds(self):
        '''
        Compare the rank energy distributions against the given distribution
        list

        Returns
        -------
        : dict
            Has len(self.dist_list) + 1.  All the distribution class names
            passed to the constructor are key words as well as 'observed' which
            references the observed data. Each keyword looks up
            a list of arrays.  Each list is len(self.ied_list) long and
            contains the predicted reds for the empirical data sets for the
            given distribution. 

        '''

        reds_dict = {}
        reds_dict['observed'] = self.ied_list
        for dist in self.dist_list:
            dist.fit((self.ied_list, self.sad_list))
            reds_dict[get_name(dist)] = dist.rad()
        return reds_dict

class CompareSED(object):
    '''
    Class compares predicted species-level energy distribution (s) with the
    observed species-level energy distribution (s)

    Attributes
    ----------
    self.sed_list : list of arrays
        Each array in this list is a species-level energy distribution (SED).
        Each SED will be compared to the distributions given in dist_list.
    self.ied_list : list of arrays
        Each array in this list is an individual energy distribution (IED).
        self.ied_list will always have the same length as self.sed_list.
        CompareSED will broadcast ied_list to sed_list if only one IED is
        given. 
    self.sad_list : list of arrays
        Each array in this list is an species abundance distribution (SAD).  
        self.sad_list should have the same length as self.ied_list and
        self.sed_list.  If only one SAD is given, it will be broadcast to the
        length of self.sed_list.   
    self.*_criteria : a list of dictionaries or None
        If not None, each dictionary specifies the divisions made on the plot
        that generated each SED, IED, and SAD and IED in self.sed_list,
        self.ied_list, and self.sad_list.  All self.*_criteria should have the
        same length.
    self.dist_list : a list of distribution objects
        Each object is a distribution to which the IEDs in self.ied_list will
        be compared.
    self.spp_names : list of strings or None
        If not None, each string in self.spp_names is a species ID which
        corresponds to an array in self.sed_list.

    '''

    def __init__(self, data_list, dist_list, patch=False):
        '''
        Parameters
        ----------
        data_list : tuple of lists or output from Patch object
            Tuple of length three. The first object is a list np.arrays that
            contains species energy distributions. The second object in the
            tuple is a list of np.arrays that are the empirical community
            energy distributions. The third object is a list of np.arrays that
            are the empirical SADs. All three lists should should be the same
            length. See patch argument in this method for information about
            Patch object output.
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
        '''

        if patch:
            # TODO: Check length of input objects!

            #Sort species energy
            self.sed_criteria = []
            self.sed_list = []
            self.spp_names = []
            for obj in data_list[0]:
                self.spp_names.append(list(obj[1].viewkeys()))
                for kw in self.spp_names[-1]:
                    self.sed_list.append(obj[1][kw])
                    self.sed_criteria.append(obj[0])

            #Sort community energy
            self.ied_criteria = []
            self.ied_list = []
            for i, obj in enumerate(data_list[1]):

                # For consistency I am copying the ied data for each species 
                num = len(self.spp_names[i])
                tcri = [obj[0] for i in xrange(num)]
                self.ied_criteria += tcri
                teng = [obj[1] for i in xrange(num)]
                self.ied_list += teng

            #Sort sad
            self.sad_criteria = []
            self.sad_list = []
            for i, obj in enumerate(data_list[2]):

                # Copy sad data for each species
                num = len(self.spp_names[i])
                tcri = [obj[0] for i in xrange(num)]
                self.sad_criteria += tcri
                tsad = [obj[1] for i in xrange(num)]
                self.sad_list += tsad

            self.spp_names = list(np.array(self.spp_names).flatten())

        else:
            lr = [len(dt) for dt in data_list] # len list for checking
            if lr[1] != lr[0]:
                if lr[1] == 1:
                    self.ied_list = [np.array(data_list[1][0]) for i in\
                                                                 xrange(lr[0])]
                else:
                    raise AttributeError("Length of IED list must be 1 or the"\
                                         + " same length as the SED list")
            else:
                self.ied_list = [np.array(data) for data in data_list[1]]

            if lr[2] != lr[0]:
                if lr[2] == 1:
                    self.sad_list = [np.array(data_list[2][0]) for i in\
                                                                 xrange(lr[0])]
                else:
                    raise AttributeError("Length of SAD list must be 1 or the"\
                                         + " same length as the SED list")
            else:
                self.sad_list = [np.array(data) for data in data_list[2]]

            self.sad_criteria = None
            self.ied_criteria = None
            self.sed_list = [np.array(data) for data in data_list[0]]
            self.sed_criteria = None
            self.spp_names = None

        # Clean SADs just in case they have not been cleaned
        self.sad_list = [data[data != 0] for data in self.sad_list]        

        self.dist_list = make_dist_list(dist_list)
        
        # Double checking lengths
        assert len(self.sad_list) == len(self.ied_list), 'sad_list and' \
                                        + ' eng_list must be the same length'
        assert len(self.sed_list) == len(self.ied_list), 'speng_list and ' \
                                        + ' eng_list must be the same length'
    
    def compare_reds(self):
        '''
        Comparison of species level energy distributions

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
            Returns self.spp_names which could be a list of lists or None.
            These names are the species names that correspond numerically with
            the arrays in within each distribution.

        '''
        reds_dict = {}
        reds_dict['observed'] = self.sed_list
        for dist in self.dist_list:
            dist.fit((self.sed_list, self.ied_list, self.sad_list))
            reds_dict[get_name(dist)] = dist.rad()

        return reds_dict, self.spp_names

class CompareASED(object):
    '''
    Compares theoretical and observed ased's

    '''

    def __init__(self, data_list, dist_list, patch=False):
        '''
        Parameters
        ----------
        data_list : tuple of lists or output from Patch object
            Tuple of length three. The first object is a list np.arrays that
            are the the average species energy distributions. 
            The second object in the tuple is a list of np.arrays that are the 
            empirical community energy distributions. The third object is a 
            list of np.arrays that are the empirical sads. See patch argument 
            in this method for information about Patch object output.
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
        '''

        if patch:

            # Unpack sad. Store spp_lists in items
            self.sad_criteria, self.sad_list, self.sad_items = \
                                                unpack(data_list[2])

            # Unpack ased 
            self.ased_items, self.ased_criteria, self.ased_list = \
                                                        unpack(data_list[0])

            # Unpack ied
            self.ied_criteria, self.ied_list = unpack(data_list[1])
        else:
            self.sad_list = [np.array(data) for data in data_list[2]]
            self.sad_items = None
            self.sad_criteria = None
            self.ased_list = [np.array(data) for data in data_list[0]]
            self.ased_items = None
            self.ased_criteria = None
            self.ied_list = [np.array(data) for data in data_list[1]]
            self.ied_criteria = None

        # Clean SADs
        self.sad_list = [data[data != 0] for data in self.sad_list]

        self.dist_list = make_dist_list(dist_list)

    def compare_rads(self):
        '''
        Comparison of average species energy distribution

        Returns
        -------
        : dict
            Has len(self.dist_list) + 1.  All the distribution class names
            passed to the constructor are key words as well as 'observed' which
            references the observed data. Each keyword looks up
            a list of arrays.  Each list is len(self.ied_list) long and
            contains the predicted reds for the empirical data sets for the
            given distribution.
        '''
        rads_dict = {}
        rads_dict['observed'] = self.ased_list
        for dist in self.dist_list:
            dist.fit((self.ied_list, self.sad_list))
            #Different Identifier?
            rads_dict[get_name(dist)] = dist.rad()

        return rads_dict

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

def make_dist_list(dist_list):
    '''
    If the dist_list is all strings, eval them.  Else return as is
    '''

    if np.all([type(dist) == str for dist in dist_list]):

        ret_dist_list = np.empty(len(dist_list), dtype=object)

        for i, dist_obj in enumerate(dist_list):

            # Clean strings
            dist_obj = dist_obj.strip()
            ret_dist_list[i] = eval(dist_obj + '()')
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
 
