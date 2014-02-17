#!/usr/bin/python

#Testing Compare Module

import unittest
from macroeco.compare import *
import numpy as np
import scipy.stats as stats
import copy
import macroeco.distributions as dist
import numpy.testing as nt

class TestCompare(unittest.TestCase):
    '''Test classes and methods in compare.py'''

    def setUp(self):
        self.sad_data = [[1,1,1,1,1,2,3,4,5,6], [2,2,2,2,2,2,2,2,2,2]]
        self.ssad_data= [[0,0,0,1,1,2,3,5,12], (0,1,1,1,2,6,12)]


    def test_CompareSAD_init(self):

        # Test that not passing in patch object object works
        sad_c = CompareSAD(self.sad_data, ['logser'])
        
        # Check that sad_data became self.observed_data
        sums = np.array([sum(x) for x in sad_c.observed_data])
        test_sums = np.array([sum(x) for x in self.sad_data])
        self.assertTrue(np.all(sums == test_sums))

        # Test that that other attributes were set correctly
        self.assertTrue(sad_c.criteria == None)
        self.assertTrue(sad_c.sad_spp_list == None)
        
        # Test that distribution object was fit
        self.assertTrue(np.all(sad_c.dist_list[0].params['tot_obs'] ==
                                                                    test_sums))
        self.assertTrue(np.all(sad_c.dist_list[0].params['n_samp'] ==
                                                            np.array([10,10])))

        # Test if patch is true!

        # Replica of patch output
        patch_true = [({'test' : 'criteria'}, np.array([1,1,1,2,3,5]),
                     np.array(['a', 'b', 'c', 'd', 'e', 'g'])), ({'test' : 
                     'criteria'}, np.array([1,1,1,2,5]), np.array(['a', 'b',
                     'c', 'd', 'g']))]
        sad_c = CompareSAD(patch_true, dist_list=['logser'], patch=True)

        # Test that the parsing happened correctly
        self.assertTrue(len(sad_c.criteria) == 2)
        self.assertTrue(len(sad_c.sad_spp_list) == 2)
        self.assertTrue(len(sad_c.observed_data) == 2)

        # Check that parameter values were fit correctly
        self.assertTrue(np.all(sad_c.dist_list[0].params['n_samp'] ==
                                                              np.array([6,5])))
        self.assertTrue(np.all(sad_c.dist_list[0].params['tot_obs'] ==
                                                           np.array([13, 10])))

        # Check that the species lists were set correctly
        self.assertTrue(np.all(sad_c.sad_spp_list[0] == 
                                    np.array(['a', 'b', 'c', 'd', 'e', 'g'])))
        self.assertTrue(np.all(sad_c.sad_spp_list[1] == 
                                    np.array(['a', 'b', 'c', 'd', 'g'])))

    def test_CompareSSAD_init(self):
        
        # Test that SSAD parses correctly when patch is False
        # Test that not passing in patch object object works
        ssad_c = CompareSSAD(self.ssad_data, ['binm'])
        
        # Check that sad_data became self.observed_data
        sums = np.array([sum(x) for x in ssad_c.observed_data])
        test_sums = np.array([sum(x) for x in self.ssad_data])
        self.assertTrue(np.all(sums == test_sums))

        # Test that that other attributes were set correctly
        self.assertTrue(ssad_c.criteria == None)
        self.assertTrue(ssad_c.sad_spp_list == None)
        
        # Test that distribution object was fit
        self.assertTrue(np.all(ssad_c.dist_list[0].params['tot_obs'] ==
                                                                    test_sums))
        self.assertTrue(np.all(ssad_c.dist_list[0].params['n_samp'] ==
                                                            np.array([9,7])))

        # Test that ssad parses correctly if patch=True
        ssad_patch = (np.array([{}, {}, {}, {}, {}]), {'spp1' :
                        np.array([0,0,1,2,4]), 'spp2' : np.array([1,1,1,1,1])})

        ssad_c = CompareSSAD(ssad_patch, dist_list = ['tgeo', 'binm'],
                    patch=True)
        
        spp_list = np.array(['spp1', 'spp2'])
        self.assertTrue(np.all(spp_list == np.sort(ssad_c.sad_spp_list)))

        # Test that distribution object was fit
        self.assertTrue(np.all(ssad_c.dist_list[0].params['tot_obs'] ==
                                                    np.array([7, 5])))
        self.assertTrue(np.all(ssad_c.dist_list[0].params['n_samp'] ==
                                                            np.array([5,5])))
        # Test that distribution object was fit
        self.assertTrue(np.all(ssad_c.dist_list[1].params['tot_obs'] ==
                                                            np.array([7,5])))
        self.assertTrue(np.all(ssad_c.dist_list[1].params['n_samp'] ==
                                                            np.array([5,5])))

        self.assertTrue(len(ssad_c.criteria) == 5)

    def test_CompareIED_init(self):
        
        # Test the CompareIED init parses correctly
        ied_data = [(np.arange(10,100), np.arange(1,40)), (np.arange(1,20),
                                                             np.arange(40,60))]
        ied_c = CompareIED(ied_data, dist_list=['psi'])

        # Check the first item in tuple became observed data
        sums = np.array([sum(x) for x in ied_c.observed_data])
        test_sums = np.array([sum(np.arange(10,100)), sum(np.arange(1,20))])
        self.assertTrue(np.all(sums == test_sums))

        self.assertTrue(ied_c.criteria == None)
        self.assertTrue(ied_c.sad_spp_list == None)

        # Test that distribution object was fit including E parameter
        self.assertTrue(np.all(ied_c.dist_list[0].params['tot_obs'] ==
                     np.array([sum(np.arange(1,40)), sum(np.arange(40,60))])))
        self.assertTrue(np.all(ied_c.dist_list[0].params['n_samp'] ==
                                                            np.array([39,20])))
        self.assertTrue(np.all(ied_c.dist_list[0].params['E'] ==
                    np.array([sum(np.arange(10,100)),sum(np.arange(1,20))])))

        # If patch is True, make sure the fit works
        patch_sad = [({'test' : 'criteria'}, np.array([1,1,1,2,3,5]),
                     np.array(['a', 'b', 'c', 'd', 'e', 'g'])), ({'test' : 
                     'criteria'}, np.array([1,1,1,2,5]), np.array(['a', 'b',
                     'c', 'd', 'g']))]

        patch_ied = [({}, np.arange(1,40), np.repeat('a', 39)), ({},
                                          np.arange(1,30), np.repeat('b', 29))]

        ied_c = CompareIED((patch_ied, patch_sad), dist_list=['nu'], patch=True)

        # Check ied_list and spp_list
        sad_spp = [np.array(['a', 'b', 'c', 'd', 'e', 'g']), 
                                          np.array(['a', 'b', 'c', 'd', 'g'])]
        bools = [np.all(a == b) for a,b in zip(np.array(ied_c.sad_spp_list),
                                                            np.array(sad_spp))]
        self.assertTrue(np.all(bools))

        ied_spp = [np.repeat('a',39), np.repeat('b',29)]
        bools = [np.all(a == b) for a,b in zip(ied_spp, ied_c.ied_spp_lists)]
        self.assertTrue(np.all(bools))

        # check criteria is right length
        self.assertTrue(len(ied_c.criteria) == 2)

        # Check that observed data is correct
        bools = [np.all(a == b) for a,b in zip(ied_c.observed_data,
                                          [np.arange(1,40),  np.arange(1,30)])]
        self.assertTrue(np.all(bools))

        # Check the fit of distribution
        self.assertTrue(np.all(ied_c.dist_list[0].params['tot_obs'] ==
                     np.array([13, 10])))
        self.assertTrue(np.all(ied_c.dist_list[0].params['n_samp'] ==
                                                            np.array([6,5])))
        self.assertTrue(np.all(ied_c.dist_list[0].params['E'] ==
                    np.array([sum(np.arange(1,40)),sum(np.arange(1,30))])))

    def test_CompareSED_init(self):

        # Test that all attributes are set correctly (sed, ied, sad)
        sed_data = [(np.arange(1,20), np.arange(1,40), np.arange(5,25)),
                           (np.arange(1,30), np.arange(5,30), np.arange(4,64))]
        
        sed_c = CompareSED(sed_data, dist_list=['theta'])

        # Did other attributes set correctly?
        self.assertTrue(sed_c.criteria == None)
        self.assertTrue(sed_c.sad_spp_list == None)

        # Check if observed sed data set correctly
        test_obs = [np.arange(1,20), np.arange(1,30)]
        bools = [np.all(a == b) for a,b in zip(sed_c.observed_data, test_obs)] 
        self.assertTrue(np.all(bools))

        # Check that distribution fit correctly
        self.assertTrue(np.all(sed_c.dist_list[0].params['tot_obs'] == 
                np.array([sum(np.arange(5,25)), sum(np.arange(4,64))]))) 
        self.assertTrue(np.all(sed_c.dist_list[0].params['n_samp'] == 
                np.array([len(np.arange(5,25)), len(np.arange(4,64))]))) 
        self.assertTrue(np.all(sed_c.dist_list[0].params['n'] == 
                np.array([len(np.arange(1,20)), len(np.arange(1,30))])))
        self.assertTrue(np.all(sed_c.dist_list[0].params['E'] == 
                np.array([sum(np.arange(1,40)), sum(np.arange(5,30))])))

        # Test if patch == True
        patch_sed = [({}, {'a' : np.arange(1,10), 'b' : np.arange(1,20), 'c':
            np.arange(1,30), 'd' : np.arange(1,40)}), ({}, 
            {'a' : np.arange(1,10), 'b' : np.arange(1,20), 'c': 
            np.arange(1,30), 'd' : np.arange(1,40)})]

        patch_sad = [({}, np.arange(1,50), np.repeat('d',20))]
        patch_ied = [({}, np.arange(4,67), np.repeat('y', 60))]

        # An error should be raised if sed,ied, and sad don't have the same
        # length
        self.assertRaises(IndexError, CompareSED, (patch_sed, patch_ied,
                                   patch_sad), dist_list=['theta'], patch=True)
        
        
        patch_sad = [({}, np.arange(1,50), np.repeat('d',20)), 
                    ({}, np.arange(1,50), np.repeat('d',20))]
        patch_ied = [({}, np.arange(4,67), np.repeat('y', 60)),
                     ({}, np.arange(4,67), np.repeat('y', 60))]

        sed_c = CompareSED((patch_sed, patch_ied, patch_sad),
                                               dist_list=['theta'], patch=True)

        # Check that observed data is set correctly
        self.assertTrue(len(sed_c.observed_data) == 8)
        test_obs = [np.arange(1,10), np.arange(1,20), np.arange(1,30),
                    np.arange(1,40)]
        test_obs += test_obs
        bools = [np.all(a == b) for a,b in zip(test_obs, sed_c.observed_data)]
        self.assertTrue(np.all(bool))

        # Check distributions fit correctly
        nt.assert_array_equal(sed_c.dist_list[0].params['n'], np.array([9,
                                                   19, 29, 39, 9, 19, 29, 39]))
        nt.assert_array_equal(sed_c.dist_list[0].params['E'],
            np.repeat(sum(np.arange(4,67)), 8))
        nt.assert_array_equal(sed_c.dist_list[0].params['tot_obs'],
            np.repeat(sum(np.arange(1,50)), 8))
        nt.assert_array_equal(sed_c.dist_list[0].params['n_samp'],
            np.repeat(len(np.arange(1,50)), 8))
        
        # Check that the species list is correct
        nt.assert_array_equal(np.array(['a', 'b', 'c', 'd', 'a', 'b', 'c',
                                        'd']), np.array(sed_c.sad_spp_list))

        # Check that criteria is correct length
        self.assertTrue(len(sed_c.criteria) == 8)

    def test_CompareASED_init(self):
        
        # Test that ased fits correctly
        
        ased_data = [(np.arange(1,10), np.arange(4,56), np.arange(1,20)),
                           (np.arange(1,34), np.arange(3,20), np.arange(1,56))]

        ased_c = CompareASED(ased_data, dist_list=['nu'])

        # Did other attributes set correctly?
        self.assertTrue(ased_c.criteria == None)
        self.assertTrue(ased_c.sad_spp_list == None)

        # Check if observed ased data set correctly
        test_obs = [np.arange(1,10), np.arange(1,34)]
        bools = [np.all(a == b) for a,b in zip(ased_c.observed_data, test_obs)] 
        self.assertTrue(np.all(bools))

        # Check that distribution fit correctly
        self.assertTrue(np.all(ased_c.dist_list[0].params['tot_obs'] == 
                np.array([sum(np.arange(1,20)), sum(np.arange(1,56))]))) 
        self.assertTrue(np.all(ased_c.dist_list[0].params['n_samp'] == 
                np.array([len(np.arange(1,20)), len(np.arange(1,56))]))) 
        self.assertTrue(np.all(ased_c.dist_list[0].params['E'] == 
                np.array([sum(np.arange(4,56)), sum(np.arange(3,20))])))

        # Test if patch == True
        patch_ased = [({}, np.arange(1,50), np.repeat('d',20)), 
                    ({}, np.arange(1,50), np.repeat('e',20))]
        patch_sad = [({}, np.arange(1,50), np.repeat('d',20)), 
                    ({}, np.arange(1,50), np.repeat('e',20))]
        patch_ied = [({}, np.arange(4,67), np.repeat('y', 60)),
                     ({}, np.arange(4,67), np.repeat('y', 60))]

        ased_c = CompareASED((patch_ased, patch_ied, patch_sad),
                                                  dist_list=['nu'], patch=True)

        # Test that species list is correct
        test_spp = [np.repeat('d', 20), np.repeat('e', 20)]
        nt.assert_array_equal(test_spp, ased_c.sad_spp_list)

        # Test that observed data is correct
        nt.assert_array_equal(ased_c.observed_data, [np.arange(1,50),
                                                            np.arange(1,50)])

        # Test that fit distribution is correct
        nt.assert_array_equal(ased_c.dist_list[0].params['tot_obs'], 
                              np.array([1225, 1225]))
        nt.assert_array_equal(ased_c.dist_list[0].params['n_samp'], 
                              np.array([49, 49]))
        nt.assert_array_equal(ased_c.dist_list[0].params['E'], 
                              np.array([sum(np.arange(4,67)),
                              sum(np.arange(4,67))]))

    def test_CompareSAR(self):
        
        # Test if patch == False
        area_list = [(np.arange(1,10), np.arange(9,18)), (np.arange(1,10),
                     np.arange(9,18))]

        full_sad = [np.arange(1,40), np.arange(1,60)]

        sar_c = CompareSAR(area_list, ['mete_sar_iter', 'logser-binm'],
                            full_sad)

        # Max area should be 1
        nt.assert_array_equal(np.array([1,1]), np.array([np.max(a) for a in
                            sar_c.a_list]))

        sar_c = CompareSAR(area_list, ['mete_sar_iter', 'logser-binm'],
                            full_sad, max_a=False)

        # Max area should be 9
        nt.assert_array_equal(np.array([9,9]), np.array([np.max(a) for a in
                            sar_c.a_list]))

        # Check species numbers
        bools = [np.all(a == b) for a,b in zip(sar_c.sar_list,
                                          [np.arange(9,18), np.arange(9,18)])]
        self.assertTrue(np.all(bools))

        # Test if patch == True

        rec_sar = np.array(zip(np.arange(1,8), np.arange(4,11)),
                  dtype=[('items', np.float), ('area', np.float)])

        sar_c = CompareSAR([(rec_sar, [])], ['mete_sar_iter'],
                           [np.arange(1,50)], max_a=False, patch=True)

        # check species numbers
        nt.assert_array_equal(np.arange(1,8), sar_c.sar_list[0])

        # Check area numbers
        nt.assert_array_equal(np.arange(4,11), sar_c.a_list[0])

        # check that error is thrown if curve is bad
        self.assertRaises(NameError, CompareSAR, [(rec_sar, [])], ['logser_binm'],
                           [np.arange(1,50)], max_a=False, patch=True)

        # Test compare_curves method
        sar_c = CompareSAR([(rec_sar, [])], ['logser-binm'],
                           [np.arange(1,50)], patch=True)

        # Test with iter_val=False and use_rad=False and all combos
        sar_c.compare_curves()
        sar_c.compare_curves(use_rad=True)
        sar_c.compare_curves(iter_vals=True, use_rad=False)
        sar_c.compare_curves(iter_vals=True, use_rad=True)

    def test_compare_mse(self):
        
        sad_c = CompareSAD(self.sad_data, ['logser', 'lognorm'])

        # Test that mse output has the appropriate formatted data
        mse = sad_c.compare_mse(mse_base='cdf')
        self.assertTrue(len(mse) == 2)
        self.assertTrue(len(mse['lognorm']) == 2 and len(mse['logser']) == 2)
        
        # Test the same thing for a rad base
        mse = sad_c.compare_mse(mse_base='rad')
        self.assertTrue(len(mse) == 2)
        self.assertTrue(len(mse['lognorm']) == 2 and len(mse['logser']) == 2)

        # Test is the the distribution has no cdf MSE is set to NaN
        sad_c = CompareSAD(self.sad_data, ['logser', 'sugihara'])
        mse = sad_c.compare_mse(mse_base='cdf')
        self.assertTrue(np.all(np.isnan(mse['sugihara'])))

        # Test that is works for if base = 'rad'
        sad_c = CompareSAD(self.sad_data, ['logser', 'sugihara'])
        mse = sad_c.compare_mse(mse_base='rad')
        self.assertTrue(type(mse['sugihara'][0] == np.float))
        
        # Test that compare mse works with ssads
        ssad_c = CompareSSAD(self.ssad_data, ['binm', 'tgeo'])
        # Test that mse output has the appropriate formatted data
        mse = ssad_c.compare_mse(mse_base='cdf')
        self.assertTrue(len(mse) == 2)
        self.assertTrue(len(mse['binm']) == 2 and len(mse['tgeo']) == 2)
        
        # Test the same thing for a rad base
        mse = ssad_c.compare_mse(mse_base='rad')
        self.assertTrue(len(mse) == 2)
        self.assertTrue(len(mse['binm']) == 2 and len(mse['tgeo']) == 2)

    def test_compare_rad_cdf(self):

        sad_c = CompareSAD(self.sad_data, ['logser'])
        
        tdist_list = copy.copy(sad_c.dist_list)
        sad_c.dist_list = []

        # Check that rad, cdf work with empty dist list
        rads = sad_c.compare_rads()
        cdfs = sad_c.compare_cdfs()
        self.assertTrue(len(rads) == 1 and len(cdfs) == 1)
        self.assertTrue('observed' in rads and 'observed' in cdfs)
        self.assertTrue(rads == sad_c.rads)
        self.assertTrue(cdfs == sad_c.cdfs)

        # Check that rad, cdf work with something in dist_list
        sad_c.dist_list = tdist_list
        sad_c.rads = None
        sad_c.cdfs = None
        rads = sad_c.compare_rads()
        cdfs = sad_c.compare_cdfs()
        self.assertTrue(len(rads) == 2 and len(cdfs) == 2)
        self.assertTrue('observed' in rads and 'logser' in rads)
        self.assertTrue('observed' in cdfs and 'logser' in cdfs)
        self.assertTrue(rads == sad_c.rads)
        self.assertTrue(cdfs == sad_c.cdfs)

        # Check that if dist doesn't have cdf empty arrays are returned
        sad_c = CompareSAD(self.sad_data, ['logser', 'sugihara'])
        cdfs = sad_c.compare_cdfs()
        self.assertTrue(len(cdfs['sugihara']) == 2)
        self.assertTrue(len(cdfs['sugihara'][0]) == 0 and
                                                 len(cdfs['sugihara'][1]) == 0)

        # check that observed rads are in the right order
        true_vals = np.array([np.all(x == np.array(y)) for x,y in 
                                            zip(rads['observed'], self.sad_data)])

        self.assertTrue(np.all(true_vals))

        # Testing that SED object returns a species list in compare_rads
        patch_sed = [({}, {'a' : np.arange(1,10), 'b' : np.arange(1,20), 'c':
            np.arange(1,30), 'd' : np.arange(1,40)}), ({}, 
            {'a' : np.arange(1,10), 'b' : np.arange(1,20), 'c': 
            np.arange(1,30), 'd' : np.arange(1,40)})]

        patch_sad = [({}, np.arange(1,50), np.repeat('d',20)), 
                    ({}, np.arange(1,50), np.repeat('d',20))]
        patch_ied = [({}, np.arange(4,67), np.repeat('y', 60)),
                     ({}, np.arange(4,67), np.repeat('y', 60))]

        sed_c = CompareSED((patch_sed, patch_ied, patch_sad),
                                               dist_list=['theta'], patch=True)
        
        # Both returns should have a species list
        cdfs = sed_c.compare_rads(return_spp=True)
        rads = sed_c.compare_cdfs(return_spp=True)
        nt.assert_array_equal(np.array(['a', 'b', 'c', 'd', 'a', 'b', 'c',
                                        'd']), np.array(cdfs[1]))
        nt.assert_array_equal(np.array(['a', 'b', 'c', 'd', 'a', 'b', 'c',
                                        'd']), np.array(cdfs[1]))
        nt.assert_array_equal(np.array(['a', 'b', 'c', 'd', 'a', 'b', 'c',
                                        'd']), np.array(rads[1]))


    def test_compare_aic(self):


        # Add another distribution and check the order of the AIC output
        sad_c = CompareSAD(self.sad_data, ['logser', 'most_even', 'nbd_lt'])
        
        aic_out = sad_c.compare_aic(crt=True)
        print aic_out

        # Most even should have the lowest AIC value for the second dataset
        self.assertTrue(aic_out[1][1] == np.min(aic_out[1]))

        aic_m = sad_c.compare_aic_measures(crt=True)

        # Most even should have the a zero delta AIC for the second dataset
        self.assertTrue(aic_m[1][1][1] == np.min(aic_m[1][1]))

        # Most even should have the highest wieght for the second dataset
        self.assertTrue(aic_m[0][1][1] == np.max(aic_m[0][1]))

        # if I don't have any distributions I should get three empty lists for
        # compare_aic_measures
        sad_c = CompareSAD(self.sad_data, [])
        aic_m = sad_c.compare_aic_measures(crt=True)
        self.assertTrue(aic_m == ([],[],[]))

        # If distribution that is passed doesn't have a pmf of pdf, check inf
        # aic values are returned
        sad_c = CompareSAD(self.sad_data, ['logser', 'sugihara'])
        aic_m = sad_c.compare_aic_measures()
        self.assertTrue(aic_m[2][0][1] == np.inf and aic_m[2][1][1] == np.inf)

    def test_compare_LRT(self):

        # Testing compare LRT with logser null model
        sad_c = CompareSAD(self.sad_data, ['nbd_lt'])

        # Is output properly formatted? 
        lrt_out = sad_c.compare_LRT(dist.logser())
        self.assertTrue(len(lrt_out) == 1 and 'logser, nbd_lt' in lrt_out)

    def test_compare_rarity(self):

        #Test compare_rarity

        sad_c = CompareSAD(self.sad_data, ['logser', 'most_even', 'nbd_lt'])
        rare = sad_c.compare_rarity(1)

        # Observed should have 5
        self.assertTrue(rare['observed'][1][0] == 5)
        
        # Most even should have 10 species <= 2
        rare = sad_c.compare_rarity((1,2))
        self.assertTrue(rare['observed'][1][0] == 5)
        self.assertTrue(rare['most_even'][2][1] == 10)

    def test_compare_moments(self):

        # Test the compare_moments output is formatted correctly
        sad_c = CompareSAD(self.sad_data, ['logser', 'nbd_lt'])
        mom = sad_c.compare_moments()
        self.assertTrue(len(mom) == 3)

        # Test that observed and all distributions are considered
        lengths = np.array([len(mom[x]) for x in mom.iterkeys()])

        self.assertTrue(np.array_equal(lengths, np.repeat(3, 3)))

    def test_summary(self):

        # Test that summary output is correct
        # Test is there are no dists in dist_list
        sad_c = CompareSAD(self.sad_data, [])
        sumry = sad_c.summary()
        # Test that there is only observed in summary dict
        self.assertTrue(len(sumry) == 1 and 'observed' in sumry)
        
        # Test if we have two distributions but one doesn't have a cdf
        sad_c = CompareSAD(self.sad_data, ['logser', 'sugihara'])
        smry = sad_c.summary()
        self.assertTrue(len(smry) == 3)

        # Logseries dict and sugihara dict should have 9 kw
        self.assertTrue(len(smry['logser']) == 9 and len(smry['sugihara']) ==
                                                                             9)

        # AIC values for sugihara should in inf
        self.assertTrue(np.all(smry['sugihara']['aic'] == np.array([np.inf,
                                                                     np.inf])))
        # IED should be able to call summary
        ied_data = [(np.arange(10,100), np.arange(1,40)), (np.arange(1,20),
                                                             np.arange(40,60))]
        ied_c = CompareIED(ied_data, dist_list=['psi'])
        smry = ied_c.summary()
        self.assertTrue(smry['observed']['balls'] == [4905, 190])

    def test_nll(self):
        
        # Test against R result: sum(dnorm(c(1,2,3,4,5), log=TRUE))
        R_res = 32.09469
        test_vals = stats.norm.pdf((1,2,3,4,5))
        lglk = nll([test_vals])[0]
        self.assertTrue(R_res == np.round(lglk, decimals=5))

    def test_empirical_cdf(self):
        
        #Test against R's ecdf function
        test_data = [1,1,1,1,2,3,4,5,6,6]
        R_res = [.4,.4,.4,.4,.5,.6,.7,.8,1,1]
        res = empirical_cdf(test_data)
        self.assertTrue(np.array_equal(R_res, res))
        
        test_data = [3,3,3,3]
        R_res = [1,1,1,1]
        res = empirical_cdf(test_data)
        self.assertTrue(np.array_equal(R_res, res))

    def test_aic(self):
        
        # Test that passing either a pmf of nll gives the same result
        test_vals = stats.norm.pdf((1,2,3,4,5,6,7,8))
        aic1 = aic([test_vals], 2, loglik=False)
        aic2 = aic(nll([test_vals]), 2, loglik=True)

        self.assertTrue(aic1[0] == aic2[0])
        # Expected AIC for test_vals
        expected = 6.837877066 # Calculated by hand
        self.assertTrue(np.round(aic1[0], decimals=9), expected)

        test_vals = stats.gamma.pdf((1,1,1,4,5,7,12),2)
        aic1 = aic([test_vals], 2, loglik=False)
        expected = 51.146902
        self.assertTrue(np.round(aic1[0], decimals=6), expected)

    def test_aicc(self):
        
        # Test that passing either a pmf of nll gives the same result
        test_vals = stats.norm.pdf((1,2,3,4,5,6,7,8))
        aic1 = aicc([test_vals], 2, loglik=False)
        aic2 = aicc(nll([test_vals]), 2, 8, loglik=True)

        self.assertTrue(aic1[0] == aic2[0])

        # Test that aicc gives the correct values
        expected = 225.10302
        self.assertTrue(expected == np.round(aic1[0], decimals=5))

        # Test Assertion error is thrown if no n param
        self.assertRaises(AssertionError, aicc, 56, 2)


    def test_aic_weights(self):
        
        vals = [1,1,1,2,3,4,7,23,78]
        aic_vals = aicc([stats.norm.pdf(vals, scale=100), stats.norm.pdf(vals,
                                                                    scale=99)],
                                                            [2,2],loglik=False)
        aicw, delta_aic = aic_weights(aic_vals)
        pred = np.array([ 0.47909787,  0.52090213])
        self.assertTrue(np.array_equal(np.round(aicw, decimals=8), pred))
         

    def test_ks_two_sample(self):
        # Unittested in scipy, testing that this function works

        d, p = ks_two_sample([1,1,2,3,4,5,6,12], [1,2,3,4,5,5,5,5,5,7,8,9])

    def test_likelihood_ratio(self):
        
        # Test against what the lrtest() R function returns
        model1 = 158.0494
        model0 = 139.806
        R_chisquare = 36.4868
        R_p = 1.537e-09

        pred_chi, pred_p = likelihood_ratio(model0, model1, 1)[0]

        self.assertTrue(np.round(pred_chi, decimals=4) == R_chisquare)
        pred_p = np.round(pred_p, decimals=12) 
        self.assertTrue(pred_p == R_p)


    def test_variance(self):
        
        # Test that I get back the correct values
        data = [[0,1,2,3,4,45,18,56,24,56], [1,1,1,1,56,78,23,23]]
        expt = []
        expt.append(np.var(data[0], ddof=1))
        expt.append(np.var(data[1], ddof=1))
        resulting_vals = variance(data)
        self.assertTrue(np.array_equal(np.array(expt),
                                                    np.array(resulting_vals)))
        # Using np.var which is optimized and unittested

    def test_skew(self):
        
        # Using the scipy.stats definition which is optimized and unittested
        data = [[0,1,2,3,4,45,18,56,24,56], [1,1,1,1,56,78,23,23]]
        expt = []
        expt.append(stats.skew(data[0]))
        expt.append(stats.skew(data[1]))
        resulting_vals = skew(data)
        self.assertTrue(np.array_equal(np.array(expt),
                                                    np.array(resulting_vals)))

    def test_kurtosis(self):
        
        # Using the scipy.stats definition which is optimized and unittested
        data = [[0,1,2,3,4,45,18,56,24,56], [1,1,1,1,56,78,23,23]]
        expt = []
        expt.append(stats.kurtosis(data[0]))
        expt.append(stats.kurtosis(data[1]))
        resulting_vals = kurtosis(data)
        self.assertTrue(np.array_equal(np.array(expt),
                                                    np.array(resulting_vals)))

    def test_mean_square_error(self):
        
        # Test against R mse function
        pred = np.arange(1,9)
        obs = np.arange(7, 15)

        comp_val = 36
        pred = mean_squared_error(pred, obs)
        self.assertEqual(pred, comp_val)

    def test_bootstrap_moment(self):
        
        data1 = np.arange(1, 31)
        data2 = np.arange(20, 50)
        # Test the return is empty if wrong keyword is given
        bs_vals = bootstrap_moment(data1, data2, ['men', 'vaiance',
                            'sew', 'kurtoss'], num_samp=100)

        self.assertTrue(len(bs_vals) == 0)
        
        # Test bootstrap moment against William Rice's (UCSB) bootstrap 
        # programs in Statistics 101.  Just testing the mean, but the
        # implementation is the same for all of them
        test_ci = np.array([-23.4, -14.6])

        bs_vals = bootstrap_moment(data1, data2, ['mean', 'variance',
                            'skew', 'kurtosis'], num_samp=50000)

        # Check that Bill Rice's and our 95% CIs match
        self.assertTrue(np.array_equal(test_ci, np.round(bs_vals['mean'][1],
            decimals=1)))
        
        # Check that the deltas match
        self.assertTrue(-19 == bs_vals["mean"][0])

        # Check that the length is right
        self.assertTrue(len(bs_vals) == 4)

if __name__ == '__main__':
    unittest.main()
