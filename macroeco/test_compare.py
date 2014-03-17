#!/usr/bin/python
"""
Tests for compare module

"""
from __future__ import division

from numpy.testing import (TestCase, assert_equal, assert_array_equal, 
                           assert_almost_equal, assert_array_almost_equal,
                           assert_allclose, assert_, assert_raises)

from macroeco.compare import *
import numpy as np
import scipy.stats as stats
import copy
import macroeco.distributions as dist
import numpy.testing as nt

class TestCompare(TestCase):
    '''Test Methods in compare.py'''

    def test_nll(self):
        
        # Test against R result: sum(dnorm(c(1,2,3,4,5), log=TRUE))
        R_res = 32.09469
        test_vals = stats.norm.pdf((1,2,3,4,5))
        lglk = get_nll(test_vals)
        assert_equal(R_res, np.round(lglk, decimals=5))

    def test_empirical_cdf(self):
        
        #Test against R's ecdf function

        # Test Case 1
        test_data = [1,1,1,1,2,3,4,5,6,6]
        R_res = [.4,.4,.4,.4,.5,.6,.7,.8,1,1]
        res = get_empirical_cdf(test_data)
        assert_array_equal(R_res, res)
        
        # Test Case 2
        test_data = [3,3,3,3]
        R_res = [1,1,1,1]
        res = get_empirical_cdf(test_data)
        assert_array_equal(R_res, res)

    def test_aic(self):
        
        test_vals = stats.norm.pdf((1,2,3,4,5,6,7,8))
        aic1 = get_AIC(test_vals, (1,1))
        expected = 222.703016531  # Calculated by hand
        assert_equal(np.round(aic1, decimals=9), expected)


        test_vals = stats.gamma.pdf((1,1,1,4,5,7,12),2)
        aic1 = get_AIC(test_vals, (1,1))
        expected = 51.146902
        assert_equal(np.round(aic1, decimals=6), expected)

    def test_aicc(self):
        
        # Test values
        test_vals = stats.norm.pdf((1,2,3,4,5,6,7,8))
        aic1 = get_AICC(test_vals, (1,1))

        # Test that aicc gives the correct values
        expected = 225.10302
        assert_equal(expected, np.round(aic1, decimals=5))

    def test_aic_weights(self):
        
        # Test values
        vals = [1,1,1,2,3,4,7,23,78]
        values = [stats.norm.pdf(vals, scale=100), stats.norm.pdf(vals,
                                                                    scale=99)]

        aic_vals = [get_AICC(tval, 1) for tval in values]
        aicw, delta_aic = get_AIC_weights(aic_vals)
        pred = np.array([ 0.47909787,  0.52090213])
        assert_array_almost_equal(aicw, pred)

    def test_gen_loss_function(self):
        
        # Test absolute value loss function
        loss_fxn = 'np.abs(obs - pred)'
        loss = gen_loss_function(loss_fxn)

        obs = np.random.randint(3, 59, 100)
        pred = np.random.randint(3, 59, 100)
        test_loss = np.sum(np.abs(obs - pred))

        pred_loss = loss.total_loss(obs, pred)
        assert_equal(pred_loss, test_loss)
        
        # Test sum of squares loss function
        test_loss = np.sum((obs - pred)**2)
        pred_loss = get_sum_of_squares(obs, pred)
        assert_equal(test_loss, pred_loss)

        # Test MSE loss function 
        loss_fxn = 'np.abs(obs - pred) / len(obs)'
        loss = gen_loss_function(loss_fxn)

        test_loss = np.sum(np.abs(obs - pred) / len(obs))
        pred_loss = loss.total_loss(obs, pred)
        assert_equal(test_loss, pred_loss)
        


#         
#
#    def test_ks_two_sample(self):
#        # Unittested in scipy, testing that this function works
#
#        d, p = ks_two_sample([1,1,2,3,4,5,6,12], [1,2,3,4,5,5,5,5,5,7,8,9])
#
#    def test_likelihood_ratio(self):
#        
#        # Test against what the lrtest() R function returns
#        model1 = 158.0494
#        model0 = 139.806
#        R_chisquare = 36.4868
#        R_p = 1.537e-09
#
#        pred_chi, pred_p = likelihood_ratio(model0, model1, 1)[0]
#
#        self.assertTrue(np.round(pred_chi, decimals=4) == R_chisquare)
#        pred_p = np.round(pred_p, decimals=12) 
#        self.assertTrue(pred_p == R_p)
#
#
#    def test_variance(self):
#        
#        # Test that I get back the correct values
#        data = [[0,1,2,3,4,45,18,56,24,56], [1,1,1,1,56,78,23,23]]
#        expt = []
#        expt.append(np.var(data[0], ddof=1))
#        expt.append(np.var(data[1], ddof=1))
#        resulting_vals = variance(data)
#        self.assertTrue(np.array_equal(np.array(expt),
#                                                    np.array(resulting_vals)))
#        # Using np.var which is optimized and unittested
#
#    def test_skew(self):
#        
#        # Using the scipy.stats definition which is optimized and unittested
#        data = [[0,1,2,3,4,45,18,56,24,56], [1,1,1,1,56,78,23,23]]
#        expt = []
#        expt.append(stats.skew(data[0]))
#        expt.append(stats.skew(data[1]))
#        resulting_vals = skew(data)
#        self.assertTrue(np.array_equal(np.array(expt),
#                                                    np.array(resulting_vals)))
#
#    def test_kurtosis(self):
#        
#        # Using the scipy.stats definition which is optimized and unittested
#        data = [[0,1,2,3,4,45,18,56,24,56], [1,1,1,1,56,78,23,23]]
#        expt = []
#        expt.append(stats.kurtosis(data[0]))
#        expt.append(stats.kurtosis(data[1]))
#        resulting_vals = kurtosis(data)
#        self.assertTrue(np.array_equal(np.array(expt),
#                                                    np.array(resulting_vals)))
#
#    def test_mean_square_error(self):
#        
#        # Test against R mse function
#        pred = np.arange(1,9)
#        obs = np.arange(7, 15)
#
#        comp_val = 36
#        pred = mean_squared_error(pred, obs)
#        self.assertEqual(pred, comp_val)
#
#    def test_bootstrap_moment(self):
#        
#        data1 = np.arange(1, 31)
#        data2 = np.arange(20, 50)
#        # Test the return is empty if wrong keyword is given
#        bs_vals = bootstrap_moment(data1, data2, ['men', 'vaiance',
#                            'sew', 'kurtoss'], num_samp=100)
#
#        self.assertTrue(len(bs_vals) == 0)
#        
#        # Test bootstrap moment against William Rice's (UCSB) bootstrap 
#        # programs in Statistics 101.  Just testing the mean, but the
#        # implementation is the same for all of them
#        test_ci = np.array([-23.4, -14.6])
#
#        bs_vals = bootstrap_moment(data1, data2, ['mean', 'variance',
#                            'skew', 'kurtosis'], num_samp=50000)
#
#        # Check that Bill Rice's and our 95% CIs match
#        self.assertTrue(np.array_equal(test_ci, np.round(bs_vals['mean'][1],
#            decimals=1)))
#        
#        # Check that the deltas match
#        self.assertTrue(-19 == bs_vals["mean"][0])
#
#        # Check that the length is right
#        self.assertTrue(len(bs_vals) == 4)
#
