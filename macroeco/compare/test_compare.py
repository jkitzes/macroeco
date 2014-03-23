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
import numpy.testing as nt


class TestCompare(TestCase):
    '''Test Methods in compare.py'''

    def test_nll(self):
        
        # Test against R result: sum(dnorm(c(1,2,3,4,5), log=TRUE))
        R_res = 32.09469
        test_vals = stats.norm.pdf((1, 2, 3, 4, 5))
        lglk = nll(test_vals)
        assert_equal(R_res, np.round(lglk, decimals=5))

    def test_empirical_cdf(self):
        
        #Test against R's ecdf function

        # Test Case 1
        test_data = [1, 1, 1, 1, 2, 3, 4, 5, 6, 6]
        R_res = [.4, .4, .4, .4, .5, .6, .7, .8, 1, 1]
        res = empirical_cdf(test_data)
        assert_array_equal(R_res, res)
        
        # Test Case 2
        test_data = [3, 3, 3, 3]
        R_res = [1, 1, 1, 1]
        res = empirical_cdf(test_data)
        assert_array_equal(R_res, res)

    def test_aic(self):
        
        test_vals = stats.norm.pdf((1, 2, 3, 4, 5, 6, 7, 8))
        aic1 = AIC(test_vals, (1, 1))
        expected = 222.703016531  # Calculated by hand
        assert_equal(np.round(aic1, decimals=9), expected)

        test_vals = stats.gamma.pdf((1, 1, 1, 4, 5, 7, 12), 2)
        aic1 = AIC(test_vals, (1, 1))
        expected = 51.146902
        assert_equal(np.round(aic1, decimals=6), expected)

    def test_aicc(self):
        
        # Test values
        test_vals = stats.norm.pdf((1, 2, 3, 4, 5, 6, 7, 8))
        aic1 = AICC(test_vals, (1, 1))

        # Test that aicc gives the correct values
        expected = 225.10302
        assert_equal(expected, np.round(aic1, decimals=5))

    def test_aic_weights(self):
        
        # Test values
        vals = [1, 1, 1, 2, 3, 4, 7, 23, 78]
        values = [stats.norm.pdf(vals, scale=100), stats.norm.pdf(vals,
                                                                    scale=99)]

        aic_vals = [AICC(tval, 1) for tval in values]
        aicw, delta_aic = AIC_weights(aic_vals)
        pred = np.array([0.47909787, 0.52090213])
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
        test_loss = np.sum((obs - pred) ** 2)
        pred_loss = sum_of_squares(obs, pred)
        assert_equal(test_loss, pred_loss)

        # Test MSE loss function
        loss_fxn = 'np.abs(obs - pred) / len(obs)'
        loss = gen_loss_function(loss_fxn)

        test_loss = np.sum(np.abs(obs - pred) / len(obs))
        pred_loss = loss.total_loss(obs, pred)
        assert_equal(test_loss, pred_loss)

    def test_r_squared(self):

        # Already unittested in scipy. Checking for functionaliity
        test_data = np.random.randint(5, 100, 100)
        rsq = r_squared(test_data, test_data)
        assert_equal(rsq, 1)

    def test_chi_squared(self):

        # Compare two distributions
        # Chi squared function itself is already unittested in scipy

        bin_max = 16
        p = 0.99
        dist1 = stats.logser(p=p).rvs(100)
        dist2 = stats.logser(p=p).rvs(100)

        bin1 = bin_data(dist1, np.max(bin_max))[0]
        bin2 = bin_data(dist2, np.max(bin_max))[0]

        res = chi_squared([bin1, bin2])

        # Check three distributions
        dist3 = stats.logser(p=p).rvs(100)
        bin3 = bin_data(dist3, np.max(bin_max))[0]

        res = chi_squared([bin1, bin2, bin3])

        # Check error is thrown with only one dist
        assert_raises(AssertionError, chi_squared, [bin1])

        # Check error is thrown if bins are different lengths
        assert_raises(AssertionError, chi_squared, [bin1, bin2[:-1]])

    def test_bin_data(self):

        # Test against R's vegan prestonfit: prestonfit(data, tiesplit=FALSE)
        # Note that vegan drops the bins with 0 values

        data = np.array([1, 1, 1, 1, 2, 2, 4, 4, 8, 16, 17.1, 89])
        vegan = np.array([4, 2, 2, 1, 1, 1, 0, 1], dtype=np.float)
        test_res = bin_data(data, max(data))[0]
        assert_array_equal(test_res, vegan)

        data = np.array([1, 1, 1, 1, 4, 5, 6, 7, 12, 34, 56])
        vegan = np.array([4, 0, 1, 3, 1, 0, 2], dtype=np.float)
        test_res = bin_data(data, max(data))[0]
        assert_array_equal(test_res, vegan)

        # Test boundary condition
        data = np.array([1, 2])
        vegan = np.array([1, 1], dtype=np.float)
        test_res = bin_data(data, max(data))[0]
        assert_array_equal(test_res, vegan)

        data = np.array([1, 1, 1])
        vegan = np.array([3], dtype=np.float)
        test_res = bin_data(data, max(data))[0]
        assert_array_equal(test_res, vegan)

        data = np.array([1, 2, 3])
        vegan = np.array([1, 1, 1], dtype=np.float)
        test_res = bin_data(data, max(data))[0]
        assert_array_equal(test_res, vegan)

    def test_lrt(self):
        
        # Test against what the lrtest() R function returns
        model1 = 158.0494
        model0 = 139.806
        R_chisquare = 36.4868
        R_p = 1.537e-09

        pred_chi, pred_p = lrt(model1, model0, 1)

        assert_almost_equal(pred_chi, R_chisquare)
        assert_almost_equal(pred_p, R_p)

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
