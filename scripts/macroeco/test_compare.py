#!/usr/bin/python

#Testing Compare Module

import unittest
from macroeco.compare import *
import numpy as np
import scipy.stats as stats

class TestCompare(unittest.TestCase):
    '''Test classes and methods in compare.py'''

    def setup(self):
        pass
        

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
        pass

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


if __name__ == '__main__':
    unittest.main()
