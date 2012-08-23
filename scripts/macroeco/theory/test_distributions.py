#!/usr/bin/python

#Testing predict sad
#NOTE: Unit test must be called with cwd macroeco/code/macroeco

import unittest
from distributions import *
import numpy as np
import os 

gcwd = os.getcwd
pd = os.path.dirname
jp = os.path.join

class TestDistributions(unittest.TestCase):
    '''Test the functions within distributions.py'''

    def setUp(self):
        
        self.abund_list = []
        self.abund_list.append(np.array([1,1,2,3,4,4,7,12]))
        self.abund_list.append(np.array([1,1,1,1,1,2]))
        self.abund_list.append(np.array([3,2,2,1,5,6,4,5,7]))
        self.abund_list.append(np.array([1,1,4,5,7,89,100]))
        self.R = [0.8453224, 1.951546, 1.040038, 0.3102524]

    def test_plognorm_pmf(self):
        for i, sad in enumerate(self.abund_list):
            log_abund = np.log(sad)
            mu = np.mean(log_abund)
            var = np.var(log_abund, ddof=1)
            sigma = var**0.5
            pmf = sum(plognorm().pmf(sad, mu, sigma))
            diff1 = np.round(self.R[i], decimals=5) - np.round(pmf, decimals=5)
            self.assertTrue(abs(diff1) == 0)
        self.assertTrue(sum(np.round(plognorm().pmf([1,2,3,4,5], -3,-3),\
                            decimals=3)) == 0)

    def test_plognorm_cdf(self):
        pass

    """def test_mete_lgsr_pmf(self):
        self.assertRaises(AssertionError, mete_lgsr_pmf, 1, 45, 45)
        self.assertRaises(AssertionError, mete_lgsr_pmf, 1, 234, 67)
        self.assertRaises(AssertionError, mete_lgsr_pmf, 1, 34, 0)
        pmf = mete_lgsr_pmf(np.arange(1, 568), 34, 567)
        self.assertTrue(len(pmf) == 567)
        #Testing that values equal values from John's book (Harte 2011)
        pmf, x = mete_lgsr_pmf(1, 4, 4 * 4, param_out=True)
        self.assertTrue(np.round(x['beta'], decimals=4) == 0.0459)
        pmf, x = mete_lgsr_pmf(1, 4, 2**4 * 4,param_out=True)
        self.assertTrue(np.round(x['beta'], decimals=5) == -0.00884)
        pmf, x = mete_lgsr_pmf(1, 4, 2**8 * 4, param_out=True)
        self.assertTrue(np.round(x['beta'], decimals=5) == -0.00161)
        pmf, x = mete_lgsr_pmf(1, 16, 2**8 * 16, param_out=True)
        self.assertTrue(np.round(x['beta'], decimals=6) == 0.000413)
        pmf, x = mete_lgsr_pmf(1, 64, 2**12 * 64, param_out=True)
        self.assertTrue(np.round(x['beta'], decimals=7) == 0.0000228)
        
    def test_mete_lgsr_approx_pmf(self):
        self.assertRaises(AssertionError, mete_lgsr_approx_pmf, 1, 45, 45)
        self.assertRaises(AssertionError, mete_lgsr_approx_pmf, 1, 234, 67)
        self.assertRaises(AssertionError, mete_lgsr_approx_pmf, 1, 34, 0)
        #Testing that values = values from John's book (Harte 2011)
        pmf, x = mete_lgsr_approx_pmf(1, 4, 4 * 4, param_out=True)
        self.assertTrue(np.round(x['beta'], decimals=3) == 0.116)
        pmf, x = mete_lgsr_approx_pmf(1, 4, 2**4 * 4, param_out=True)
        self.assertTrue(np.round(x['beta'], decimals=4) == 0.0148)
        pmf, x = mete_lgsr_approx_pmf(1, 4, 2**8 * 4, param_out=True)
        self.assertTrue(np.round(x['beta'], decimals=6) == 0.000516)
        pmf, x = mete_lgsr_approx_pmf(1, 16, 2**8 * 16, param_out=True)
        self.assertTrue(np.round(x['beta'], decimals=6) == 0.000516)
        pmf, x = mete_lgsr_approx_pmf(1, 64, 2**12 * 64, param_out=True)
        self.assertTrue(np.round(x['beta'], decimals=7) == 0.0000229 or\
                        np.round(x['beta'], decimals=7) == 0.0000228)

    def test_lgsr_pmf(self):
        self.assertRaises(AssertionError, lgsr_pmf, 1, 45, 45)
        self.assertRaises(AssertionError, lgsr_pmf, 1, 234, 67)
        self.assertRaises(AssertionError, lgsr_pmf, 1, 34, 0)
        #Test against value in Fisher's paper Fisher et al. 1943
        pmf, x = lgsr_pmf(1, 240, 15609, param_out=True)
        self.assertTrue(np.round(x['x'], decimals=4) == 0.9974)

    #Kind of a lame test...
    def test_geo_pmf(self):
        geo = geo_pmf(np.arange(1, 101), 23, 100)
        self.assertTrue(np.array_equal(neg_binom_pmf(np.arange(1, 101), 23,\
                        100, 1), geo))
        self.assertTrue(len(geo) == 100)

    def test_sugihara_rank_abund(self):
        #Testing against A. J. Baczkowski 1997 paper values
        #The function does not agree with Backzkowskis values for S = 4
        ps = sugihara_rank_abund(3, 100, sample_size=50000) / 100
        self.assertTrue(np.round(ps[0], decimals=2) == 0.65)
        self.assertTrue(np.round(ps[1], decimals=2) == 0.23)
        self.assertTrue(np.round(ps[2], decimals=2) == 0.10)"""

if __name__ == '__main__':
    unittest.main()
