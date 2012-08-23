#!/usr/bin/python

#Testing predict sad

'''
References for Unit Tests
-------------------------
Baczkowski, A. J. 1997. Sugihara's Sequential Model for Species Abundances.
Internal Review.

Magurran, A. E. 1988. Ecological Diversity and Its Measuremnt. Princeton
University Press.

May, R. M. 1975. Patterns of species abundance and diversity. In Ecology and
Evolution of Communities (eds M. L. Cody and J. M. Diamond), Harvard University
Press.
'''

import unittest
from sad_distr import *
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
        #R output from function dpolono in VGAM package
        self.R = [0.8453224, 1.951546, 1.040038, 0.3102524]
        #Simulated values in A. J. Baczkowski 1997 paper 
        self.sugi = [.4761, .1962, .1180, .0751, .0499, .0337, .0226, .0148,\
                     .0090, .0047]

    def test_lgsr_pmf(self):
        self.assertRaises(AssertionError, lgsr().pmf, 1, 45, 45)
        self.assertRaises(AssertionError, lgsr().pmf, 1, 234, 67)
        self.assertRaises(AssertionError, lgsr().pmf, 1, 34, 0)
        #Test against value in Fisher's paper Fisher et al. 1943
        pmf, x = lgsr().pmf(1, 240, 15609, param_out=True)
        self.assertTrue(np.round(x['x'], decimals=4) == 0.9974)

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

    def test_mete_lgsr(self):
        self.assertRaises(AssertionError, mete_lgsr().pmf, 1, 45, 45)
        self.assertRaises(AssertionError, mete_lgsr().pmf, 1, 234, 67)
        self.assertRaises(AssertionError, mete_lgsr().pmf, 1, 34, 0)
        pmf = mete_lgsr().pmf(np.arange(1, 568), 34, 567)
        self.assertTrue(len(pmf) == 567)
        #Testing that values equal values from John's book (Harte 2011)
        pmf, x = mete_lgsr().pmf(1, 4, 4 * 4, param_out=True)
        self.assertTrue(np.round(x['beta'], decimals=4) == 0.0459)
        pmf, x = mete_lgsr().pmf(1, 4, 2**4 * 4,param_out=True)
        self.assertTrue(np.round(x['beta'], decimals=5) == -0.00884)
        pmf, x = mete_lgsr().pmf(1, 4, 2**8 * 4, param_out=True)
        self.assertTrue(np.round(x['beta'], decimals=5) == -0.00161)
        pmf, x = mete_lgsr().pmf(1, 16, 2**8 * 16, param_out=True)
        self.assertTrue(np.round(x['beta'], decimals=6) == 0.000413)
        pmf, x = mete_lgsr().pmf(1, 64, 2**12 * 64, param_out=True)
        self.assertTrue(np.round(x['beta'], decimals=7) == 0.0000228)
        
    def test_mete_lgsr_approx_pmf(self):
        self.assertRaises(AssertionError, mete_lgsr_approx().pmf, 1, 45, 45)
        self.assertRaises(AssertionError, mete_lgsr_approx().pmf, 1, 234, 67)
        self.assertRaises(AssertionError, mete_lgsr_approx().pmf, 1, 34, 0)
        #Testing that values = values from John's book (Harte 2011)
        pmf, x = mete_lgsr_approx().pmf(1, 4, 4 * 4, param_out=True)
        self.assertTrue(np.round(x['beta'], decimals=3) == 0.116)
        pmf, x = mete_lgsr_approx().pmf(1, 4, 2**4 * 4, param_out=True)
        self.assertTrue(np.round(x['beta'], decimals=4) == 0.0148)
        pmf, x = mete_lgsr_approx().pmf(1, 4, 2**8 * 4, param_out=True)
        self.assertTrue(np.round(x['beta'], decimals=6) == 0.000516)
        pmf, x = mete_lgsr_approx().pmf(1, 16, 2**8 * 16, param_out=True)
        self.assertTrue(np.round(x['beta'], decimals=6) == 0.000516)
        pmf, x = mete_lgsr_approx().pmf(1, 64, 2**12 * 64, param_out=True)
        self.assertTrue(np.round(x['beta'], decimals=7) == 0.0000229 or\
                        np.round(x['beta'], decimals=7) == 0.0000228)

    def test_sugihara_rank_abund(self):
        #Testing against A. J. Baczkowski 1997 paper values
        #Passes with 2% error regularly. Being conservative with 5% error
        ps = sugihara().rad(10, 400, sample_size=20000) / 400
        error = .05 * ps
        diff = np.array(self.sugi) - ps
        ind = np.abs(diff) <= error
        self.assertTrue(np.all(ind))

    def test_geo_series(self):
        #Using data from Magurran (1998)
        obs_sad = [370,210,120,66,35,31,15,9,3,2,1]
        mag_pred_sad = [387.6, 213.8, 117.8, 64.5, 35.5, 19.8, 10.7, 6.2, 3.3,\
                        1.8, 1.0]
        gk = geo_series().fit_k(obs_sad)
        self.assertTrue(np.round(gk, decimals=3) == 0.449)
        geo_sad = np.round(geo_series().rad(len(obs_sad), sum(obs_sad), .449),\
                           decimals=1)
        diff = np.floor(np.array(mag_pred_sad)) - np.floor(geo_sad)
        self.assertTrue(np.all(diff == 0))

    def test_broken_stick(self):
        #Using data from Magurran (1998)
        self.assertRaises(AssertionError, broken_stick().pmf, 112, 12, 111)
        obs_sad = [103,115,13,2,67,36,51,8,6,61,10,21,7,65,4,49,92,37,16,6,23,\
                   9,2,6,5,4,1,3,1,9,2]
        expt = [1.077, 1.040, 1.004, .970, .937, .904, .873, .843, .814,\
                    .786, .759, .732, .707, .683, .659, .636]
        S = len(obs_sad)
        N = sum(obs_sad)
        bs = np.round(broken_stick().pmf(np.arange(1, 17), S, N) * S,\
                      decimals=3)
        diff = np.array(expt) - bs
        self.assertTrue(np.all(diff == 0))

    def test_lognorm_pmf(self):
        r_output = [0.1210, .0806, .0601, 0.0476, 0.0391, .0331,  0.0285,\
                    0.0249, 0.0221, 0.0197]
        lnorm = np.round(lognorm().pmf(np.arange(1,11), 2, 2), decimals=4)
        diff = r_output - lnorm
        self.assertTrue(np.all(diff == 0))



if __name__ == '__main__':
    unittest.main()
