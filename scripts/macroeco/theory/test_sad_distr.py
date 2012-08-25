#!/usr/bin/python

#Testing sad distributions

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

class TestSADDistributions(unittest.TestCase):
    '''Test the functions within sad_distr.py'''

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
        self.assertRaises(AssertionError, lgsr(S=45, N=45).pmf, 1)
        self.assertRaises(AssertionError, lgsr(S=234, N=67).pmf, 1)
        self.assertRaises(AssertionError, lgsr(S=34, N=0).pmf, 1)
        #Test against value in Fisher's paper Fisher et al. 1943
        pmf, x = lgsr(S=240, N=15609).pmf(1, param_out=True)
        self.assertTrue(np.round(x['x'], decimals=4) == 0.9974)
        cdf = np.round(lgsr(S=45, N=1200).cdf(1200), decimals=1)
        self.assertTrue(cdf == 1)
    
    #Should test against ethans psolver
    def test_plognorm(self):
        for i, sad in enumerate(self.abund_list):
            log_abund = np.log(sad)
            mu = np.mean(log_abund)
            var = np.var(log_abund, ddof=1)
            sigma = var**0.5
            pmf = sum(plognorm(mu=mu, sigma=sigma).pmf(sad))
            diff1 = np.round(self.R[i], decimals=5) - np.round(pmf, decimals=5)
            self.assertTrue(abs(diff1) == 0)
        self.assertTrue(sum(np.round(plognorm(mu=-3,sigma=-3).\
                                     pmf([1,2,3,4,5]), decimals=3)) == 0) 
        plognorm().fit(self.abund_list[0])
        N = sum(self.abund_list[0]); S = len(self.abund_list[0])
        plognorm(S=S, N=N, mu=2, sigma=2).cdf(5)

    def test_trun_lognorm(self):
        N = sum(self.abund_list[0]); S = len(self.abund_list[0])
        trun_plognorm(S=S, N=N, mu=2, sigma=2).cdf(5)
        trun_plognorm(mu=2, sigma=2).pmf([2,3,4,5,23])
        trun_plognorm().fit(self.abund_list[0])
        trun_plognorm(mu=10, sigma=1).cdf(45)
        EW_fit = {'mu' :.90, 'sigma' : 2.18}
        sad = [1,1,1,1,2,3,5,6,12,13,15,23,45,67,112]
        pred = trun_plognorm().fit(sad)
        self.assertTrue(EW_fit['mu'] == np.round(pred['mu'], decimals=2))
        self.assertTrue(EW_fit['sigma'] == np.round(pred['sigma'], decimals=2))

    def test_mete_lgsr(self):
        self.assertRaises(AssertionError, mete_lgsr(S=45, N=45).pmf, 1)
        self.assertRaises(AssertionError, mete_lgsr(S=234, N=67).pmf, 1)
        self.assertRaises(AssertionError, mete_lgsr(S=34, N=0).pmf, 1)
        pmf = mete_lgsr(S=34, N=567).pmf(np.arange(1, 568))
        self.assertTrue(len(pmf) == 567)
        #Testing that values equal values from John's book (Harte 2011)
        pmf, x = mete_lgsr(S=4, N=4 * 4).pmf(1, param_out=True)
        self.assertTrue(np.round(x['beta'], decimals=4) == 0.0459)
        pmf, x = mete_lgsr(S=4, N=2**4 * 4).pmf(1 ,param_out=True)
        self.assertTrue(np.round(x['beta'], decimals=5) == -0.00884)
        pmf, x = mete_lgsr(S=4, N=2**8 * 4).pmf(1, param_out=True)
        self.assertTrue(np.round(x['beta'], decimals=5) == -0.00161)
        pmf, x = mete_lgsr(S=16, N=2**8 * 16).pmf(1, param_out=True)
        self.assertTrue(np.round(x['beta'], decimals=6) == 0.000413)
        pmf, x = mete_lgsr(S=64, N=2**12 * 64).pmf(1, param_out=True)
        self.assertTrue(np.round(x['beta'], decimals=7) == 0.0000228)
        mete_lgsr(S=64, N=1000).rad()
        mete_lgsr(S=64, N=1000).cdf((1,1,2,4,5,7,12))

 
    def test_mete_lgsr_approx_pmf(self):
        self.assertRaises(AssertionError, mete_lgsr_approx(S=45, N=45).pmf, 1)
        self.assertRaises(AssertionError, mete_lgsr_approx(S=234, N=67).pmf, 1)
        self.assertRaises(AssertionError, mete_lgsr_approx(S=34, N=0).pmf, 1)
        #Testing that values = values from John's book (Harte 2011)
        pmf, x = mete_lgsr_approx(S=4, N=4 * 4).pmf(1, param_out=True)
        self.assertTrue(np.round(x['beta'], decimals=3) == 0.116)
        pmf, x = mete_lgsr_approx(S=4, N=2**4 * 4).pmf(1, param_out=True)
        self.assertTrue(np.round(x['beta'], decimals=4) == 0.0148)
        pmf, x = mete_lgsr_approx(S=4, N=2**8 * 4).pmf(1, param_out=True)
        self.assertTrue(np.round(x['beta'], decimals=6) == 0.000516)
        pmf, x = mete_lgsr_approx(S=16, N=2**8 * 16).pmf(1, param_out=True)
        self.assertTrue(np.round(x['beta'], decimals=6) == 0.000516)
        pmf, x = mete_lgsr_approx(S=64, N=2**12 * 64).pmf(1, param_out=True)
        self.assertTrue(np.round(x['beta'], decimals=7) == 0.0000229 or\
                        np.round(x['beta'], decimals=7) == 0.0000228)
        mete_lgsr_approx(S=64, N=1000).rad()
        mete_lgsr_approx(S=64, N=1000).cdf((1,1,2,4,5,7,12))
    
    
    def test_sugihara_rank_abund(self):
        #Testing against A. J. Baczkowski 1997 paper values
        #Passes with 2% error regularly. Being conservative with 5% error
        ps = sugihara(S=10, N=400).rad(sample_size=20000) / 400
        error = .05 * ps
        diff = np.array(self.sugi) - ps
        ind = np.abs(diff) <= error
        self.assertTrue(np.all(ind))
    

    def test_geo_series(self):
        #Using data from Magurran (1998)
        obs_sad = [370,210,120,66,35,31,15,9,3,2,1]
        mag_pred_sad = [387.6, 213.8, 117.8, 64.5, 35.5, 19.8, 10.7, 6.2, 3.3,\
                        1.8, 1.0]
        gk = geo_series().fit(obs_sad)['k']
        self.assertTrue(np.round(gk, decimals=3) == 0.449)
        geo_sad = np.round(geo_series(S=len(obs_sad), N=sum(obs_sad), k=.449).\
                           rad(), decimals=1)
        diff = np.floor(np.array(mag_pred_sad)) - np.floor(geo_sad)
        self.assertTrue(np.all(diff == 0))

    def test_broken_stick(self):
        #Using data from Magurran (1998)
        self.assertRaises(AssertionError, broken_stick(S=12, N=111).pmf, 112)
        obs_sad = [103,115,13,2,67,36,51,8,6,61,10,21,7,65,4,49,92,37,16,6,23,\
                   9,2,6,5,4,1,3,1,9,2]
        expt = [1.077, 1.040, 1.004, .970, .937, .904, .873, .843, .814,\
                    .786, .759, .732, .707, .683, .659, .636]
        S = len(obs_sad)
        N = sum(obs_sad)
        bs = np.round(broken_stick(S=S, N=N).pmf(np.arange(1, 17)) * S,\
                      decimals=3)
        diff = np.array(expt) - bs
        self.assertTrue(np.all(diff == 0))
        broken_stick(S=23, N=500).cdf([1,2,500])
        broken_stick(S=23, N=500).rad()

    def test_lognorm_pmf(self):
        r_output = [0.1210, .0806, .0601, 0.0476, 0.0391, .0331,  0.0285,\
                    0.0249, 0.0221, 0.0197]
        lnorm = np.round(lognorm(mu=2, sigma=2).pmf(np.arange(1,11)), \
                                                                    decimals=4)
        diff = r_output - lnorm
        self.assertTrue(np.all(diff == 0))
        lognorm().fit(self.abund_list[0])
        N=sum(self.abund_list[0]); S=len(self.abund_list[0])
        lognorm(N=N, S=S, mu=2, sigma=5).rad()
        lognorm(mu=1.5, sigma=3.45).cdf([1,1,4,5,12])

if __name__ == '__main__':
    unittest.main()
