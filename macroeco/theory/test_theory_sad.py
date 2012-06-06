#!/usr/bin/python


#Testing predict sad
#NOTE: Unit test must be called with cwd macroeco/code/macroeco

import unittest
from theory_sad import *
import numpy as np
from rpy2.robjects.packages import importr
import rpy2.robjects as robjects
import os 

gcwd = os.getcwd
pd = os.path.dirname
jp = os.path.join



class TestPredictSAD(unittest.TestCase):
    '''Test the functions within predict_sad.py'''

    def setUp(self):
        
        self.vgam = importr('VGAM')
        self.abund_list = []
        self.abund_list.append(np.array([1,1,2,3,4,4,7,12]))
        self.abund_list.append(np.array([1,1,1,1,1,2]))
        self.abund_list.append(np.array([3,2,2,1,5,6,4,5,7]))
        self.abund_list.append(np.array([1,1,4,5,7,89,100]))

    def test_plognorm_pmf(self):
        for sad in self.abund_list:
            log_abund = np.log(sad)
            mean = np.mean(log_abund)
            var = np.var(log_abund, ddof=1)
            sd = var**0.5
            x = robjects.IntVector(np.arange(1,sum(sad) + 1))
            rdpolono = np.array(self.vgam.dpolono(x, meanlog=mean, sdlog=sd))
            pmf = plognorm_pmf(np.arange(1, sum(sad) + 1), mean, var)
            diff1 = np.round(rdpolono - pmf, decimals=5)
            self.assertTrue(np.sum(diff1) == 0)
        self.assertTrue(sum(np.round(plognorm_pmf([1,2,3,4,5], -3,-3),\
                            decimals=3)) == 0) 

    #Using Ethan White's functions to test
    def test_trun_plognorm_pmf(self):
        for sad in self.abund_list:
            log_abund = np.log(sad)
            mean = np.mean(log_abund)
            var = np.var(log_abund, ddof=1)
            sd = var**0.5
            ew = pln_ll(mean, sd, sad)[0]
            tplog = sum(np.log(trun_plognorm_pmf(sad, mean, var)))
            self.assertTrue(np.round(ew, decimals=2) == np.round(tplog, decimals=2))
    
    def test_mete_lgsr_pmf(self):
        self.assertRaises(AssertionError, mete_lgsr_pmf, 1, 45, 45)
        self.assertRaises(AssertionError, mete_lgsr_pmf, 1, 234, 67)
        self.assertRaises(AssertionError, mete_lgsr_pmf, 1, 34, 0)
        pmf = mete_lgsr_pmf(np.arange(1, 568), 34, 567)
        self.assertTrue(len(pmf) == 567)
        #Testing that values equal values from John's book (Harte 2011)
        pmf, x = mete_lgsr_pmf(1, 4, 4 * 4, testing=True)
        self.assertTrue(np.round(-np.log(x), decimals=4) == 0.0459)
        pmf, x = mete_lgsr_pmf(1, 4, 2**4 * 4,testing=True)
        self.assertTrue(np.round(-np.log(x), decimals=5) == -0.00884)
        pmf, x = mete_lgsr_pmf(1, 4, 2**8 * 4, testing=True)
        self.assertTrue(np.round(-np.log(x), decimals=5) == -0.00161)
        pmf, x = mete_lgsr_pmf(1, 16, 2**8 * 16, testing=True)
        self.assertTrue(np.round(-np.log(x), decimals=6) == 0.000413)
        pmf, x = mete_lgsr_pmf(1, 64, 2**12 * 64, testing=True)
        self.assertTrue(np.round(-np.log(x), decimals=7) == 0.0000228)
        


    def test_mete_lgsr_approx_pmf(self):
        self.assertRaises(AssertionError, mete_lgsr_approx_pmf, 1, 45, 45)
        self.assertRaises(AssertionError, mete_lgsr_approx_pmf, 1, 234, 67)
        self.assertRaises(AssertionError, mete_lgsr_approx_pmf, 1, 34, 0)
        #Testing that values = values from John's book (Harte 2011)
        pmf, x = mete_lgsr_approx_pmf(1, 4, 4 * 4, testing=True)
        self.assertTrue(np.round(-np.log(x), decimals=3) == 0.116)
        pmf, x = mete_lgsr_approx_pmf(1, 4, 2**4 * 4, testing=True)
        self.assertTrue(np.round(-np.log(x), decimals=4) == 0.0148)
        pmf, x = mete_lgsr_approx_pmf(1, 4, 2**8 * 4, testing=True)
        self.assertTrue(np.round(-np.log(x), decimals=6) == 0.000516)
        pmf, x = mete_lgsr_approx_pmf(1, 16, 2**8 * 16, testing=True)
        self.assertTrue(np.round(-np.log(x), decimals=6) == 0.000516)
        pmf, x = mete_lgsr_approx_pmf(1, 64, 2**12 * 64, testing=True)
        self.assertTrue(np.round(-np.log(x), decimals=7) == 0.0000229 or\
                        np.round(-np.log(x), decimals=7) == 0.0000228)

    
    def test_lgsr_pmf(self):
        self.assertRaises(AssertionError, lgsr_pmf, 1, 45, 45)
        self.assertRaises(AssertionError, lgsr_pmf, 1, 234, 67)
        self.assertRaises(AssertionError, lgsr_pmf, 1, 34, 0)
        #Test against value in Fisher's paper Fisher et al. 1943
        pmf, x = lgsr_pmf(1, 240, 15609, testing=True)
        self.assertTrue(np.round(x, decimals=4) == 0.9974)

    #Kind of a lame test...
    def test_geo_pmf(self):
        geo = geo_pmf(np.arange(1, 101), 23, 100)
        self.assertTrue(np.array_equal(neg_binom_pmf(np.arange(1, 101), 23,\
                        100, 1), geo))
        self.assertTrue(len(geo) == 100)
    '''
    def test_macroeco_pmf(self):
        self.assertRaises(AssertionError, macroeco_pmf, 23, 100, 'geom')
        self.assertRaises(AssertionError, macroeco_pmf, 23, 100, 'neg_binom')
        self.assertTrue(np.array_equal(macroeco_pmf(23, 100, 'mete', sad=self.abund_list[0]),\
                                       mete_lgsr_pmf(23, 100, pmf_ret=True)))
        self.assertTrue(np.array_equal(macroeco_pmf(23, 100, 'mete_approx', sad=self.abund_list[0]),\
                                       mete_lgsr_approx_pmf(23, 100, pmf_ret=True)))
        self.assertTrue(np.array_equal(macroeco_pmf(23, 100, 'geo', sad=self.abund_list[0]),\
                                       geo_pmf(23, 100, pmf_ret=True)))

    def test_nll(self):
        self.assertRaises(TypeError, nll, (1,1,1,2,3,4,8))
        self.assertRaises(AssertionError, nll, (1,1,1,2,3,4,8), '')
        self.assertRaises(AssertionError, nll, (1,1,1,2,3,4,8), 'metes')
        self.assertRaises(AssertionError, nll, (1,1,1,2,3,4,8), 'Mete')
        sad = self.abund_list[0]
        negll = nll(tuple(sad), 'mete')
        eqnll = -sum(np.log(mete_lgsr_pmf(len(sad), sum(sad), abundances=tuple(sad))))
        self.assertTrue(negll == eqnll)

    def test_make_rank_abund(self):
        sad_pmf = mete_lgsr_pmf(12, 112, pmf_ret=True)
        rank_abund = make_rank_abund(tuple(sad_pmf), 12)
        self.assertTrue(len(rank_abund) == 12)
    '''



if __name__ == '__main__':
    unittest.main()
