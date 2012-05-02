#!/usr/bin/python


#Testing predict sad
#NOTE: Unit test must be called with cwd macroeco/code/macroeco

import unittest
from predict_sad import *
import empirical
import numpy as np
from rpy2.robjects.packages import importr
import rpy2.robjects as robjects
import os 
import glob

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
            pmf = plognorm_pmf(mean, var, sad, pmf_ret=True)
            diff1 = np.round(rdpolono - pmf, decimals=5)
            self.assertTrue(np.sum(diff1) == 0)

    def test_trun_plognorm_pmf(self):
        for sad in self.abund_list:
            mean = np.mean(np.log(sad))
            var = np.var(np.log(sad), ddof=1)
            trunll = -sum(np.log(trun_plognorm_pmf(mean, var, sad)))
            ew = -1 * float(pln_ll(mean, var**0.5, sad))
            import pdb; pdb.set_trace()
            self.assertTrue(round(trunll, ndigits=5) == round(ew, ndigits=5))

            

    def test_mete_lgsr_pmf(self):
        self.assertRaises(AssertionError, mete_lgsr_pmf, 45, 45)
        self.assertRaises(AssertionError, mete_lgsr_pmf, 234, 67)
        self.assertRaises(AssertionError, mete_lgsr_pmf, 34, 0)
        pmf = mete_lgsr_pmf(34, 567, pmf_ret=True)
        self.assertTrue(len(pmf) == 567)
        #TODO: More testing

    def test_mete_lgsr_approx_pmf(self):
        self.assertRaises(AssertionError, mete_lgsr_pmf, 45, 45)
        self.assertRaises(AssertionError, mete_lgsr_pmf, 234, 67)
        self.assertRaises(AssertionError, mete_lgsr_pmf, 34, 0)
        #TODO: More testing


if __name__ == '__main__':
    unittest.main()
