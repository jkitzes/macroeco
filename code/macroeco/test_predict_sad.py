#Testing predict_sad.py

import unittest
from predict_sad import *
import numpy as np
from rpy2.robjects.packages import importr
import rpy2.robjects as robjects

class TestPredictSAD(unittest.TestCase):
    '''Test the functions within predict_sad.py'''

    def setUp(self):
        
        self.vgam = importr('VGAM')
        self.abund1 = np.array([1,1,2,3,4,4,7,12])
        

    def test_plognorm_pmf(self):
        log_abund = np.log(self.abund1)
        mean = np.mean(log_abund)
        sd = np.var(log_abund, ddof=1)**0.5
        x = robjects.IntVector(np.arange(1,sum(self.abund1) + 1))
        rdpolono = np.array(self.vgam.dpolono(x, meanlog=mean, sdlog=sd))
        pmf = plognorm_pmf(self.abund1)
        diff = np.round(rdpolono - pmf, decimals=9)
        self.assertTrue(np.sum(diff) == 0)

    def test_mete_lgsr_pmf(self):
        self.assertRaises(AssertionError, mete_lgsr_pmf, 45, 45)
        self.assertRaises(AssertionError, mete_lgsr_pmf, 234, 67)
        self.assertRaises(AssertionError, mete_lgsr_pmf, 34, 0)
        pmf = mete_lgsr_pmf(34, 567)
        self.assertTrue(len(pmf) == 567)
        #TODO: More testing

    def test_mete_lgsr_approx_pmf(self):
        self.assertRaises(AssertionError, mete_lgsr_pmf, 45, 45)
        self.assertRaises(AssertionError, mete_lgsr_pmf, 234, 67)
        self.assertRaises(AssertionError, mete_lgsr_pmf, 34, 0)
        #TODO: More testing


if __name__ == '__main__':
    unittest.main()
