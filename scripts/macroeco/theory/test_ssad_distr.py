#!/usr/bin/python

#Testing ssad distributions

import unittest
from ssad_distr import *
import numpy as np

class Test_SSAD_Distributions(unittest.TestCase):

    def setUp(self):
        self.x = 2
        
    def test_bin(self):
        self.assertTrue(type(binm().pmf(0, 34, .5)) == type(np.array([1])))
        self.assertTrue(type(binm().cdf(0, 34, .50)) == type(np.array([1])))
        self.assertTrue(binm().cdf(34, 34, .5) == 1)
        self.assertRaises(AssertionError, binm().pmf, 35, 34, .5)
        self.assertRaises(AssertionError, binm().cdf, 35, 34, .5)
        self.assertRaises(AssertionError, binm().pmf, 35, 34, 1.1)

    def test_pois(self):
        self.assertTrue(np.round(pois().cdf(34, 34, .5), decimals=1) == 1)
        self.assertRaises(AssertionError, pois().pmf, 35, 34, .5)
        self.assertRaises(AssertionError, pois().cdf, 35, 34, .5)
        self.assertRaises(AssertionError, pois().pmf, 35, 34, 1.1)

    def test_nbd(self):
        self.assertTrue(np.round(nbd().cdf(800, 800, .1, .3), decimals=1) == 1)
        self.assertRaises(AssertionError, nbd().pmf, 35, 34, .5, .3)
        self.assertRaises(AssertionError, nbd().cdf, 35, 34, .5, .3)
        self.assertRaises(AssertionError, nbd().pmf, 35, 34, 1.1, .3)



if __name__ == '__main__':
    unittest.main()
        

