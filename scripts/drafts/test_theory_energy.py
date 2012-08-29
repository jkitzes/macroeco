#!/usr/bin/python

'''
Unit tests for predict_energy.py
'''

from __future__ import division
import unittest
import numpy as np
from predict_energy import *

class TestEnergy(unittest.TestCase):
    '''Test the functions within predict_energy.py'''

    def setUp(self):
        self.x = 2

    #More testing should be done
    def test_energy_theta_pdf(self):
        self.assertRaises(AssertionError, mete_energy_theta_pdf, 56, 34, 34, 12)
        self.assertRaises(AssertionError, mete_energy_theta_pdf, 4, 23, 12, 16)
        self.assertRaises(AssertionError, mete_energy_theta_pdf, 4, 23, 56, 23)
        pdf, lambda2 = mete_energy_theta_pdf(4, 4*4, 4*(4*4), 10, testing=True)
        self.assertTrue(np.round(lambda2, decimals=3) == 0.083)
        S = 16
        N = S * (2**8)
        E = N * (2**10)
        pdf, lambda2 = mete_energy_theta_pdf(S, N, E, N - 1, testing=True)
        self.assertTrue(np.round(lambda2, decimals=8) == 3.82e-6)

    def test_mete_energy_nu_pdf(self):
        self.assertRaises(AssertionError, mete_energy_nu_pdf, 23, 12, 45)
        self.assertRaises(AssertionError, mete_energy_nu_pdf, 4, 15, 4)
        self.assertRaises(AssertionError, mete_energy_nu_pdf, 4, 0, 56)
        S = 4
        N = S * (2**4)
        E = N * (2**2)
        pdf, lambda2, beta = mete_energy_nu_pdf(S, N, E, testing=True)
        self.assertTrue(np.round(beta - lambda2, decimals=3) == -0.030)
        self.assertTrue(np.round(lambda2, decimals=3) == 0.021)
        S = 16
        N = S * (2**8)
        E = N * (4)
        pdf, lambda2, beta = mete_energy_nu_pdf(S, N, E, testing=True)
        self.assertTrue(np.round(beta - lambda2, decimals=5) == -0.00089)
        self.assertTrue(np.round(lambda2, decimals=4) == 0.0013)


    def test_mete_energy_psi_pdf(self):
        self.assertRaises(AssertionError, mete_energy_psi_pdf, 23, 12, 45)
        self.assertRaises(AssertionError, mete_energy_psi_pdf, 4, 15, 4)
        self.assertRaises(AssertionError, mete_energy_psi_pdf, 4, 0, 56)
        S = 4
        N = S * (2**4)
        E = N * (2**2)
        pdf, lambda2, beta = mete_energy_psi_pdf(S, N, E, testing=True)
        self.assertTrue(np.round(beta - lambda2, decimals=3) == -0.030)
        self.assertTrue(np.round(lambda2, decimals=3) == 0.021)
        S = 16
        N = S * (2**8)
        E = N * (4)
        pdf, lambda2, beta = mete_energy_psi_pdf(S, N, E, testing=True)
        self.assertTrue(np.round(beta - lambda2, decimals=5) == -0.00089)
        self.assertTrue(np.round(lambda2, decimals=4) == 0.0013)



        


if __name__ == '__main__':
    unittest.main()
