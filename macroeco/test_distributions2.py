"""
Tests for distributions2 module

"""

from __future__ import division

from numpy.testing import (TestCase, assert_equal, assert_array_equal, 
                           assert_almost_equal, assert_array_almost_equal,
                           assert_allclose, assert_, assert_raises)

import numpy as np
from decimal import Decimal
import macroeco.distributions2 as dist2
from macroeco.distributions2 import *
import matplotlib.pyplot as plt

class TestGeom(TestCase):

    def test_pmf(self):
        vals = geom.pmf([0,1,2], 0.5)
        assert_array_almost_equal(vals, [0.5,0.25,0.125])

    def test_mean(self):
        mu1 = geom.mean(0.5)
        assert_almost_equal(mu1, 1)

        mu2 = geom.mean(0.25)
        assert_almost_equal(mu2, 3)

    def test_cdf(self):
        vals = geom.cdf([0,1,2], 0.5)
        assert_array_almost_equal(vals, [0.5,0.75,0.875])
    
    def test_translate_args(self):
        ps = geom.translate_args([10, 20])
        assert_array_almost_equal(ps, [1/11, 1/21])

    def test_fit2(self):
        p = geom.fit2([1,2,4,5])
        assert_almost_equal(p, 0.25)


class TestGeomUptrunc(TestCase):

    def test_pmf(self):
        # Expected values are regular geo cdf divided by cdf at b
        vals = geom_uptrunc.pmf([0,1,2], 0.25, 2)
        assert_array_almost_equal(vals,
                                  np.array([0.25,0.1875,0.140625])/0.578125)

    def test_cdf(self):
        # Expected values are regular geom cdf divided by cdf at b
        vals = geom_uptrunc.cdf([0,1,2], 0.5, 2)
        assert_array_almost_equal(vals, np.array([0.5,0.75,0.875])/0.875)

    def test_cdf_x_len_1(self):
        # cdf should be not throw error even if x is len 1
        vals = geom_uptrunc.cdf(0, 0.5, 2)
        assert_almost_equal(vals, 0.5/0.875)
    
    def test_mean(self):
        mu1 = geom_uptrunc.mean(0.801, 32)
        assert_almost_equal(mu1, 4, decimal=2)

    def test_translate_args_harte_16(self):
        # TODO: The Harte figures appear to be inaccurate, generate better
        # canonical test case for next two tests and for test_fit2 and
        # test_mean

        # From Harte 2011, Oxford U Press, Tab 7.4, n0=16 row, Eq 7.50
        b = 16
        mu = np.array([2, 1])  # A0/8, A0/16
        expected = np.array([1-0.669, 1-0.500])
        ps, _ = geom_uptrunc.translate_args(mu, b)
        assert_almost_equal(ps, expected, decimal=3)

    def test_translate_args_harte_32(self):
        # From Harte 2011, Oxford U Press, Tab 7.4, n0=32 row, Eq 7.50
        b = 32
        mu = np.array([4, 2])  # A0/8, A0/16
        expected = np.array([1-0.801, 1-0.667])
        ps, _ = geom_uptrunc.translate_args(mu, b)
        assert_almost_equal(ps, expected, decimal=3)

    def test_translate_args_mqwilber_hand_calc(self):
        # TODO: Confirm last 4 of tests, which more accurate
        b = np.array([60, 340, 34])
        mu = np.array([60*.1, 340*.6, 34*.9])
        expected = np.array([1-.8572, 1-1.0036, 1-1.2937])
        ps, _ = geom_uptrunc.translate_args(mu, b)
        assert_almost_equal(ps, expected, decimal=3)

    def test_translate_args_with_sum_of_pmf(self):
        p1, b1 = geom_uptrunc.translate_args(341/4, 341)  # Issue 33
        assert_array_almost_equal(1,np.sum(geom_uptrunc.pmf(range(342),p1,b1)))

        p2, b2 = geom_uptrunc.translate_args(120, 200)  # Arbitrary
        assert_array_almost_equal(1,np.sum(geom_uptrunc.pmf(range(201),p2,b2)))

    def test_fit2(self):
        p1, _ = geom_uptrunc.fit2([0,10], 10)
        assert_almost_equal(p1, 0)

        p2, _ = geom_uptrunc.fit2([1,3], 16)
        assert_almost_equal(p2, 1-0.669, decimal=2)


class TestNbinom(TestCase):

    def test_pmf(self):
        #> dnbinom(c(0,1,2), 3, mu=5)
        #[1] 0.05273438 0.09887695 0.12359619
        vals = nbinom.pmf([0,1,2], 5, 3)
        assert_array_almost_equal(vals, [0.05273438, 0.09887695, 0.12359619])

    def test_cdf(self):
        #> pnbinom(c(0,1,2),2,mu=30)
        #[1] 0.00390625 0.01123047 0.02153015
        vals = nbinom.cdf([0,1,2], 30, 2)
        assert_array_almost_equal(vals, [0.00390625, 0.01123047, 0.02153015])

    def test_mean_var(self):
        mu1, var1 = nbinom.stats(20, 2, moments='mv')
        assert_array_almost_equal([mu1, var1], [20, 20+(20**2)/2])

    def test_get_p_from_mu(self):
        assert_almost_equal(nbinom._get_p_from_mu(10, 2), 2/12)

    def test_fit2_with_rvs(self):
        np.random.seed(8)
        x = nbinom.rvs(20, 10, size=100)
        mu, k = nbinom.fit2(x)
        assert_array_almost_equal([mu, k], [20, 10], decimal=0)

    def test_fit2_with_R(self):
        #> library(MASS)
        #> fitdistr(seq(49), "negative binomial")
        x = np.array(range(1,50))
        mu, k = nbinom.fit2(x)
        assert_array_almost_equal([mu, k], [25, 2.4337345], decimal=1)

    def test_fit2_with_manual_calc(self):
        x = np.array([6,17,14,12,8,10,4,9,3,12,4,2,12,8,14,16,9,10,8,5,6])
        mu, k = nbinom.fit2(x, k_range=(0.01,10,0.01))
        assert_array_almost_equal([mu, k], [9, 8.54], decimal=2)
        

class TestExpon(TestCase):
    pass


class TestExponUptrunc(TestCase):
    pass

