"""
Tests for distributions2 module

"""

from __future__ import division

from numpy.testing import (TestCase, assert_equal, assert_array_equal,
                           assert_almost_equal, assert_array_almost_equal,
                           assert_allclose, assert_, assert_raises)

import numpy as np
from decimal import Decimal
from macroeco.models import *
from macroeco.models._distributions import _trunc_logser_solver
import matplotlib.pyplot as plt
import scipy.stats as stats

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

    def test_fit_mle(self):
        p = geom.fit_mle([1,2,4,5])
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
        # canonical test case for next two tests and for test_fit_mle and
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

    def test_fit_mle(self):
        p1, _ = geom_uptrunc.fit_mle([0,10], 10)
        assert_almost_equal(p1, 0)

        p2, _ = geom_uptrunc.fit_mle([1,3], 16)
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

    def test_fit_mle_with_rvs(self):
        np.random.seed(8)
        x = nbinom.rvs(20, 10, size=100)
        mu, k = nbinom.fit_mle(x)
        assert_array_almost_equal([mu, k], [20, 10], decimal=0)

    def test_fit_mle_with_R(self):
        #> library(MASS)
        #> fitdistr(seq(49), "negative binomial")
        x = np.array(range(1,50))
        mu, k = nbinom.fit_mle(x)
        assert_array_almost_equal([mu, k], [25, 2.4337345], decimal=1)

    def test_fit_mle_with_manual_calc(self):
        x = np.array([6,17,14,12,8,10,4,9,3,12,4,2,12,8,14,16,9,10,8,5,6])
        mu, k = nbinom.fit_mle(x, k_range=(0.01,10,0.01))
        assert_array_almost_equal([mu, k], [9, 8.54], decimal=2)

class TestCnbinom(TestCase):

    def test_pmf(self):
        # Test pmf sums to one
        pmf = cnbinom.pmf(np.arange(0, 101), 20, 1, 100)
        assert_almost_equal(np.sum(pmf), 1)

    def test_cdf(self):
        # Test cdf is one at appropriate value
        cdf = cnbinom.cdf(100, 20, 1, 100)
        assert_almost_equal(cdf, 1)

    def test_fit_of_vector(self):
        # Test fit of vector from Issue #3 (github.com/jkitzes/macroeco)
        data = np.array([3,2,1,0,0,0,0,0,0,0,0,0,0,0,0])
        k_fit = cnbinom.fit_mle(data)[0]
        assert_equal(False, k_fit == -0.26)

    def test_zillio_plots(self):
        """ Test the cnbinom function replicated the Zillio and He plots

        References
        ----------
        Zillio, T and He, F. 2010. Modeling spatial aggregation of finite
        populations. Ecology, 91, 3698-3706

        """

        # Define Preliminary a and k to test
        a = np.array([0.1, .3, .8])
        k = np.array([.1, 1, 10])
        fnbd_vec = []
        nbd_vec = []
        binm_vec = []
        descrip = []

        # Get data
        for ta in a:
            for tk in k:

                fnbd_vec.append(cnbinom.pmf(np.arange(1, 101),
                                                ta * 100, tk, 100))
                nbd_vec.append(nbinom.pmf(np.arange(1, 101), ta * 100, tk))
                binm_vec.append(stats.binom.pmf(np.arange(1, 101), 100, ta))

                descrip.append("a=%s, k=%s" % (ta, tk))

        # Loop through the data and plot it
        fig, axes = plt.subplots(3, 3, sharex=True)
        axes = axes.flatten()

        for i, ax in enumerate(axes):
            ax.plot(np.arange(1, 101), fnbd_vec[i])
            ax.plot(np.arange(1, 101), nbd_vec[i], '--')
            ax.plot(np.arange(1, 101), binm_vec[i], '.-')
            ax.legend(('fnbd', 'nbd', 'binm'), loc='best')
            ax.set_xlabel('abundance')
            ax.set_ylabel('P(x)')
            ax.text(0.6, 0.3, descrip[i], transform=ax.transAxes)

        # plt.tight_layout()
        # Uncomment to see save figure
        # fig.savefig("test_cbinom")


class TestLogserUptrunc(TestCase):

    def test_pmf(self):
        # import macroeco_distributions as mac
        # mac.trunc_logser(.8, 100).pmf(4)
        test_val = logser_uptrunc(.8, 100).pmf(4)
        assert_almost_equal(test_val, 0.063624697299)

        # import macroeco_distributions as mac
        # mac.trunc_logser(.45, 3).pmf(3)
        test_val = logser_uptrunc(.45, 3).pmf(3)
        assert_almost_equal(test_val, 0.052224371373307543)

    def test_cdf(self):
        # import macroeco_distributions as mac
        # mac.trunc_logser(.8, 100).cdf(4)
        test_val = logser_uptrunc(.8, 100).cdf(4)
        assert_almost_equal(test_val, 0.86556098617469057)

        # import macroeco_distributions as mac
        # mac.trunc_logser(.45, 3).cdf(2)
        test_val = logser_uptrunc(.45, 3).cdf(2)
        assert_array_almost_equal(test_val, 0.9477756286266924)

    def test_mean(self):
        # Expected mean is N / S

        N = 500
        S = 30.
        p = logser_uptrunc.translate_args(N / S, N)[0]
        mean = logser_uptrunc.stats(p, N)[0]
        assert_almost_equal(mean, N / S, decimal=5)

    def test_fit_mle(self):
        # Should return same result as translate args
        data = np.arange(1, 40)
        N = np.sum(data)
        S = len(data)

        fits = logser_uptrunc.fit_mle(data)
        assert_array_almost_equal(fits,
                            logser_uptrunc.translate_args(N / S, N),
                            decimal=5)

    def test_translate_args(self):
        # Test that values equal values from John's book (Harte 2011)

        lg = logser_uptrunc.translate_args(4 * 4 / 4, 4 * 4)[0]
        assert_almost_equal(-np.log(lg), 0.0459, decimal=4)

        lg = logser_uptrunc.translate_args(2 ** 4 * 4 / 4, 2 ** 4 * 4)[0]
        assert_almost_equal(-np.log(lg), -0.00884, decimal=5)

        lg = logser_uptrunc.translate_args(2 ** 8 * 4 / 4, 2 ** 8 * 4)[0]
        assert_almost_equal(-np.log(lg), -0.00161, decimal=5)

        lg = logser_uptrunc.translate_args(2 ** 8 * 16 / 16, 2 ** 8 * 16)[0]
        assert_almost_equal(-np.log(lg), 0.000413, decimal=6)

        lg = logser_uptrunc.translate_args(2 ** 12 * 64 / 64, 2 ** 12 * 64)[0]
        assert_almost_equal(-np.log(lg), 0.0000228, decimal=7)

        lg = logser_uptrunc.translate_args(20 / 20, 20)[0]
        assert_equal(0, 0)


    def test_n_close_to_s(self):
        # Test the solver doesn't fail when N is very close to S

        _trunc_logser_solver(2, 3)
        _trunc_logser_solver(3, 4)
        _trunc_logser_solver(100, 101)

class TestExpon(TestCase):
    pass


class TestExponUptrunc(TestCase):
    pass

