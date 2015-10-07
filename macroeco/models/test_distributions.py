from __future__ import division

from numpy.testing import (TestCase, assert_equal, assert_array_equal,
                           assert_almost_equal, assert_array_almost_equal,
                           assert_allclose, assert_, assert_raises)

import numpy as np
from decimal import Decimal
from macroeco.models import *
from macroeco.models._distributions import _trunc_logser_solver
import matplotlib.pyplot as plt
import scipy as sp
import scipy.stats as stats


class TestGeom(TestCase):

    def test_pmf(self):
        vals = geom.pmf([0,1,2], 0.25)
        assert_array_almost_equal(vals, np.array([0.25, 0.1875, 0.140625]))

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
                                  np.array([0.25,0.1875,0.140625]) / 0.578125)

    def test_cdf(self):
        # Expected values are regular geom cdf divided by cdf at b
        vals = geom_uptrunc.cdf([0,1,2], 0.5, 2)
        assert_array_almost_equal(vals, np.array([0.5,0.75,0.875]) / 0.875)

    def test_cdf_x_len_1(self):
        # cdf should be not throw error even if x is len 1
        vals = geom_uptrunc.cdf(0, 0.5, 2)
        assert_almost_equal(vals, 0.5 / 0.875)

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
        mu, k = nbinom.fit_mle(x, k_array=np.arange(0.01,10,0.01))
        assert_array_almost_equal([mu, k], [9, 8.54], decimal=2)

    def test_alternative_rvs(self):
        rand_alt = nbinom.rvs_alt(5, 1, l=0, size=10000)
        rand = nbinom.rvs(5, 1, size=10000)

        alt_k = nbinom.fit_mle(rand_alt, k_array=np.arange(0.5, 1.5, 0.01))
        k = nbinom.fit_mle(rand, k_array=np.arange(0.5, 1.5, 0.01))

        assert_almost_equal(alt_k, k, decimal=1)


class TestNbinom_ztrunc(TestCase):

    def test_pmf(self):
        # Test pmf gives back expected mean
        tpmf = nbinom_ztrunc.pmf(np.arange(1, 500), 4, 1)
        tmean = np.sum(np.arange(1, 500) * tpmf)
        assert_almost_equal(tmean, 4)

        # Test pmf of 0 is 0
        tpmf = nbinom_ztrunc.pmf(0, 1, 1)
        assert_equal(tpmf, 0)

    def test_cdf(self):

        # Test cdf and pmf agree!
        tpmf = np.sum(nbinom_ztrunc.pmf(np.arange(1, 20), 20, 10))
        tcdf = nbinom_ztrunc.cdf(19, 20, 10)
        assert_almost_equal(tpmf, tcdf)

    def test_get_p_from_mu(self):

        # Test the fit p values are equal to those given in He and Legendre
        # 2002
        test_values = [205.9878, 410.9853, 794.7613, 1210.0497,
                            1945.9970, 3193.8362]
        test_ks = [2, 1, 0.5, 0.3, 0.1363, 0.01]

        ps = np.array([nbinom_ztrunc.translate_args(335356 / 814., tk,
            return_p=True)[0] for tk in test_ks])

        assert_array_almost_equal(ps, test_values, decimal=0)

    def test_fit_mle(self):

        # Test fit returns something close the input
        rvs_data = nbinom_ztrunc(10, 1).rvs(size=1000)
        ml_mean, ml_k = nbinom_ztrunc.fit_mle(rvs_data)
        assert_almost_equal(ml_mean, np.mean(rvs_data))
        assert_almost_equal(ml_k, 1, decimal=0)

        rvs_data = nbinom_ztrunc(20, 10).rvs(size=1000)
        ml_mean, ml_k = nbinom_ztrunc.fit_mle(rvs_data)
        assert_almost_equal(ml_mean, np.mean(rvs_data))
        assert_almost_equal(ml_k, 10, decimal=0)


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
        fig, axes = plt.subplots(3, 3, sharex=True, figsize=(15, 7))
        axes = axes.flatten()

        for i, ax in enumerate(axes):
            ax.plot(np.arange(1, 101), fnbd_vec[i])
            ax.plot(np.arange(1, 101), nbd_vec[i], '--')
            ax.plot(np.arange(1, 101), binm_vec[i], '.-')
            ax.legend(('fnbd', 'nbd', 'binm'), loc='best')
            ax.set_xlabel('abundance')
            ax.set_ylabel('P(x)')
            ax.text(0.6, 0.3, descrip[i], transform=ax.transAxes)

        #Uncomment to save figure
        #fig.savefig("test_cnbinom")


class TestDgamma(TestCase):

    def test_pmf(self):
        # import macroeco_distribution as mac
        # mac.dis_gamma_ll([1,1,2,5,6,7], 5, .3)
        test_val = -32.3085384957
        pred_val = np.sum(dgamma.logpmf([1, 1, 2, 5, 6, 7], 5, .3))
        assert_almost_equal(test_val, pred_val)

        # ab = [1, 1, 1, 1, 2, 4, 4, 4, 4, 4, 45, 267]
        # mac.dis_gamma_ll(ab, 0.1, 200)
        test_val = -39.889246913391531
        ab = [1, 1, 1, 1, 2, 4, 4, 4, 4, 4, 45, 267]
        pred_val = np.sum(dgamma.logpmf(ab, 0.1, 200))
        assert_almost_equal(test_val, pred_val)

    def test_cdf(self):
        # Test that cdf gets close to one
        assert_almost_equal(dgamma.cdf(1000, 4, .9), 1)

    def test_fit_mle(self):
        # mac.dis_gamma_solver([1,1,2,5,6,7])
        fit_alpha = 1.1324749
        fit_theta = 2.86753
        alpha, theta = dgamma.fit_mle([1, 1, 2, 5, 6, 7])
        assert_almost_equal(fit_alpha, alpha, decimal=3)
        assert_almost_equal(fit_theta, theta, decimal=3)

    def test_rank(self):
        # When alpha is almost zero should be similar to logseries with p =
        # e^(-1 / theta)
        logseries_rank = logser_uptrunc.rank(10, np.exp(-1 / 3), 1000)
        dgamma_rank = dgamma.rank(10, 0.0001, 3)

        assert_array_equal(logseries_rank, dgamma_rank)

class TestLogser(TestCase):

    def test_pmf(self):

        # Testing against values in Williams 1944,
        # Some applications of the logarithmic series and the index of
        # diversity to ecological problems, pg. 18.

        # Acridiidae: S = 826, p = 0.92964 (There seems to be an error in
        # their data at 3 -> should be 83.3 not 88.3)
        test_vals = np.array([289.3, 134.5, 83.3, 58.1, 43.2, 33.5, 26.7, 21.7,
         17.9, 15., 12.7, 10.8, 9.3, 8., 6.9, 6.1, 5.3, 4.6, 4.1, 3.6])

        pred_pmf = logser.pmf(np.arange(1, 21), 0.92964)
        pred_vals = np.round(pred_pmf * 826, decimals=1)
        assert_array_equal(test_vals, pred_vals)

        # Mantidae: S = 209, p = 0.89781
        test_vals = np.array([82.3, 36.9, 22.1, 14.9, 10.7, 8., 6.2, 4.8, 3.9,
         3.1, 2.5, 2.1, 1.7, 1.4, 1.2, 1., 0.9, 0.7, 0.6, 0.5])

        pred_pmf = logser.pmf(np.arange(1, 21), 0.89781)
        pred_vals = np.round(pred_pmf * 209, decimals=1)
        assert_array_equal(test_vals, pred_vals)

        # Blattidae: S = 197, p = 0.96476
        test_vals = np.array([56.8, 27.4, 17.6, 12.8, 9.8, 7.9, 6.5, 5.5, 4.7,
         4.1, 3.6, 3.2, 2.8, 2.5, 2.3, 2.1, 1.9, 1.7,
         1.6, 1.4, 1.3, 1.2, 1.1, 1., 1., 0.9, 0.8,
         0.8, 0.7, 0.7])

        pred_pmf = logser.pmf(np.arange(1, 31), 0.96476)
        pred_vals = np.round(pred_pmf * 197, decimals=1)
        assert_array_equal(test_vals, pred_vals)

    def test_translate_args(self):

        # Using values from Williams 1994
        test_vals = [0.92964, 0.89781, 0.96476, 0.97003]
        data = [4112 / 826., 805. / 209, 1612. / 197, 480. / 52]

        pred_vals = [logser.translate_args(td) for td in data]

        assert_array_almost_equal(test_vals, pred_vals, decimal=5)

    def test_fit_mle(self):

        test_val = .97003  # Value from Williams 1944
        x = np.arange(1, 53.)
        norm_x = x / sum(x)
        data = norm_x * (480)
        pred_val = logser.fit_mle(data)
        assert_almost_equal(test_val, pred_val, decimal=5)


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

    def test_rank(self):
        # Test rank against values generated by hand
        exp_vals = np.array([1., 1., 2., 3., 4., 7., 11., 18., 31., 62.])

        # Test values generated
        test_vals = logser_uptrunc.rank(10, .99, 100)

        assert_array_equal(exp_vals, test_vals)

    def test_rvs(self):

        # Make sure random number generator is returning what is expected
        res1 = logser_uptrunc.rvs(.9, 100)
        assert_equal(1, len(np.atleast_1d(res1)))

        res2 = lognorm.rvs(.9, 100, size=5)  # Should be length 5
        assert_equal(5, len(res2))


class TestLognorm(TestCase):

    def test_pmf(self):
        # R pmf: dlnorm(c(1:10), 2, 2)
        r_output = np.array([0.1210, .0806, .0601, 0.0476, 0.0391, .0331,
            0.0285, 0.0249, 0.0221, 0.0197])

        test1 = lognorm.pdf(np.arange(1, 11), 2, 2)
        assert_array_almost_equal(test1, r_output, decimal=4)

        # R pmf: dlnorm(5, -3, 5)
        r_ans = 0.0104333
        test2 = lognorm.pdf(5, -3, 5)
        assert_almost_equal(test2, r_ans)

    def test_cdf(self):
        # R cdf: plnorm(c(1,1,4,5,12), 1.2, 3.45)
        r_output = np.array([0.3639854, 0.3639854, 0.5215318, 0.5472346,
                                        0.6452161])

        test = lognorm.cdf([1, 1, 4, 5, 12], 1.2, 3.45)
        assert_array_almost_equal(test, r_output, decimal=7)

    def test_translate_args(self):

        mean = 67; sigma = 2
        mu, sigma = lognorm.translate_args(mean, sigma)

        # Expected mu: np.log(mean) - (sigma**2 / 2)
        exp_mu = 2.2046926
        assert_almost_equal(mu, exp_mu)

    def test_fit_mle(self):
        '''
        # R code
        pmf <- function(x, N, S, sigma){
            mu = log(N / S) - (sigma^2 / 2)
            dlnorm(x, meanlog=mu, sdlog=sigma)
        }

        mle <- function(sdlog, x, N, S){
            -sum(log(pmf(x, N, S, sdlog)))
        }

        params <- function(x){
            N = sum(x);
            S = length(x);
        optimize(mle, interval=c(0,5), x, N, S)
        }

        data = # some data
        params(data)'''

        data1 = [1, 1, 1, 1, 1, 2, 2, 3, 3, 4, 5, 6, 123, 456]
        data2 = [2, 2, 2, 4, 67, 34, 152, 9]

        r_fits = [2.07598, 1.59213]  # data1, data2

        testfit1 = lognorm.fit_mle(data1, fix_mean=True)[1]
        testfit2 = lognorm.fit_mle(data2, fix_mean=True)[1]

        assert_almost_equal(r_fits[0], testfit1, decimal=5)
        assert_almost_equal(r_fits[1], testfit2, decimal=5)

        # Scipy code: stats.lognorm.fit(data1, floc=0)
        scipy_ans = 1.79518287
        test1 = lognorm.fit_mle(data1)[1]
        assert_almost_equal(scipy_ans, test1)

    def test_rvs(self):

        # Test that multiple random numbers can be returned without error
        res1 = lognorm.rvs(5, 5)  # Should be length 1
        assert_equal(1, len(np.atleast_1d(res1)))

        res2 = lognorm.rvs(5, 5, size=5)  # Should be length 5
        assert_equal(5, len(res2))

    def test_stats(self):

        # Test that stats returns the correct stats
        mu, sigma = lognorm.translate_args(50, 2)
        mean, sigma = lognorm.stats(mu, sigma, moments="mv")
        assert_almost_equal(50, mean)

        res = lognorm.stats(mu, sigma, moments="mvsk")
        assert_equal(len(res), 4)


class TestPlnorm(TestCase):

    def test_pmf(self):

        # Test against R VGAM fxn: dpolono(c(1:10), -1, 3)
        r_res = [0.121392844, 0.057692006, 0.035586652, 0.024863530,
            0.018681089, 0.014721035, 0.011998072, 0.010027588, 0.008545518,
            0.007396607]

        test = plnorm.pmf(np.arange(1, 11), -1, 3)
        assert_array_almost_equal(r_res, test)

        # Test against macroeco_distributions.pln:
        # pln.pmf([0, 50, 1000], 2.34, 5, 0)

        md_res = np.array([2.86468926e-01, 1.51922299e-03, 5.25717609e-05])
        test = plnorm.pmf([0, 50, 1000], 2.34, 5)
        assert_array_almost_equal(md_res, test)

        # Unit test from test_macroeco_distributions

        # Test values for Poisson lognomal are chosen from Table 1 and Table 2
        # in Grundy Biometrika 38:427-434.
        # In Table 1 the values are deducted from 1 which give p(0).
        pln_table1 = [[-2.0, 2, '0.9749'],
                      [-2.0, 8, '0.9022'],
                      [-2.0, 16, '0.8317'],
                      [0.5, 2, '0.1792'],
                      [0.5, 8, '0.2908'],
                      [0.5, 16, '0.3416'],
                      [3, 2, '0.0000'],
                      [3, 8, '0.0069'],
                      [3, 16, '0.0365']]

        pln_table2 = [[-2.0, 2, '0.0234'],
                      [-2.0, 8, '0.0538'],
                      [-2.0, 16, '0.0593'],
                      [0.5, 2, '0.1512'],
                      [0.5, 8, '0.1123'],
                      [0.5, 16, '0.0879'],
                      [3, 2, '0.0000'],
                      [3, 8, '0.0065'],
                      [3, 16, '0.0193']]

        for vals in pln_table1:
            test = plnorm.pmf(0, np.log(10 ** vals[0]), vals[1] ** .5)
            assert_almost_equal(test, float(vals[2]), decimal=4)

        for vals in pln_table2:
            test = plnorm.pmf(1, np.log(10 ** vals[0]), vals[1] ** .5)
            assert_almost_equal(test, float(vals[2]), decimal=4)


    def test_cdf(self):

        # Test against R VGAM fxn: ppolono(c(0, 15, 10000), .1, 2)
        r_res = [0.3954088, 0.9048902, 0.9999973]
        test = plnorm.cdf([0, 15, 10000], .1, 2)
        assert_array_almost_equal(r_res, test, decimal=5)

        # Test against macroeco_distributions:
        # pln.cdf([1,2,3], 20, 4, 0)

        md_res = np.array([7.34761277e-07, 1.18860746e-06, 1.67083480e-06])
        test = plnorm.cdf([1, 2, 3], 20, 4)
        assert_array_almost_equal(md_res, test, decimal=5)

    def test_fit_mle(self):

        # Test against R poilog: poilogMLE(data, zTrune=FALSE)
        data = np.array([1,1,1,1,1,2,2,2,3,3,4,4,5,5,6,6,12,45,67])
        Rfits = (1.31928, 1.18775)
        fits = plnorm.fit_mle(data)
        assert_array_almost_equal(Rfits, fits, decimal=3)

        # Test against macroeco_distributions
        # pln_solver(data, lower_trunc=False)
        md_res = (1.3195580310886075, 1.1876019842774048)
        assert_array_almost_equal(md_res, fits, decimal=4)

    def test_rank(self):
        # This should be a slow test!

        # Test against ppf.
        # >>> n = 50
        # >>> vals = (np.arange(1, n+1) - 0.5) / n
        # >>> plnorm.ppf(vals, 1, 1)
        test_case = np.array([  0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,
            1.,   1.,   1., 1.,   1.,   1.,   1.,   1.,   1.,   2.,   2.,   2.,
            2.,   2., 2.,   2.,   3.,   3.,   3.,   3.,   3.,   3.,   4.,   4.,
            4., 4.,   5.,   5.,   5.,   6.,   6.,   6.,   7.,   7.,   8.,   9.,
            10.,  11.,  13.,  15.,  19.,  29.])

        pred_res = plnorm.rank(50, 1, 1, crit=0.5, upper=40)

        # Test the values are within one
        diff = np.abs(pred_res - test_case)
        zeros = np.sum(diff == 0)
        ones = np.sum(diff == 1)
        assert_equal(zeros + ones, len(diff))


class TestPlnormZtrunc(TestCase):

    def test_pmf(self):

        # Test against macroeco_distributions:
        # pln.pmf([0, 50, 1000], 2.34, 5, 1)
        md_res = np.array([0, 2.12916164e-03, 7.36783061e-05])
        test = plnorm_ztrunc.pmf([0, 50, 1000], 2.34, 5)

        assert_array_almost_equal(md_res, test)

    def test_cdf(self):

        # Test against dpolonorm
        # ppolono(c(1,2,3), 4.3, 100) / (1 - ppolono(0, 4.3, 100))
        r_res = [0.007670365, 0.011507417, 0.014065948]

        test = plnorm_ztrunc.cdf(np.arange(1, 4), 4.3, 100)
        assert_array_almost_equal(r_res, test)

    def test_fit_mle(self):

        data = np.array([1,1,1,4,4,4,4,5,5,5,12,44,55,112])

        # macroeco_distributions fit: pln_solver(data)
        md_fits = (1.068510556981163, 1.8800439687956865)
        test = plnorm_ztrunc.fit_mle(data)
        assert_array_almost_equal(test, md_fits, decimal=4)

        # R poilog: poilogMLE(data)
        r_fits = (1.067620, 1.880646)
        assert_array_almost_equal(test, r_fits, decimal=3)

    def test_rank(self):

        # TODO: Can't test this against ppf because ppf is too slow

        # Make sure it is working when crit = 0
        test = [ 1., 1., 2., 2., 2., 2., 2., 3., 3.,
                4., 5., 5., 6., 6., 7., 7., 8., 11., 14., 22.]
        rad = plnorm_ztrunc.rank(20, 1, 1, crit=0, upper=40)
        assert_array_equal(test, rad)

class TestExpon(TestCase):

    def test_pdf(self):
        vals = expon.pdf([0,1,2], 2.5)
        assert_almost_equal(vals, [2.5, 0.205212497, 0.016844867])

    def test_mean(self):
        mu1 = expon.mean(0.5)
        assert_almost_equal(mu1, 2)

        mu2 = expon.mean(0.25)
        assert_almost_equal(mu2, 4)

    def test_cdf(self):
        vals = expon.cdf([0,1,2], 0.5)
        assert_array_almost_equal(vals, [0, 0.39346934, 0.632120559])

    def test_translate_args(self):
        assert_almost_equal(1/13, expon.translate_args(13))

    def test_fit_mle(self):
        assert_almost_equal(1/8, expon.fit_mle([6,7,9,10]))


class TestExponUptrunc(TestCase):

    def test_pdf(self):
        vals = expon_uptrunc.pdf([0,1,2], 0.2, 10)
        assert_almost_equal(vals, [0.231303529, 0.189375312, 0.155047392])

    def test_pdf_lambda_equal_zero_is_uniform(self):
        vals = expon_uptrunc.pdf([0,1,2], 0.0000001, 10)
        assert_almost_equal(vals, [0.1, 0.1, 0.1])

    def test_pdf_integrates_to_one(self):
        val1 = sp.integrate.quad(expon_uptrunc.pdf, 0, 10, (0.2, 10))
        assert_almost_equal(val1[0], 1)

        val2 = sp.integrate.quad(expon_uptrunc.pdf, 0, 100, (.000000001, 100))
        assert_almost_equal(val2[0], 1)

        val3 = sp.integrate.quad(expon_uptrunc.pdf, 0, 100, (-5, 100))
        assert_almost_equal(val3[0], 1)

    def test_mean_lambda_equal_zero(self):
        # If lam zero (uniform distribution), mean should be 1/2 b
        assert_almost_equal(expon_uptrunc.mean(0.0000001, 10), 5, 5)

    def test_mean(self):
        def integrand(x, lam, b):
            return x * expon_uptrunc.pdf(x, lam, b)

        for lam in [2, 4.5]:
            val = sp.integrate.quad(integrand, 0, 5, args=(lam, 10))[0]
            assert_almost_equal(expon_uptrunc.mean(lam, 5), val, 4)

    def test_cdf(self):
        vals = expon_uptrunc.cdf([0,1,2], 0.2, 10)
        assert_array_almost_equal(vals, [0, 0.209641082, 0.381280683])

    def test_translate_args_uniform_case(self):
        lam = expon_uptrunc.translate_args(5, 10)
        assert_almost_equal(lam[0], 0)

    def test_translate_args(self):
        # mean -> lambda -> mean comparison
        lam = expon_uptrunc.translate_args(3, 10)
        assert_almost_equal(expon_uptrunc.mean(lam, 10), 3)

    def test_fit_mle_uniform_case(self):
        data = [5,5,5]
        mean = np.mean(data)
        lam = expon_uptrunc.fit_mle(data, 10)[0]
        assert_almost_equal(expon_uptrunc.mean(lam, 10), 5, 4)

    def test_fit_mle(self):
        data = [4,5,7,8]
        mean = np.mean(data)
        lam = expon_uptrunc.fit_mle(data, 10)[0]
        assert_almost_equal(expon_uptrunc.mean(lam, 10), 6)

