#!/usr/bin/python

#Testing distributions

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

Weecology functions from https://github.com/weecology/macroecotools.
'''

import unittest
from macroeco.distributions import *
import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt

# TODO: Need to add fit functions to tests with new fit functions.

# TODO: Do we need to test rad's? Against what?


class TestDistributions(unittest.TestCase):
    '''Test the functions within sad_distr.py'''

    def setUp(self):
        self.abund_list = []
        self.abund_list.append(np.array([1,1,2,3,4,4,7,12]))
        self.abund_list.append(np.array([1,1,1,1,1,2]))
        self.abund_list.append(np.array([3,2,2,1,5,6,4,5,7]))
        self.abund_list.append(np.array([1,1,4,5,7,89,100]))

        self.sar = ([0.0004, 0.0025, 0.005, 0.01, 0.25, 1], [2.2356, 10.0175,
        15.87, 24.32, 101.25, 155])
        self.sad = np.arange(1, 156)


    def test_logser(self):
        # Test error raising
        self.assertRaises(AssertionError, logser(n_samp=234, tot_obs=67).pmf,
                                                                             1)
        self.assertRaises(AssertionError, logser(n_samp=34, tot_obs=0).pmf, 1)

        # Testing the logser raises error when given data with zeros
        self.assertRaises(ValueError, logser().fit, [[0,1,2,3,4]])

        # Test pmf against value in Fisher's paper Fisher et al. 1943
        lgser = logser(n_samp=240, tot_obs=15609)
        pmf = lgser.pmf(1)
        self.assertTrue(np.round(lgser.var['p'][0], decimals=4) == 0.9974)

        # Test cdf reaches 1
        cdf = np.round(logser(n_samp=45, tot_obs=1200).cdf(1200)[0][0],
                                                                    decimals=1)
        self.assertTrue(cdf == 1)

        # Test known value of cdf when p = 0.90335
        known_cdf = 0.737623 # From scipy.stats.logser
        pred_cdf = logser(n_samp=14, tot_obs=56).cdf(4)[0][0]
        self.assertTrue(np.round(pred_cdf, decimals = 6) == known_cdf)


    def test_logser_ut(self):
        # Test error raising
        self.assertRaises(AssertionError, logser_ut(n_samp=234, tot_obs=67).pmf, 1)
        self.assertRaises(AssertionError, logser_ut(n_samp=34, tot_obs=0).pmf, 1)

        # Test that pmf is correct length
        pmf = logser_ut(n_samp=34, tot_obs=567).pmf(np.arange(1, 568))[0]
        self.assertTrue(len(pmf) == 567)

        # Test that values equal values from John's book (Harte 2011)
        lg = logser_ut(n_samp=4, tot_obs=4 * 4)
        pmf = lg.pmf(1)
        self.assertTrue(np.round(-np.log(lg.var['x'][0]), decimals=4) == 0.0459)
        lg = logser_ut(n_samp=4, tot_obs=2**4 * 4)
        pmf = lg.pmf(1)
        self.assertTrue(np.round(-np.log(lg.var['x'][0]), decimals=5) == -0.00884)
        lg = logser_ut(n_samp=4, tot_obs=2**8 * 4)
        pmf = lg.pmf(1)
        self.assertTrue(np.round(-np.log(lg.var['x'][0]), decimals=5) == -0.00161)
        lg = logser_ut(n_samp=16, tot_obs=2**8 * 16)
        pmf  = lg.pmf(1)
        self.assertTrue(np.round(-np.log(lg.var['x'][0]), decimals=6) == 0.000413)
        lg = logser_ut(n_samp=64, tot_obs=2**12 * 64)
        pmf  = lg.pmf(1)
        self.assertTrue(np.round(-np.log(lg.var['x'][0]), decimals=7) == 0.0000228)

        # Check that they don't fail
        logser_ut(n_samp=64, tot_obs=1000).rad()
        logser_ut(n_samp=64, tot_obs=1000).cdf((1,1,2,4,5,7,12))

        # Test correct answer when n_samp == tot_obs
        lg = logser_ut(n_samp=31, tot_obs=31)
        pmf  = lg.pmf([1,2,3,4,5])
        self.assertTrue(lg.var['x'][0] == 0)
        self.assertTrue(np.array_equal(pmf[0], np.array([1,0,0,0,0])))


    def test_logser_ut_appx(self):
        # Test error raising
        self.assertRaises(AssertionError, logser_ut_appx(n_samp=234,
                                                            tot_obs=67).pmf, 1)
        self.assertRaises(AssertionError, logser_ut_appx(n_samp=34,
                                                             tot_obs=0).pmf, 1)

        # Test that values equal values from John's book (Harte 2011)
        lg = logser_ut_appx(n_samp=4, tot_obs=4 * 4)
        pmf = lg.pmf(1)
        self.assertTrue(np.round(-np.log(lg.var['x'][0]), decimals=3) == 0.116)
        lg = logser_ut_appx(n_samp=4, tot_obs=2**4 * 4)
        pmf  = lg.pmf(1)
        self.assertTrue(np.round(-np.log(lg.var['x'][0]), decimals=4) == 0.0148)
        lg = logser_ut_appx(n_samp=4, tot_obs=2**8 * 4)
        pmf  = lg.pmf(1)
        self.assertTrue(np.round(-np.log(lg.var['x'][0]), decimals=6) == 0.000516)
        lg = logser_ut_appx(n_samp=16, tot_obs=2**8 * 16)
        pmf  = lg.pmf(1)
        self.assertTrue(np.round(-np.log(lg.var['x'][0]), decimals=6) == 0.000516)
        lg = logser_ut_appx(n_samp=64, tot_obs=2**12 * 64)
        pmf  = lg.pmf(1)
        self.assertTrue(np.round(-np.log(lg.var['x'][0]), decimals=7) == 0.0000229
                        or
                        np.round(-np.log(lg.var['x'][0]), decimals=7) == 0.0000228)

        lg = logser_ut_appx(n_samp=31, tot_obs=31)
        pmf = lg.pmf([1,2,3,4,5])
        self.assertTrue(lg.var['x'][0] == 0)
        self.assertTrue(np.array_equal(pmf[0], np.array([1,0,0,0,0])))

        # Test that they don't fail
        logser_ut_appx(n_samp=64, tot_obs=1000).rad()
        logser_ut_appx(n_samp=64, tot_obs=1000).cdf((1,1,2,4,5,7,12))


    def test_plognorm(self):
        # TODO: Should test against Ethans psolver

        # R output from function dpolono in VGAM package
        R = [0.8453224, 1.951546, 1.040038, 0.3102524]

        # Test known values of pmf from R against our results
        for i, sad in enumerate(self.abund_list):
            log_abund = np.log(sad)
            mu = np.mean(log_abund)
            var = np.var(log_abund, ddof=1)
            sigma = var**0.5
            pmf = sum(plognorm(mu=mu, sigma=sigma).pmf(sad)[0])
            diff1 = np.round(R[i], decimals=5) - np.round(pmf, decimals=5)
            self.assertTrue(abs(diff1) == 0)

        # Test pmf is zero when mu or sigma negative
        self.assertTrue(sum(np.round(plognorm(mu=-3,sigma=3).\
                                     pmf([1,2,3,4,5])[0], decimals=3)) == 0)
        self.assertTrue(sum(np.round(plognorm(mu=3,sigma=-3).\
                                     pmf([1,2,3,4,5])[0], decimals=3)) == 0)

        # Test that MLE fit matches R package poilog
        Rmu = 1.31928; Rsigma = 1.18775
        test_vec1 = np.array([1,1,1,1,1,2,2,2,3,3,4,4,5,5,6,6,12,45,67])
        test_plog = plognorm().fit([test_vec1])
        print test_plog.params
        self.assertTrue(np.round(test_plog.params['mu'][0], decimals = 5) ==
                            Rmu)
        self.assertTrue(np.round(test_plog.params['sigma'][0], decimals = 5) ==
                            Rsigma)

        # Test that these don't fail
        plognorm().fit([self.abund_list[0]])
        plognorm(mu=2, sigma=2).cdf(5)


    def test_plognorm_lt(self):

        #Test our pmf against R's poilog
        R_zero_trun = [0.11620, 0.07216, 0.05201, 0.04049, 0.02783, 0.02398,
                       0.00686]
        pred_plog = plognorm_lt(mu=2, sigma=3).pmf([1,2,3,4,6,7,23])[0]
        self.assertTrue(np.array_equal(R_zero_trun, np.round(pred_plog,
                                                                  decimals=5)))

        # Test that cdf sums appropriately
        pred_pmf_plog = plognorm_lt(mu=1.2, sigma=1.5).\
                                                   pmf([1,2,3,4,5,6,7,8])[0]
        pred_cdf_plog = plognorm_lt(mu=1.2, sigma=1.5).cdf(8)[0][0]
        self.assertTrue(np.sum(pred_pmf_plog) == pred_cdf_plog)

        # Test fit against Ethan White results and poilog
        EW_fit = {'mu' :.90, 'sigma' : 2.18}
        R_fit = {'mu' : .904, 'sigma' : 2.184}
        sad = [1,1,1,1,2,3,5,6,12,13,15,23,45,67,112]
        dist = plognorm_lt().fit([sad])
        mu = dist.params['mu'][0]
        sigma = dist.params['sigma'][0]
        self.assertTrue(EW_fit['mu'] == np.round(mu, decimals=2))
        self.assertTrue(EW_fit['sigma'] == np.round(sigma, decimals=2))
        self.assertTrue(R_fit['mu'] == np.round(mu, decimals=3))
        self.assertTrue(R_fit['sigma'] == np.round(sigma, decimals=3))

        # Test that these don't fail
        plognorm_lt(mu=[2,3], sigma=[2,3]).cdf(5)
        plognorm_lt(mu=2, sigma=2).pmf([2,3,4,5,23])
        plognorm_lt().fit([self.abund_list[0]])
        plognorm_lt(mu=10, sigma=1).cdf(45)


    def test_lognorm(self):

        # Test pmf against R output
        r_output = [0.1210, .0806, .0601, 0.0476, 0.0391, .0331,  0.0285,\
                    0.0249, 0.0221, 0.0197]

        r_output2 = [0.1522, 0.1326, 0.1048, 0.0827, 0.0662, 0.0538, 0.0443,
                     0.0198, 0.0012]

        lnorm = np.round(lognorm(tot_obs=np.exp(4), n_samp=1, sigma=2).\
                                    pmf(np.arange(1,11))[0], decimals=4)
        diff = r_output - lnorm
        self.assertTrue(np.all(diff == 0))

        lnorm = np.round(lognorm(tot_obs = np.exp(1.5 + (1.2**2 / 2)) * 50,
                         n_samp=50,sigma=1.2).pmf([1,2,3,4,5,6,7,12,45])[0],
                         decimals=4)
        diff = r_output2 - lnorm
        self.assertTrue(np.all(diff == 0))

        # Test cdf against R cdf
        rcdf = np.array([0.3319, 0.3319, 0.4869, 0.5127, 0.6124])
        pycdf = np.round(lognorm(tot_obs=np.exp(1.5 + (3.45**2 / 2)), n_samp=1,
                            sigma=3.45).cdf([1,1,4,5,12])[0], decimals=4)
        diff = rcdf - pycdf
        self.assertTrue(np.all(diff == 0))

        # Test fit function returns expected result: R code
        '''
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
        }'''
        fit_array1 = [1,1,1,1,1,2,2,3,3,4,5,6,123,456]
        fit_array2 = [2,2,2,4,67,34,152,9]
        r_lognorm_fits = np.array([2.07598, 1.59213])
        pyfit1 = lognorm().fit([fit_array1]).params['sigma'][0]
        pyfit2 = lognorm().fit([fit_array2]).params['sigma'][0]
        diff = r_lognorm_fits - np.round([pyfit1, pyfit2], decimals=5)
        self.assertTrue(np.all(diff == 0))

        # Test that these don't fail
        lognorm().fit([self.abund_list[0]])
        tot_obs=sum(self.abund_list[0])
        n_samp=len(self.abund_list[0])
        lognorm(tot_obs=tot_obs, n_samp=n_samp, mu=2, sigma=5).rad()

        # Test fit parameter length is correct
        dist = lognorm().fit(self.abund_list)
        dist.pmf(3)
        dist.pmf([[3],[4],[5],[6]])
        self.assertTrue(len(dist.params['tot_obs']) == 4)


    def test_geo_ser(self):
        # TODO: Test pmf.
        # Visually, the CDF should be a straight line on a log(abundance) vs.
        # cumulative density plot.  This is true.  However, not sure that my
        # normalization of the continuous geometric series provided by May is
        # ok...

        # Data from Magurran (1998)
        obs_sad = [370,210,120,66,35,31,15,9,3,2,1]
        mag_pred_sad = [387.6, 213.8, 117.8, 64.5, 35.5, 19.8, 10.7, 6.2, 3.3,\
                        1.8, 1.0]

        # Test fit
        dist = geo_ser().fit([obs_sad])
        gk = dist.params['k'][0]
        self.assertTrue(np.round(gk, decimals=3) == 0.449)

        # Test rad
        geo_sad = np.round(geo_ser(n_samp=len(obs_sad), tot_obs=sum(obs_sad),
                                                k=.449).rad(), decimals=1)[0]
        print geo_sad
        diff = np.floor(np.sort(np.array(mag_pred_sad))) - np.floor(geo_sad)
        self.assertTrue(np.all(diff == 0))

        # Test length of var is correct
        dist = geo_ser().fit(self.abund_list)
        self.assertTrue(len(dist.params['k']) == 4)


    def test_broken_stick(self):
        # Test that n_except throws approriate error if length n_samp and tot_obs are not
        # the same as length pmf
        self.assertRaises(TypeError, broken_stick(n_samp=12, tot_obs=111).
                                                                pmf,[[2],[3]])

        # Data from Magurran (1998)
        obs_sad = [103,115,13,2,67,36,51,8,6,61,10,21,7,65,4,49,92,37,16,6,23,\
                   9,2,6,5,4,1,3,1,9,2]
        expt = [1.077, 1.040, 1.004, .970, .937, .904, .873, .843, .814,\
                    .786, .759, .732, .707, .683, .659, .636]

        # Test pmf
        n_samp = len(obs_sad)
        tot_obs = sum(obs_sad)
        bs = np.round(broken_stick(n_samp=n_samp, tot_obs=tot_obs).
                            pmf(np.arange(1, 17))[0] * n_samp, decimals=3)
        diff = np.array(expt) - bs
        self.assertTrue(np.all(diff == 0))

        # Test that these don't fail
        broken_stick(n_samp=23, tot_obs=500).cdf([1,2,500])
        broken_stick(n_samp=23, tot_obs=500).rad()

        # Test basic fit
        check = broken_stick().fit(self.abund_list)
        self.assertTrue(np.all(check.params['tot_obs'] == np.array([sum(ab) for
                        ab in self.abund_list])))
        self.assertTrue(np.all(check.params['n_samp'] == np.array([len(ab) for
                        ab in self.abund_list])))

    def test_dgamma(self):

        # Don't have any good published graphs to test it against. Test
        # everything is working

        obs_sad = [103,115,13,2,67,36,51,8,6,61,10,21,7,65,4,49,92,37,16,6,23,\
                   9,2,6,5,4,1,3,1,9,2]
        dg = dgamma().fit([obs_sad])

        # Check that the parameters are in vars
        self.assertTrue('alpha' in dg.var)
        self.assertTrue('theta' in dg.var)

        # Check that the distribution sums to one.
        pmf = dg.pmf(np.arange(1, sum(obs_sad)))[0]
        self.assertTrue(np.round(sum(pmf), decimals=1) == 1)

    def test_sugihara(self):
        # Data from A. J. Baczkowski 1997
        sugi = [.4761, .1962, .1180, .0751, .0499, .0337, .0226, .0148,\
                     .0090, .0047]

        # Test rad
        # Passes with 2% error regularly. Being conservative with 5% error
        ps = sugihara(n_samp=10, tot_obs=400).rad(sample_size=20000)[0] / 400
        error = .05 * ps
        diff = np.sort(np.array(sugi)) - ps
        ind = np.abs(diff) <= error
        self.assertTrue(np.all(ind))

        # Test that error is raised for cdf, pdf, pmf methods
        self.assertRaises(NotImplementedError, sugihara().pmf, 67)
        self.assertRaises(NotImplementedError, sugihara().cdf, 34)
        self.assertRaises(NotImplementedError, sugihara().pdf, 23)


    def test_binm(self):
        # Using scipy.binom which is already unit tested.

        # Check that pdf and cdf give correct answers
        dist = binm(tot_obs=8123, n_samp=10)
        self.assertTrue(dist.cdf(8123)[0][0] == 1)
        self.assertTrue(dist.cdf(0) == dist.pmf(0))
        self.assertTrue(dist.cdf(1)[0][0] == (sum(dist.pmf([0,1])[0])))

        # Check that appropriate errors are raised
        self.assertRaises(TypeError, dist.pmf, [[1], [1]])
        self.assertRaises(TypeError, binm(g=3, N=45).pmf, (23,24))

        # Check that fit works
        dist = binm().fit(self.abund_list)

    def test_pois(self):
        # Using scipy.poisson which is already unit tested

        # Check that pdf and cdf give correct answers
        dist = pois(tot_obs=112, n_samp=20)
        self.assertTrue(dist.cdf(112)[0][0] == 1)
        a = np.round(dist.cdf(0)[0][0], decimals=12)
        b = np.round(dist.pmf(0)[0][0], decimals=12)
        self.assertTrue(a == b)
        a = np.round(dist.cdf(23)[0], decimals=12)

        # Make sure that fit and rad work
        dist = pois().fit(self.abund_list)
        rads = dist.rad()
        self.assertTrue(len(rads) == 4)

    def test_nbd(self):

        # Test TypeError if k not given
        dist = nbd(tot_obs=2300, n_samp=45)
        self.assertRaises(TypeError, dist.pmf, 45)

        # Test that cdf is about 1 at tot_obs if tot_obs is large
        dist = nbd(tot_obs=2300, n_samp=24, k=2)
        self.assertTrue(np.round(dist.cdf(2300)[0][0], decimals=1) == 1.0)

        # Test the nbd fits a geometric distribution with k=1
        np.random.seed(345)
        p = 0.001
        geo_data = np.random.geometric(p, size=10000)
        dist = nbd().fit([geo_data])
        self.assertTrue(np.round(dist.params['k'][0], decimals=1) == 1)

    def test_nbd_lt(self):
        # TODO: test pmf

        # Test that cdf is about one
        dist = nbd_lt(tot_obs=2300, n_samp=45, k=3)
        d = dist.cdf(2300)[0][0]
        self.assertTrue(np.round(d, decimals=1) == 1.0)

        # Multiple entries both yield cdf with 1
        dist = nbd_lt(tot_obs=[400, 600], n_samp=[30, 23], k=[3,2])
        cdf = dist.cdf([[400], [600]])
        a = np.round(cdf[0][0], decimals=1)
        b = np.round(cdf[0][0], decimals=1)
        self.assertTrue(a == b)

        # Test the fit p values are equal to those given in He and Legendre
        # 2002
        # I am rounding to the nearest whole number, those I have confirmed
        # that the decimals are very close too
        he_values = np.round([205.9878, 410.9853, 794.7613, 1210.0497, 1945.9970,
                                3193.8362], decimals=0)
        he_ks = [2, 1, 0.5, 0.3, 0.1363, 0.01]
        tnbd = nbd_lt(tot_obs=335356, n_samp=814, k=he_ks)
        tnbd.pmf(1)
        pred = np.round(tnbd.var['p'], decimals=0)
        print pred
        print he_values
        self.assertTrue(np.array_equal(he_values, pred))

        # Test that fixing the bias leads to the proper mean
        ks = np.linspace(.01, 5, num=100)
        vals = np.arange(1,1000)
        for k in ks:
            ob = nbd_lt(tot_obs=500, n_samp=20, k=k)
            pred_vals = ob.pmf(vals)[0]
            bmean = sum(vals * pred_vals)
            self.assertTrue(np.round(bmean, decimals=0) == 500 / 20.)

    def test_fnbd(self):

        # Test that no error is thrown if a zero is passed
        fnbd().fit([[0,1,2,3,4,5,6]])

        # TypeError if k is not given
        dist = fnbd(tot_obs=2300, n_samp=20)
        self.assertRaises(TypeError, dist.pmf, 45)

        # Test that cdf sums to one
        dist = fnbd(tot_obs=2300, n_samp=15, k=2)
        self.assertTrue(np.round(dist.cdf(2300)[0][0], decimals=1) == 1.0)

        # Test that cdf and pdf at 0 give the same answer
        a = np.round(dist.pmf(0)[0][0], decimals=12)
        b = np.round(dist.cdf(0)[0][0], decimals=12)
        self.assertTrue(a == b)

        # Test that fit of of geometric data gives k=1
        # NOTE: Testing the fit of FNBD against Zillio and He values is
        # difficult because the abundance data is not given.
        np.random.seed(345)
        p = 0.001
        geo_data = np.random.geometric(p, size=10000)
        dist = fnbd().fit([geo_data])
        self.assertTrue(np.round(dist.params['k'][0], decimals=1) == 1)

        # Test against published data in Zillio and He 2010
        # Generated plots match exactly with plots in Zillio and He, 2010
        # Code to generate plots: Unquote and run nosetest if you want to see
        # the plots

        a = np.array([0.1, .3, .8])
        k = np.array([.1, 1, 10])
        fnbd_vec = []
        nbd_vec = []
        binm_vec = []
        descrip = []
        for ta in a:
            for tk in k:
                fnbd_vec.append(fnbd(tot_obs=100, n_samp = 1./ta,
                                            k=tk).pmf(np.arange(1,101))[0])
                nbd_vec.append(nbd(tot_obs=100, n_samp = 1./ta,
                                              k=tk).pmf(np.arange(1,101))[0])
                binm_vec.append(binm(tot_obs=100, n_samp = 1./ta
                                              ).pmf(np.arange(1,101))[0])

                descrip.append("a=%s, k=%s" % (ta, tk))
        for i in xrange(len(fnbd_vec)):
            plt.clf()
            plt.plot(np.arange(1,101), fnbd_vec[i])
            plt.plot(np.arange(1,101), nbd_vec[i], '--')
            plt.plot(np.arange(1,101), binm_vec[i], '.-')
            plt.legend(('fnbd', 'nbd', 'binm'), loc='best')
            plt.xlabel('abundance')
            plt.ylabel('P(x)')
            plt.text(plt.xlim()[1] * 0.6, plt.ylim()[1] * 0.8, descrip[i])
            #plt.show()
            plt.clf()

        # Based on Zillio and He 2010, Calculating a few pmf values by hand.
        # Going to test the fnbd against these values.

    def test_geo(self):
        # This is just a wrapper function for nbd. Already tested. Will just
        # confirm that everything runs.

        # Test fit work and returns expected results
        test = geo().fit([[0,0,0,1,4,67], [1,1,3,5,23]])
        self.assertTrue(np.all(test.params['tot_obs'] == np.array([72, 33])))
        self.assertTrue(np.all(test.params['n_samp'] == np.array([6,5])))

        # Test that tot_obs is broadcast
        test  = geo(tot_obs=456, n_samp = [34,56,12])
        test.pmf(0)
        self.assertTrue(np.all(test.params['tot_obs'] == np.array(
                                                               [456,456,456])))

        # Test that error is thrown if broadcasting can't be done
        test = geo(tot_obs=[456,320], n_samp = [34,56,12])
        self.assertRaises(ValueError, test.pmf, 0)

        # Test that cdf and pdf return same value as nbd with k = 1
        test_geo = geo(tot_obs=[50, 34], n_samp=2).pmf([0,1,2])
        test_nbd = nbd(tot_obs=[50,34], n_samp=2, k=1).pmf([0,1,2])
        self.assertTrue((len(test_geo) == len(test_nbd)) and (len(test_geo) ==
                                                                            2))
        for i in xrange(len(test_geo)):
            self.assertTrue(np.array_equal(test_geo[i], test_nbd[i]))

    def test_fgeo(self):

        # Test fit work and returns expected results
        test = fgeo().fit([[0,0,0,1,4,67], [1,1,3,5,23]])
        self.assertTrue(np.all(test.params['tot_obs'] == np.array([72, 33])))
        self.assertTrue(np.all(test.params['n_samp'] == np.array([6,5])))

        # Test that tot_obs is broadcast
        test  = fgeo(tot_obs=456, n_samp = [34,56,12])
        test.pmf(0)
        self.assertTrue(np.all(test.params['tot_obs'] == np.array(
                                                               [456,456,456])))

        # Test that error is thrown if broadcasting can't be done
        test = fgeo(tot_obs=[456,320], n_samp = [34,56,12])
        self.assertRaises(ValueError, test.pmf, 0)

        # Test that cdf and pdf return same value as fnbd with k = 1
        test_fgeo = fgeo(tot_obs=[50, 34], n_samp=2).pmf([0,1,2])
        test_fnbd = fnbd(tot_obs=[50,34], n_samp=2, k=1).pmf([0,1,2])
        self.assertTrue((len(test_fgeo) == len(test_fnbd)) and (len(test_fgeo) ==
                                                                            2))
        for i in xrange(len(test_fgeo)):
            self.assertTrue(np.array_equal(test_fgeo[i], test_fnbd[i]))

    def test_tgeo(self):

        # Test against values from Harte 2011
        x_vals = [0.333, 0.434, .568, .707, .823, .901]
        tg = tgeo(tot_obs=[1,2,4,8,16,32], n_samp=4)
        tg.pmf(0)
        pred_vals = np.round(tg.var['x'], 3)
        print pred_vals
        self.assertTrue(np.array_equal(x_vals, pred_vals))

        # In Harte 2011 .143 is given as .125, but this is a mistake.  Every
        # other value is exactly as expected from teh the book
        x_vals = [0.143, .220, .344, .505, .669, .801]
        tg = tgeo(tot_obs=[1,2,4,8,16,32], n_samp=8)
        tg.pmf(0)
        pred_vals = np.round(tg.var['x'], 3)
        print pred_vals
        self.assertTrue(np.array_equal(x_vals, pred_vals))

        x_vals = [0.067, .115, .201, .334, .5, .667]
        tg = tgeo(tot_obs=[1,2,4,8,16,32], n_samp=16)
        tg.pmf(0)
        pred_vals = np.round(tg.var['x'], 3)
        print pred_vals
        self.assertTrue(np.array_equal(x_vals, pred_vals))

        # Test tgeo cdf is one
        dist = tgeo(n_samp=10, tot_obs=2345)
        self.assertTrue(np.round(dist.cdf(2345)[0][0], decimals=1) == 1.0)

        # When n_samp < 2 weird things happen
        # Testing Lagrange multiplier against values generated by hand
        # [(n=60, a=.1), (n=340, a=.6), (n=34, a=.9), (n=12, a=.9), (n=2, .9),
        # (n=1, a=.1),(n=1, a=0.0001),
        x_vals = np.array([.8572, 1.0036, 1.2937, 1.8298, 5.6056, 0.1111])
        tg = tgeo(tot_obs=[60,340,34,12, 2, 1],
                         n_samp=(1./.1, 1/.6, 1/.9, 1/.9, 1/.9, 1/.1))
        tg.pmf(0)
        pred_vals = np.round(tg.var['x'], decimals=4)
        self.assertTrue(np.array_equal(x_vals, pred_vals))

        x_vals = np.array([1.0e-4, 1.0e-5])
        tg = tgeo(tot_obs=[1,1], n_samp=[1/.0001, 1/.00001])
        tg.pmf(0)
        pred_vals = np.round(tg.var['x'], decimals=6)
        self.assertTrue(np.array_equal(x_vals, pred_vals))

        # Optimizer is starting to round. Tried brentq, bisect and fsolve
        x_vals = np.array([9, 11])
        tg = tgeo(tot_obs=[1,10], n_samp=[1/.9, 1/.99])
        tg.pmf(0)
        pred_vals = np.round(tg.var['x'], decimals=4)
        print pred_vals
        self.assertTrue(np.array_equal(x_vals, pred_vals))

        # Test a case that was failing for Erica Newman
        x_val = [.9896]
        tg = tgeo(tot_obs=341, n_samp=4)
        tg.pmf(0)
        print tg.var['x']
        self.assertTrue(np.round(tg.var['x'], 4) == x_val[0])

        # Test that pdf and cdf give correct values
        check = dist.pmf([1,1,2,3,4,5,12,34,65])
        self.assertTrue(dist.cdf(0)[0][0] == dist.pmf(0)[0][0])
        self.assertTrue(dist.cdf(23)[0][0] ==
                    np.sum(dist.pmf(np.arange(0,24))[0]))

        # Test that fit provides the correct number of tot_obs. Already have
        # tested generic fit method.
        dist = tgeo().fit(self.abund_list)
        self.assertTrue(len(dist.params['tot_obs']) == 4)


    def test_mete_sar_iter(self):

        # Check mete sar against EW values
        EWsar_down = np.array([8.79, 12.37, 16.71, 21.81, 27.59, 34])
        #S = 23, N=3400, anchor_area=123, target_area=2000)
        EWsar_up = np.array([23, 26.47, 30.11, 33.92, 37.89, 42.01])
        sar = mete_sar_iter(n_samp=34, tot_obs=1122).iter_vals([(2 / 45.0)])
        spp = np.round(sar['items'], decimals=2)
        error = 0.001 * spp
        diff = np.abs(spp - EWsar_down)
        self.assertTrue(np.all(diff <= error))
        sar = mete_sar_iter(n_samp=23, tot_obs=3400).iter_vals([2000 / 123.0])
        spp = np.round(sar['items'], decimals=2)
        error = 0.005 * spp
        diff = np.abs(spp - EWsar_up)
        self.assertTrue(np.all(diff <= error))

        # Check that Exception is raises if you downscale too far
        self.assertRaises(Exception, mete_sar_iter(n_samp=12, tot_obs=100).iter_vals
                                                , downscale=8)

        # Check that the return has the correct length when a_list is None
        sar = mete_sar_iter(n_samp=34, tot_obs=1000).iter_vals(None, upscale=4
                                                                , downscale=6)
        self.assertTrue(len(sar) == 11)

        # Check that only halving or doubling results are returned when
        # non_iter=True
        sar = mete_sar_iter(n_samp=34, tot_obs=1000).iter_vals([1,2,.5,.25,5,.4],
                                                                 non_iter=True)
        self.assertTrue(len(sar) == 4)

        # Check errors are thrown
        sar = mete_sar_iter(n_samp=34, tot_obs=1000)

        # Check that fit method fits correctly with two arguments passed
        sar = mete_sar_iter().fit(self.sad, self.sar)
        self.assertTrue(sar.params['n_samp'] == 155)
        self.assertTrue(sar.params['tot_obs'] == sum(np.arange(1, 156)))

        # Check that fit method fits correctly with one argument passed
        sar = mete_sar_iter().fit(self.sad)
        self.assertTrue(sar.params['n_samp'] == 155)
        self.assertTrue(sar.params['tot_obs'] == sum(np.arange(1, 156)))

        # Test values
        #N/S:  1000/100, 5000/100, 200/10
        EWslop = np.array([.3887, .2612, 0.3140])

        # Test universal SAR curve values
        answ = []
        answ.append(mete_sar_iter(n_samp=100, tot_obs=1000).\
                                                univ_curve(num_iter=0)[0]['z'][0])
        answ.append(mete_sar_iter(n_samp=100, tot_obs=5000).\
                                                univ_curve(num_iter=0)[0]['z'][0])
        answ.append(mete_sar_iter(n_samp=10, tot_obs=200).\
                                                univ_curve(num_iter=0)[0]['z'][0])
        answ = np.array(answ)
        #Using different methods to calculate so use error
        error = 0.05 * answ
        diff = np.abs(EWslop - answ)
        self.assertTrue(np.all(diff <= error))

        # Check correct errors thrown
        sar = mete_sar_iter(n_samp=45, tot_obs=5000)
        self.assertRaises(ValueError, sar.univ_curve, direction='hello')

        # That vals method is not implemented
        self.assertRaises(NotImplementedError, sar.vals, 4)

    def test_power_law(self):

        # Check that fit produces correct result. Predicted species at 1 should
        # be greater than observed.
        sar = powerlaw().fit(self.sad, self.sar)
        g = sar.vals([1])
        self.assertTrue(np.round(g['items'][0], decimals=0) == 200)

        # Check that c and z exist and check values of other parameters.
        sar.params['c']; sar.params['z']
        self.assertTrue(sar.params['n_samp'] == 155)
        self.assertTrue(sar.params['tot_obs'] == sum(np.arange(1, 156)))

        # Check universal curve has constant values
        uni, sar_stuff = sar.univ_curve()
        self.assertTrue(len(np.unique(np.round(uni['z'], decimals=5))) == 1)

        # Check that passing in param changes multiplier and correct error are
        # thrown
        self.assertRaises(ValueError, sar.univ_curve, direction='asf')
        self.assertRaises(AssertionError, sar.univ_curve, param='apple')
        res1, sar_st = sar.univ_curve(num_iter=5, param='tot_obs')
        res2, st_stf = sar.univ_curve(num_iter=5, param=None)
        self.assertTrue(not(np.array_equal(res1['x_over_y'],
                                                            res2['x_over_y'])))
        self.assertTrue((np.array_equal(res1['z'], res2['z'])))

    def test_gen_sar(self):
        '''Testing that this actually works'''

        # Testing that gen_sar actually runs.  Not sure what values to test it
        # against.

        # Test that mete_sar_iter and mete_sar iterated return similar values.
        # They will be slightly different because approximations are used in
        # mete_sar_iter.
        msi = mete_sar_iter(tot_obs=600, n_samp=40).iter_vals([1,2,.8,.2,.1])
        mete_sar = gen_sar(logser_ut(), tgeo(), tot_obs=600, n_samp=40)
        ms = mete_sar.iter_vals([1,2,.8,.2,.1])
        self.assertTrue(len(msi) == len(ms))
        error = 0.001 * msi['items']
        diff = np.abs(msi['items'] - ms['items'])
        self.assertTrue(np.all(diff <= error))

        # Test that changing the base changes the value of iter_vals
        gnsar = gen_sar(broken_stick(), binm(), n_samp=40, tot_obs=600)
        base2 = gnsar.iter_vals([1,2,.8,.2,.3], base=2)
        base3 = gnsar.iter_vals([1,2,.8,.2,.3], base=3)
        self.assertTrue(not(np.array_equal(base2['area'], base3['area'])))


        # Test that non_iter=False, returns only areas that match a_list
        a_list1 = [1,2,.5,.25,.1]
        a_list2 = [1,2.1,.4,.1]
        a_list3 = [3,(1/3.), .35]
        areas = gnsar.iter_vals(a_list1, base=2, non_iter=True)['area']
        self.assertTrue(set(areas) == set([1,2,.5,.25]))
        areas = gnsar.iter_vals(a_list2, base=2, non_iter=True)['area']
        self.assertTrue(set(areas) == set([1]))
        areas = gnsar.iter_vals(a_list3, base=3, non_iter=True)['area']
        self.assertTrue(set(areas) == set([3, (1/3.)]))

        # Test if I don't specify an a_list I can still upscale and down scale
        sar_arr = gnsar.iter_vals(downscale=1, upscale=1)
        self.assertTrue(len(sar_arr) == 3)
        # Non_iter should be overridden
        sar_arr = gnsar.iter_vals(downscale=1, upscale=1, non_iter=True)
        self.assertTrue(len(sar_arr) == 3)


        # Test vals and fit performs properly.  Test that fit ignores all args
        # but the first one too.
        sar = gen_sar(logser(), geo()).fit(self.sad, 4, [456,678])
        g = sar.vals([.001,.04, .5, 1])
        self.assertTrue(np.round(g['items'][3], decimals=0) == 155)
        sar2 = gen_sar(logser(), geo())
        sar2.params['n_samp'] = sar.params['n_samp']
        sar2.params['tot_obs'] = sar.params['tot_obs']
        s1 = sar.vals([.5]); s2 = sar2.vals([0.5])
        self.assertTrue(s1['items'][0] == s2['items'][0])

        # Test universal curve can iterate and not iterate
        itera = gnsar.univ_curve(num_iter=3, iterative=True)
        non_iter = gnsar.univ_curve(num_iter=3)
        self.assertTrue(len(itera) == len(non_iter))

    #More testing should be done
    def test_theta(self):

        # Testing assertions
        self.assertRaises(AssertionError, theta(n_samp=34, tot_obs=300, E=4000,
                                                                n=301).pdf, 1)

        # Test lambda values
        tht = theta(n_samp=4, tot_obs=4*4, E=4*4*4, n=10)
        pdf = tht.pdf(1)
        self.assertTrue(np.round(tht.var['lambda_2'][0], decimals=3) == 0.083)
        S = 16
        N = S * (2**8)
        E = N * (2**10)
        tht = theta(n_samp=S, tot_obs=N, E=E, n=N - 1)
        pdf = tht.pdf(1)
        self.assertTrue(np.round(tht.var['lambda_2'][0], decimals=8) == 3.82e-6)

        # Test cdf
        self.assertTrue(tht.cdf(E)[0][0] == 1)
        self.assertTrue(tht.cdf(1)[0][0] == 0)

        # Test rad doesn't throw error
        tht.rad()

        # TODO: test fit

    def test_psi(self):

        # Test lagrange values

        S = 4
        N = S * (2**4)
        E = N * (2**2)
        ps = psi(n_samp=S, tot_obs=N, E=E)
        pdf  = ps.pdf(1)
        self.assertTrue(np.round(ps.var['beta'][0] - ps.var['lambda_2'][0],
                                                        decimals=3) == -0.030)
        S = 16
        N = S * (2**8)
        E = N * (4)
        ps = psi(n_samp=S, tot_obs=N, E=E)
        pdf = ps.pdf(1)
        self.assertTrue(np.round(ps.var['beta'][0] - ps.var['lambda_2'][0],
                                                    decimals=5) == -0.00089)
        self.assertTrue(np.round(ps.var['lambda_2'], decimals=4) == 0.0013)

        # Test cdf
        self.assertTrue(np.round(ps.cdf(E)[0][0], decimals=0) == 1)
        self.assertTrue(ps.cdf(1)[0][0] == 0)

        # Test rad doesn't throw an error
        ps.rad()

    def test_nu(self):

        # Test error is raised when pdf called
        self.assertRaises(NotImplementedError, nu(n_samp=30, tot_obs=400,
                            E=5000).pdf, 0)


        # Check that pmf sums to one with appropriate values of e
        E = 5000
        nudist = nu(tot_obs=500, n_samp=50, E=E)
        pmf = nudist.pmf(np.arange(1.1, E + .1, step=.1))[0]
        self.assertTrue(np.round(sum(.1 * pmf), decimals=1) == 1.0)

        # Value with no support should equal 0
        self.assertTrue(nudist.pmf(1)[0][0] == 0)
        self.assertTrue(nudist.cdf(1)[0][0] == 0)

        #Check that the last value in cdf is 1
        self.assertTrue(np.round(nudist.cdf(E)[0][0], decimals=1) == 1)

        # Max support should equal 1
        l2 = nudist.var['lambda_2'][0]
        e_max = 1 + (1 / l2)
        self.assertTrue(np.round(nudist.cdf(e_max), decimals=1) == 1)


        # Test that rad works
        g = nudist.rad()
        self.assertTrue((len(g[0]) == 50))

        # Test fit
        g = nu().fit([([1,2,3,4,5,6,7], [1,2,3,4,5,6,7])])
        self.assertTrue(g.params['tot_obs'][0] == 28)
        self.assertTrue(g.params['n_samp'][0] == 7)
        self.assertTrue(g.params['E'][0] == 28)

if __name__ == '__main__':
    unittest.main()











