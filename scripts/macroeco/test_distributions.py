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
        pmf, var = logser(n_samp=240, tot_obs=15609).pmf(1)
        self.assertTrue(np.round(var['p'][0], decimals=4) == 0.9974)

        # Test cdf reaches 1
        cdf = np.round(logser(n_samp=45, tot_obs=1200).cdf(1200)[0][0][0], 
                                                                    decimals=1)
        self.assertTrue(cdf == 1)

        # Test known value of cdf when p = 0.90335
        known_cdf = 0.737623 # From scipy.stats.logser
        pred_cdf = logser(n_samp=14, tot_obs=56).cdf(4)[0][0][0]
        self.assertTrue(np.round(pred_cdf, decimals = 6) == known_cdf)


    def test_logser_ut(self):
        # Test error raising
        self.assertRaises(AssertionError, logser_ut(n_samp=234, tot_obs=67).pmf, 1)
        self.assertRaises(AssertionError, logser_ut(n_samp=34, tot_obs=0).pmf, 1)

        # Test that pmf is correct length
        pmf = logser_ut(n_samp=34, tot_obs=567).pmf(np.arange(1, 568))[0][0]
        self.assertTrue(len(pmf) == 567)

        # Test that values equal values from John's book (Harte 2011)
        pmf, x = logser_ut(n_samp=4, tot_obs=4 * 4).pmf(1)
        self.assertTrue(np.round(-np.log(x['x'][0]), decimals=4) == 0.0459)
        pmf, x = logser_ut(n_samp=4, tot_obs=2**4 * 4).pmf(1)
        self.assertTrue(np.round(-np.log(x['x'][0]), decimals=5) == -0.00884)
        pmf, x = logser_ut(n_samp=4, tot_obs=2**8 * 4).pmf(1)
        self.assertTrue(np.round(-np.log(x['x'][0]), decimals=5) == -0.00161)
        pmf, x = logser_ut(n_samp=16, tot_obs=2**8 * 16).pmf(1)
        self.assertTrue(np.round(-np.log(x['x'][0]), decimals=6) == 0.000413)
        pmf, x = logser_ut(n_samp=64, tot_obs=2**12 * 64).pmf(1)
        self.assertTrue(np.round(-np.log(x['x'][0]), decimals=7) == 0.0000228)
        
        # Check that they don't fail
        logser_ut(n_samp=64, tot_obs=1000).rad()
        logser_ut(n_samp=64, tot_obs=1000).cdf((1,1,2,4,5,7,12))
        
        # Test correct answer when n_samp == tot_obs
        pmf, x = logser_ut(n_samp=31, tot_obs=31).pmf([1,2,3,4,5])
        self.assertTrue(x['x'][0] == 0)
        self.assertTrue(np.array_equal(pmf[0], np.array([1,0,0,0,0])))


    def test_logser_ut_appx(self):
        # Test error raising
        self.assertRaises(AssertionError, logser_ut_appx(n_samp=234, 
                                                            tot_obs=67).pmf, 1)
        self.assertRaises(AssertionError, logser_ut_appx(n_samp=34, 
                                                             tot_obs=0).pmf, 1)

        # Test that values equal values from John's book (Harte 2011)
        pmf, x = logser_ut_appx(n_samp=4, tot_obs=4 * 4).pmf(1)
        self.assertTrue(np.round(-np.log(x['x'][0]), decimals=3) == 0.116)
        pmf, x = logser_ut_appx(n_samp=4, tot_obs=2**4 * 4).pmf(1)
        self.assertTrue(np.round(-np.log(x['x'][0]), decimals=4) == 0.0148)
        pmf, x = logser_ut_appx(n_samp=4, tot_obs=2**8 * 4).pmf(1)
        self.assertTrue(np.round(-np.log(x['x'][0]), decimals=6) == 0.000516)
        pmf, x = logser_ut_appx(n_samp=16, tot_obs=2**8 * 16).pmf(1)
        self.assertTrue(np.round(-np.log(x['x'][0]), decimals=6) == 0.000516)
        pmf, x = logser_ut_appx(n_samp=64, tot_obs=2**12 * 64).pmf(1)
        self.assertTrue(np.round(-np.log(x['x'][0]), decimals=7) == 0.0000229
                        or
                        np.round(-np.log(x['x'][0]), decimals=7) == 0.0000228)

        pmf, x = logser_ut_appx(n_samp=31, tot_obs=31).pmf([1,2,3,4,5])
        self.assertTrue(x['x'][0] == 0)
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
            pmf = sum(plognorm(mu=mu, sigma=sigma).pmf(sad)[0][0])
            diff1 = np.round(R[i], decimals=5) - np.round(pmf, decimals=5)
            self.assertTrue(abs(diff1) == 0)

        # Test pmf is zero when mu or sigma negative
        self.assertTrue(sum(np.round(plognorm(mu=-3,sigma=3).\
                                     pmf([1,2,3,4,5])[0][0], decimals=3)) == 0) 
        self.assertTrue(sum(np.round(plognorm(mu=3,sigma=-3).\
                                     pmf([1,2,3,4,5])[0][0], decimals=3)) == 0)

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
        # TODO: No test below - should test pmf and cdf, at minimum

        #Test our pmf against R's poilog
        R_zero_trun = [0.11620, 0.07216, 0.05201, 0.04049, 0.02783, 0.02398,
                       0.00686]
        pred_plog = plognorm_lt(mu=2, sigma=3).pmf([1,2,3,4,6,7,23])[0][0] 
        self.assertTrue(np.array_equal(R_zero_trun, np.round(pred_plog,
                                                                  decimals=5)))

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

        rcdf = np.array([0.3319, 0.3319, 0.4869, 0.5127, 0.6124])
        lnorm = np.round(lognorm(tot_obs=np.exp(3), n_samp=1, sigma=2).\
                                    pmf(np.arange(1,11))[0][0], decimals=4)
        diff = r_output - lnorm
        self.assertTrue(np.all(diff == 0))

        # Test cdf against R cdf
        pycdf = np.round(lognorm(tot_obs=np.exp(1.5 + (3.45 / 2)), n_samp=1, 
                            sigma=3.45).cdf([1,1,4,5,12])[0][0], decimals=4)
        diff = rcdf - pycdf
        self.assertTrue(np.all(diff == 0))
        
        # Test that these don't fail
        lognorm().fit([self.abund_list[0]])
        tot_obs=sum(self.abund_list[0])
        n_samp=len(self.abund_list[0])
        lognorm(tot_obs=tot_obs, n_samp=n_samp, mu=2, sigma=5).rad()

        # Test fit parameter length is correct
        # TODO: Test fit
        dist = lognorm().fit(self.abund_list)
        dist.pmf(3)
        dist.pmf([[3],[4],[5],[6]])
        self.assertTrue(len(dist.params['tot_obs']) == 4) 

      
    def test_geo_ser(self):
        # TODO: Test pmf, cdf

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
        diff = np.floor(np.array(mag_pred_sad)) - np.floor(geo_sad)
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
                            pmf(np.arange(1, 17))[0][0] * n_samp, decimals=3)
        diff = np.array(expt) - bs
        self.assertTrue(np.all(diff == 0))

        # Test that these don't fail 
        broken_stick(n_samp=23, tot_obs=500).cdf([1,2,500])
        broken_stick(n_samp=23, tot_obs=500).rad()
        broken_stick().fit(self.abund_list)
    
    
    def test_sugihara(self):
        # Data from A. J. Baczkowski 1997
        sugi = [.4761, .1962, .1180, .0751, .0499, .0337, .0226, .0148,\
                     .0090, .0047]

        # Test rad
        # Passes with 2% error regularly. Being conservative with 5% error
        ps = sugihara(n_samp=10, tot_obs=400).rad(sample_size=20000)[0] / 400
        error = .05 * ps
        diff = np.array(sugi) - ps
        ind = np.abs(diff) <= error
        self.assertTrue(np.all(ind))

        # Test that error is raised for cdf and pdf methods
        self.assertRaises(NotImplementedError, sugihara().pmf, 67)
        self.assertRaises(NotImplementedError, sugihara().cdf, 34)

     
    def test_binm(self):
        
        # Check that pdf and cdf give correct answers
        dist = binm(tot_obs=8123, n_samp=10)
        self.assertTrue(dist.cdf(8123)[0][0][0] == 1)
        self.assertTrue(dist.cdf(0) == dist.pmf(0))
        self.assertTrue(dist.cdf(1)[0][0][0] == (sum(dist.pmf([0,1])[0][0])))

        # Check that appropriate errors are raised
        self.assertRaises(TypeError, dist.pmf, [[1], [1]])
        self.assertRaises(TypeError, binm(g=3, N=45).pmf, (23,24))

        # Check that fit works
        dist = binm().fit(self.abund_list)
     
    def test_pois(self):

        # Check that pdf and cdf give correct answers
        dist = pois(tot_obs=112, n_samp=20)
        self.assertTrue(dist.cdf(112)[0][0][0] == 1)
        a = np.round(dist.cdf(0)[0][0][0], decimals=12)
        b = np.round(dist.pmf(0)[0][0][0], decimals=12)
        self.assertTrue(a == b)
        a = np.round(dist.cdf(23)[0][0], decimals=12)

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
        self.assertTrue(np.round(dist.cdf(2300)[0][0][0], decimals=1) == 1.0)

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
        self.assertTrue(np.round(dist.cdf(2300)[0][0][0], decimals=1) == 1.0)

        # Check that k of length one is extended to length 2
        dist = nbd_lt(tot_obs=[400, 600], n_samp=[30, 23], k=[3])
        pmf, var = dist.pmf(1)
        self.assertTrue(np.array_equal(var['k'], np.array([3,3]))) 

        # Multiple entries both yield cdf with 1
        dist = nbd_lt(tot_obs=[400, 600], n_samp=[30, 23], k=[3,2])
        cdf = dist.cdf([[400], [600]])
        a = np.round(cdf[0][0][0], decimals=1)
        b = np.round(cdf[0][0][0], decimals=1)
        self.assertTrue(a == b)


    def test_fnbd(self):

        # Test that no error is thrown if a zero is passed
        fnbd().fit([[0,1,2,3,4,5,6]])
        
        # TypeError if k is not given
        dist = fnbd(tot_obs=2300, n_samp=20)
        self.assertRaises(TypeError, dist.pmf, 45)

        # Test that cdf sums to one
        dist = fnbd(tot_obs=2300, n_samp=15, k=2)
        self.assertTrue(np.round(dist.cdf(2300)[0][0][0], decimals=1) == 1.0)

        # Test that cdf and pdf at 0 give the same answer
        a = np.round(dist.pmf(0)[0][0][0], decimals=12)
        b = np.round(dist.cdf(0)[0][0][0], decimals=12)
        self.assertTrue(a == b)

        # Test that fit of of geometric data gives k=1
        np.random.seed(345)
        p = 0.001
        geo_data = np.random.geometric(p, size=10000)
        dist = fnbd().fit([geo_data])
        self.assertTrue(np.round(dist.params['k'][0], decimals=1) == 1)

    
    def test_tgeo(self):

        # Test tgeo cdf is one
        dist = tgeo(n_samp=10, tot_obs=2345)
        self.assertTrue(np.round(dist.cdf(2345)[0][0][0], decimals=1) == 1.0)

        #It would be good to test against values in Harte book.

        # Test that pdf and cdf give correct values
        check = dist.pmf([1,1,2,3,4,5,12,34,65])
        self.assertTrue(dist.cdf(0)[0][0][0] == dist.pmf(0)[0][0][0])
        self.assertTrue(dist.cdf(23)[0][0][0] == 
                    np.sum(dist.pmf(np.arange(0,24))[0][0]))

        # Test that fit provides the correct number of tot_obs
        dist = tgeo().fit(self.abund_list)
        self.assertTrue(len(dist.params['tot_obs']) == 4)
    
    
    def test_mete_sar_iter(self):
        
        # Check mete sar against EW values
        EWsar_down = np.array([8.79, 12.37, 16.71, 21.81, 27.59, 34])
        #S = 23, N=3400, anchor_area=123, target_area=2000)
        EWsar_up = np.array([23, 26.47, 30.11, 33.92, 37.89, 42.01])
        sar = mete_sar_iter(n_samp=34, tot_obs=1122).vals([(2 / 45.0)])
        spp = np.round(sar['items'], decimals=2)
        error = 0.001 * spp
        diff = np.abs(spp - EWsar_down)
        self.assertTrue(np.all(diff <= error))
        sar = mete_sar_iter(n_samp=23, tot_obs=3400).vals([2000 / 123.0])
        spp = np.round(sar['items'], decimals=2)
        error = 0.005 * spp
        diff = np.abs(spp - EWsar_up)
        self.assertTrue(np.all(diff <= error))

        # Check that Exception is raises if you downscale too far
        self.assertRaises(Exception, mete_sar_iter(n_samp=12, tot_obs=100).vals
                                                , None, downcale=8)

        # Check that the return has the correct length when a_list is None
        sar = mete_sar_iter(n_samp=34, tot_obs=1000).vals(None, upscale=4
                                                                , downscale=6)
        self.assertTrue(len(sar) == 11)

        # Check that only halving or doubling results are returned when 
        # non_iter=True
        sar = mete_sar_iter(n_samp=34, tot_obs=1000).vals([1,2,.5,.25,5,.4],
                                                                 non_iter=True)
        self.assertTrue(len(sar) == 4) 

        # Check errors are thrown
        sar = mete_sar_iter(n_samp=34, tot_obs=1000)
        self.assertRaises(TypeError, sar.vals, .1)

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
                                                univ_curve(num_iter=0)['z'][0])
        answ.append(mete_sar_iter(n_samp=100, tot_obs=5000).\
                                                univ_curve(num_iter=0)['z'][0])
        answ.append(mete_sar_iter(n_samp=10, tot_obs=200).\
                                                univ_curve(num_iter=0)['z'][0])
        answ = np.array(answ)
        #Using different methods to calculate so use error
        error = 0.05 * answ
        diff = np.abs(EWslop - answ)
        self.assertTrue(np.all(diff <= error))

        # Check correct errors thrown
        sar = mete_sar_iter(n_samp=45, tot_obs=5000)
        self.assertRaises(ValueError, sar.univ_curve, direction='hello')

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
        uni = sar.univ_curve()
        self.assertTrue(len(np.unique(np.round(uni['z'], decimals=5))) == 1)

        # Check that passing in param changes multiplier and correct error are
        # thrown
        self.assertRaises(ValueError, sar.univ_curve, direction='asf')
        self.assertRaises(AssertionError, sar.univ_curve, param='apple')
        res1 = sar.univ_curve(num_iter=5, param='tot_obs')
        res2 = sar.univ_curve(num_iter=5, param=None)
        self.assertTrue(not(np.array_equal(res1['x_over_y'],
                                                            res2['x_over_y'])))
        self.assertTrue((np.array_equal(res1['z'], res2['z'])))
        
    def test_gen_sar(self):
        '''Testing that this actually works'''
        
        # Testing that gen_sar actually runs.  Not sure what values to test it
        # against.
        sar = gen_sar(logser(), geo()).fit(self.sad)
        g = sar.vals([.001,.04, .5, 1])
        self.assertTrue(np.round(g['items'][3], decimals=0) == 155)
        sar2 = gen_sar(logser(), geo())
        sar2.params['n_samp'] = sar.params['n_samp']
        sar2.params['tot_obs'] = sar.params['tot_obs']
        sar2.params['sad_pmf'] = sar.params['sad_pmf']
        s1 = sar.vals([.5]); s2 = sar2.vals([0.5])
        self.assertTrue(s1['items'][0] == s2['items'][0])

        # Test univ_curve upscales and downscales
        #sar3 = gen_sar(broken_stick(), binm()).fit(self.sad)
        #sar3.univ_curve(num_iter=2, direction='up')
        #sar3.univ_curve(num_iter=2, direction='down')

    #More testing should be done
    def test_theta(self):

        # Testing assertions
        self.assertRaises(AssertionError, theta(n_samp=34, tot_obs=300, E=4000,
                                                                n=301).pdf, 1)

        # Test lambda values
        tht = theta(n_samp=4, tot_obs=4*4, E=4*4*4, n=10)
        pdf, var = tht.pdf(1)
        self.assertTrue(np.round(var['l2'][0], decimals=3) == 0.083)
        S = 16
        N = S * (2**8)
        E = N * (2**10)
        tht = theta(n_samp=S, tot_obs=N, E=E, n=N - 1)
        pdf, var = tht.pdf(1)
        self.assertTrue(np.round(var['l2'][0], decimals=8) == 3.82e-6)

        # Test cdf
        self.assertTrue(tht.cdf(E)[0][0][0] == 1)
        self.assertTrue(tht.cdf(1)[0][0][0] == 0)

        # Test rad doesn't throw error
        tht.rad()

        # TODO: test fit

    def test_psi(self):

        # Test lagrange values

        S = 4
        N = S * (2**4)
        E = N * (2**2)
        ps = psi(n_samp=S, tot_obs=N, E=E)
        pdf, var = ps.pdf(1)
        self.assertTrue(np.round(var['beta'][0] - var['l2'][0], decimals=3) 
                                                                    == -0.030)
        S = 16
        N = S * (2**8)
        E = N * (4)
        ps = psi(n_samp=S, tot_obs=N, E=E)
        pdf, var = ps.pdf(1)
        self.assertTrue(np.round(var['beta'][0] - var['l2'][0], decimals=5) 
                                                                == -0.00089)
        self.assertTrue(np.round(var['l2'], decimals=4) == 0.0013)

        # Test cdf
        print ps.cdf(E)[0][0][0]
        self.assertTrue(np.round(ps.cdf(E)[0][0][0], decimals=0) == 1)
        self.assertTrue(ps.cdf(1)[0][0][0] == 0)

        # Test rad doesn't throw an error
        ps.rad()
        



        



        







  
