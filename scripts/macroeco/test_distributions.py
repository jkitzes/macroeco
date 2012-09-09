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
from distributions import *
import numpy as np

#TODO:  Need to rerun tests with new fit functions. 

class TestDistributions(unittest.TestCase):
    '''Test the functions within sad_distr.py'''

    def setUp(self):
        
        self.abund_list = []
        self.abund_list.append(np.array([1,1,2,3,4,4,7,12]))
        self.abund_list.append(np.array([1,1,1,1,1,2]))
        self.abund_list.append(np.array([3,2,2,1,5,6,4,5,7]))
        self.abund_list.append(np.array([1,1,4,5,7,89,100]))
        #R output from function dpolono in VGAM package
        self.R = [0.8453224, 1.951546, 1.040038, 0.3102524]
        #Simulated values in A. J. Baczkowski 1997 paper 
        self.sugi = [.4761, .1962, .1180, .0751, .0499, .0337, .0226, .0148,\
                     .0090, .0047]
        #Testing against weecology function: N=1122, S=34, target_area=2
        #Anchor_area=45
        #There will be some differences because we are using approximations
        #to make things faster
        self.EWsar_down = np.array([8.79, 12.37, 16.71, 21.81, 27.59, 34])
        #S = 23, N=3400, anchor_area=123, target_area=2000)
        self.EWsar_up = np.array([23, 26.47, 30.11, 33.92, 37.89, 42.01])
        #N/S:  1000/100, 5000/100, 200/10
        self.EWslop = np.array([.3887, .2612, 0.3140])
        self.sar = ([0.0004, 0.0025, 0.005, 0.01, 0.25, 1], [2.2356, 10.0175,
        15.87, 24.32, 101.25, 155])
        self.sad = np.arange(1, 156)
    
    def test_logser(self):
        self.assertRaises(AssertionError, logser(S=45, N=45).pmf, 1)
        self.assertRaises(AssertionError, logser(S=234, N=67).pmf, 1)
        self.assertRaises(AssertionError, logser(S=34, N=0).pmf, 1)
        #Test against value in Fisher's paper Fisher et al. 1943
        pmf, var = logser(S=240, N=15609).pmf(1)
        self.assertTrue(np.round(var['p'][0], decimals=4) == 0.9974)
        cdf = np.round(logser(S=45, N=1200).cdf(1200)[0][0][0], decimals=1)
        self.assertTrue(cdf == 1)
    
     
    #Should test against ethans psolver
    def test_plognorm(self):
        for i, sad in enumerate(self.abund_list):
            log_abund = np.log(sad)
            mu = np.mean(log_abund)
            var = np.var(log_abund, ddof=1)
            sigma = var**0.5
            pmf = sum(plognorm(mu=mu, sigma=sigma).pmf(sad)[0][0])
            diff1 = np.round(self.R[i], decimals=5) - np.round(pmf, decimals=5)
            self.assertTrue(abs(diff1) == 0)
        self.assertTrue(sum(np.round(plognorm(mu=-3,sigma=-3).\
                                     pmf([1,2,3,4,5])[0][0], decimals=3)) == 0) 
        plognorm().fit([self.abund_list[0]])
        N = sum(self.abund_list[0]); S = len(self.abund_list[0])
        plognorm(S=S, N=N, mu=2, sigma=2).cdf(5)

    
    def test_plognorm_lt(self):
        N = sum(self.abund_list[0]); S = len(self.abund_list[0])
        plognorm_lt(S=[S,S+2], N=[N, N+2], mu=[2,3], sigma=[2,3]).cdf(5)
        plognorm_lt(mu=2, sigma=2).pmf([2,3,4,5,23])
        plognorm_lt().fit([self.abund_list[0]])
        plognorm_lt(mu=10, sigma=1).cdf(45)
        EW_fit = {'mu' :.90, 'sigma' : 2.18}
        sad = [1,1,1,1,2,3,5,6,12,13,15,23,45,67,112]
        dist = plognorm_lt().fit([sad])
        mu = dist.params['mu'][0]; sigma = dist.params['sigma'][0]
        self.assertTrue(EW_fit['mu'] == np.round(mu, decimals=2))
        self.assertTrue(EW_fit['sigma'] == np.round(sigma, decimals=2))
    
    
    def test_logser_ut(self):
        self.assertRaises(AssertionError, logser_ut(S=45, N=45).pmf, 1)
        self.assertRaises(AssertionError, logser_ut(S=234, N=67).pmf, 1)
        self.assertRaises(AssertionError, logser_ut(S=34, N=0).pmf, 1)
        pmf = logser_ut(S=34, N=567).pmf(np.arange(1, 568))[0][0]
        self.assertTrue(len(pmf) == 567)
        #Testing that values equal values from John's book (Harte 2011)
        pmf, x = logser_ut(S=4, N=4 * 4).pmf(1)
        self.assertTrue(np.round(-np.log(x['x'][0]), decimals=4) == 0.0459)
        pmf, x = logser_ut(S=4, N=2**4 * 4).pmf(1)
        self.assertTrue(np.round(-np.log(x['x'][0]), decimals=5) == -0.00884)
        pmf, x = logser_ut(S=4, N=2**8 * 4).pmf(1)
        self.assertTrue(np.round(-np.log(x['x'][0]), decimals=5) == -0.00161)
        pmf, x = logser_ut(S=16, N=2**8 * 16).pmf(1)
        self.assertTrue(np.round(-np.log(x['x'][0]), decimals=6) == 0.000413)
        pmf, x = logser_ut(S=64, N=2**12 * 64).pmf(1)
        self.assertTrue(np.round(-np.log(x['x'][0]), decimals=7) == 0.0000228)
        logser_ut(S=64, N=1000).rad()
        logser_ut(S=64, N=1000).cdf((1,1,2,4,5,7,12))

    def test_logser_ut_appx(self):
        self.assertRaises(AssertionError, logser_ut_appx(S=45, N=45).pmf, 1)
        self.assertRaises(AssertionError, logser_ut_appx(S=234, N=67).pmf, 1)
        self.assertRaises(AssertionError, logser_ut_appx(S=34, N=0).pmf, 1)
        #Testing that values = values from John's book (Harte 2011)
        pmf, x = logser_ut_appx(S=4, N=4 * 4).pmf(1)
        self.assertTrue(np.round(-np.log(x['x'][0]), decimals=3) == 0.116)
        pmf, x = logser_ut_appx(S=4, N=2**4 * 4).pmf(1)
        self.assertTrue(np.round(-np.log(x['x'][0]), decimals=4) == 0.0148)
        pmf, x = logser_ut_appx(S=4, N=2**8 * 4).pmf(1)
        self.assertTrue(np.round(-np.log(x['x'][0]), decimals=6) == 0.000516)
        pmf, x = logser_ut_appx(S=16, N=2**8 * 16).pmf(1)
        self.assertTrue(np.round(-np.log(x['x'][0]), decimals=6) == 0.000516)
        pmf, x = logser_ut_appx(S=64, N=2**12 * 64).pmf(1)
        self.assertTrue(np.round(-np.log(x['x'][0]), decimals=7) == 0.0000229 or\
                        np.round(-np.log(x['x'][0]), decimals=7) == 0.0000228)
        logser_ut_appx(S=64, N=1000).rad()
        logser_ut_appx(S=64, N=1000).cdf((1,1,2,4,5,7,12))
    
     
    def test_sugihara_rank_abund(self):
        #Testing against A. J. Baczkowski 1997 paper values
        #Passes with 2% error regularly. Being conservative with 5% error
        ps = sugihara(S=10, N=400).rad(sample_size=20000)[0] / 400
        error = .05 * ps
        diff = np.array(self.sugi) - ps
        ind = np.abs(diff) <= error
        self.assertTrue(np.all(ind))
    
    def test_geo_ser(self):
        #Using data from Magurran (1998)
        obs_sad = [370,210,120,66,35,31,15,9,3,2,1]
        mag_pred_sad = [387.6, 213.8, 117.8, 64.5, 35.5, 19.8, 10.7, 6.2, 3.3,\
                        1.8, 1.0]
        dist = geo_ser().fit([obs_sad])
        gk = dist.params['k'][0]
        self.assertTrue(np.round(gk, decimals=3) == 0.449)
        geo_sad = np.round(geo_ser(S=len(obs_sad), N=sum(obs_sad), k=.449).\
                           rad(), decimals=1)[0]
        diff = np.floor(np.array(mag_pred_sad)) - np.floor(geo_sad)
        self.assertTrue(np.all(diff == 0))
        dist = geo_ser().fit(self.abund_list)
        self.assertTrue(len(dist.params['k']) == 4)
    
    def test_broken_stick(self):
        #Using data from Magurran (1998)
        self.assertRaises(TypeError, broken_stick(S=12, N=111).pmf,[[2],[3]])
        obs_sad = [103,115,13,2,67,36,51,8,6,61,10,21,7,65,4,49,92,37,16,6,23,\
                   9,2,6,5,4,1,3,1,9,2]
        expt = [1.077, 1.040, 1.004, .970, .937, .904, .873, .843, .814,\
                    .786, .759, .732, .707, .683, .659, .636]
        S = len(obs_sad)
        N = sum(obs_sad)
        bs = np.round(broken_stick(S=S, N=N).pmf(np.arange(1, 17))[0][0] * S,\
                      decimals=3)
        diff = np.array(expt) - bs
        self.assertTrue(np.all(diff == 0))
        broken_stick(S=23, N=500).cdf([1,2,500])
        broken_stick(S=23, N=500).rad()
        broken_stick().fit(self.abund_list)

    def test_lognorm(self):
        r_output = [0.1210, .0806, .0601, 0.0476, 0.0391, .0331,  0.0285,\
                    0.0249, 0.0221, 0.0197]
        lnorm = np.round(lognorm(mu=2, sigma=2).pmf(np.arange(1,11))[0][0],
                                                                    decimals=4)
        diff = r_output - lnorm
        self.assertTrue(np.all(diff == 0))
        lognorm().fit([self.abund_list[0]])
        N=sum(self.abund_list[0]); S=len(self.abund_list[0])
        lognorm(N=N, S=S, mu=2, sigma=5).rad()
        lognorm(mu=1.5, sigma=3.45).cdf([1,1,4,5,12])
        dist = lognorm().fit(self.abund_list)
        dist.pmf(3)
        dist.pmf([[3],[4],[5],[6]])
        self.assertTrue(len(dist.params['mu']) == 4)
     
    def test_binm(self):
        dist = binm(tot_obs=8123, n_samp=10)
        self.assertTrue(dist.cdf(8123)[0][0][0] == 1)
        self.assertTrue(dist.cdf(0) == dist.pmf(0))
        self.assertTrue(dist.cdf(1)[0][0][0] == (sum(dist.pmf([0,1])[0][0])))
        self.assertRaises(TypeError, dist.pmf, [[1], [1]])
        self.assertRaises(AssertionError, binm(g=3, N=45).pmf, (23,24))
        dist = binm().fit(self.abund_list)
     
    def test_pois(self):
        dist = pois(N=112, n_samp=20)
        self.assertTrue(dist.cdf(112)[0][0][0] == 1)
        a = np.round(dist.cdf(0)[0][0][0], decimals=12)
        b = np.round(dist.pmf(0)[0][0][0], decimals=12)
        self.assertTrue(a == b)
        a = np.round(dist.cdf(23)[0][0], decimals=12)
        dist = pois().fit(self.abund_list)
        rads = dist.rad()
        self.assertTrue(len(rads) == 4)

    def test_nbd(self):
        dist = nbd(N=2300, n_samp=45)
        self.assertRaises(AssertionError, dist.pmf, 45)
        dist = nbd(N=2300, n_samp=24, k=2)
        self.assertTrue(np.round(dist.cdf(2300)[0][0][0], decimals=1) == 1.0)
        np.random.seed(345)
        p = 0.001
        geo_data = np.random.geometric(p, size=10000)
        dist = nbd().fit([geo_data])
        self.assertTrue(np.round(dist.params['k'][0], decimals=1) == 1)

    def test_fnbd(self):
        dist = fnbd(N=2300, n_samp=20)
        self.assertRaises(AssertionError, dist.pmf, 45)
        dist = fnbd(N=2300, n_samp=15, k=2)
        self.assertTrue(np.round(dist.cdf(2300)[0][0][0], decimals=1) == 1.0)
        a = np.round(dist.pmf(0)[0][0][0], decimals=12)
        b = np.round(dist.cdf(0)[0][0][0], decimals=12)
        self.assertTrue(a == b)
        np.random.seed(345)
        p = 0.001
        geo_data = np.random.geometric(p, size=10000)
        dist = fnbd().fit([geo_data])
        self.assertTrue(np.round(dist.params['k'][0], decimals=1) == 1)

    
    def test_tgeo(self):
        dist = tgeo(n_samp=10, tot_obs=2345)
        self.assertTrue(np.round(dist.cdf(2345)[0][0][0], decimals=1) == 1.0)
        #It would be good to test against values in Harte book.
        check = dist.pmf([1,1,2,3,4,5,12,34,65])
        self.assertTrue(dist.cdf(0)[0][0][0] == dist.pmf(0)[0][0][0])
        self.assertTrue(dist.cdf(23)[0][0][0] == 
                    np.sum(dist.pmf(np.arange(0,24))[0][0]))
        dist = tgeo().fit(self.abund_list)
        self.assertTrue(len(dist.params['tot_obs']) == 4)

    
    def test_mete_sar_iter(self):
        sar = mete_sar_iter(S=34, N=1122).vals([(2 / 45.0)])
        spp = np.round(sar['items'], decimals=2)
        error = 0.001 * spp
        diff = np.abs(spp - self.EWsar_down)
        self.assertTrue(np.all(diff <= error))
        sar = mete_sar_iter(S=23, N=3400).vals([2000 / 123.0])
        spp = np.round(sar['items'], decimals=2)
        error = 0.005 * spp
        diff = np.abs(spp - self.EWsar_up)
        self.assertTrue(np.all(diff <= error))
        self.assertRaises(Exception, mete_sar_iter(S=12, N=100).vals,
                                                None, downcale=8)
        sar = mete_sar_iter(S=34, N=1000).vals(None, upscale=4, downscale=6)
        self.assertTrue(len(sar) == 11)
        sar = mete_sar_iter().fit(self.sar, self.sad)
        self.assertTrue(sar.params['S'] == 155)
        self.assertTrue(sar.params['N'] == sum(np.arange(1, 156)))
        sar.vals([.4,3,1,.8])


    def test_power_law(self):
        sar = powerlaw().fit(self.sar, self.sad)
        g = sar.vals([1])
        self.assertTrue(np.round(g['items'][0], decimals=3) == 199.808)
        sar.params['c']; sar.params['z']
        self.assertTrue(sar.params['S'] == 155)
        self.assertTrue(sar.params['N'] == sum(np.arange(1, 156)))

    def test_universal_sar(self):
        '''Interesting test results.  EW and our functions agree that using
        sar_method1 METE is not quite universal.  The z that results from
        different N and S with the same ratio differe slightly.  The amount
        they differ might just be a result of approximations.  Need to look
        into this'''
        
        answ = []
        answ.append(mete_sar_iter(S=100, N=1000).univ_curve(num_iter=0)['z'][0])
        answ.append(mete_sar_iter(S=100, N=5000).univ_curve(num_iter=0)['z'][0])
        answ.append(mete_sar_iter(S=10, N=200).univ_curve(num_iter=0)['z'][0])
        answ = np.array(answ)
        #Using different methods to calculate so use error
        error = 0.05 * answ
        diff = np.abs(self.EWslop - answ)
        self.assertTrue(np.all(diff <= error))

    def test_gen_sar(self):
        '''Testing that this actually works'''

        sar = gen_sar(logser(), geo()).fit(self.sar, self.sad)
        g = sar.vals([.001,.04, .5, 1])
        self.assertTrue(np.round(g['items'][3], decimals=0) == 155)
        sar2 = gen_sar(logser(), geo())
        sar2.params['S'] = sar.params['S']
        sar2.params['N'] = sar.params['N']
        sar2.params['sad_pmf'] = sar.params['sad_pmf']
        s1 = sar.vals([.5]); s2 = sar2.vals([0.5])
        self.assertTrue(s1['items'][0] == s2['items'][0])

