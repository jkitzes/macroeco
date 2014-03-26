#!/usr/bin/python
"""
Tests for compare module

"""
from __future__ import division

from numpy.testing import (TestCase, assert_equal, assert_array_equal,
                           assert_almost_equal, assert_array_almost_equal,
                           assert_allclose, assert_, assert_raises)

from macroeco.compare import *
import numpy as np
import scipy.stats as stats
import numpy.testing as nt


class TestNLL(TestCase):
    '''Test NLL in compare'''

    def test_nll(self):

        # Test against R result: sum(dnorm(c(1,2,3,4,5), log=TRUE))
        R_res = 32.09469
        data = np.array([1, 2, 3, 4, 5])
        model = stats.norm(loc=0, scale=1)
        lglk = nll(data, model)
        assert_almost_equal(R_res, lglk, decimal=5)


class TestAIC(TestCase):
    """Test AIC function"""

    def test_aic_basic(self):
        """Testing basic functionality of AIC"""

        # Test case 1
        model = stats.norm(loc=0, scale=1)
        data = np.arange(1, 9)
        aic1 = AIC(data, model, corrected=False)
        expected = 222.703016531  # Calculated by hand
        assert_almost_equal(aic1, expected)

        # Test case 2
        model = stats.gamma(a=2)
        data = [1, 1, 1, 2, 4, 5, 7, 12]
        aic1 = AIC(data, model, corrected=False)
        expected = 51.760607494
        assert_almost_equal(aic1, expected, decimal=6)

        # Test case 3
        model = stats.gamma(a=2, loc=0)
        aic1 = AIC(data, model, corrected=False)
        expected = 53.760607494
        assert_almost_equal(aic1, expected, decimal=6)

    def test_aic_given_params(self):
        """ Test AIC if params are given """

        # Test case 1
        model = stats.norm()
        data = np.arange(1, 9)
        aic1 = AIC(data, model, corrected=False, params=2)

        # statsmodel.tools.eval_measures.aic: aic(L, 8, 2)
        expected = 222.703016531
        assert_almost_equal(aic1, expected)

        # Test case 2
        model = stats.gamma(2)
        data = [1, 1, 1, 2, 4, 5, 7, 12]
        aic1 = AIC(data, model, corrected=False, params=1)

        # statsmodel.tools.eval_measures.aic: aic(L, 8, 1)
        expected = 51.760607494
        assert_almost_equal(aic1, expected, decimal=6)

        # Test case 3
        model = stats.gamma(2, 0)
        aic1 = AIC(data, model, corrected=False, params=2)

        # statsmodel.tools.eval_measures.aic: aic(L, 8, 2)
        expected = 53.760607494
        assert_almost_equal(aic1, expected, decimal=6)

    def test_aicc(self):
        """ Test AICC gives expected results"""

        # Test values
        model = stats.norm()
        data = np.arange(1, 9)
        aic1 = AIC(data, model, corrected=True, params=2)
        expected = 225.10302  # Calculated by hand
        assert_almost_equal(expected, aic1, decimal=5)


class TestAICWeights(TestCase):

    def test_aic_weights(self):

        # Test values
        data = [1, 1, 1, 2, 3, 4, 7, 23, 78]
        models = [stats.norm(scale=100), stats.norm(scale=99)]
        aic_vals = [AIC(data, tm) for tm in models]

        aicw, delta_aic = AIC_weights(aic_vals)

        # Calculated by hand
        pred = np.array([0.47909787, 0.52090213])
        assert_array_almost_equal(aicw, pred)


class TestRsquared(TestCase):

    def test_basic_r_squared(self):

        # Already unittested in scipy. Checking for functionaliity
        test_data = np.random.randint(5, 100, 100)
        rsq = r_squared(test_data, test_data)
        assert_equal(rsq, 1)

    def test_one_to_one_rsq(self):

        # Identical data should lead to an R^2 of 1
        test_data = np.random.randint(5, 100, 100)
        rsq = r_squared(test_data, test_data, one_to_one=True)
        assert_equal(rsq, 1)

        # Test against R^2 from fixed slope linear regression in R
        # Calculate by hand?


class TestBinData(TestCase):

    def test_bin_data_functionality(self):

        # Test against R's vegan prestonfit: prestonfit(data, tiesplit=FALSE)
        # Note that vegan drops the bins with 0 values

        data = np.array([1, 1, 1, 1, 2, 2, 4, 4, 8, 16, 17.1, 89])
        vegan = np.array([4, 2, 2, 1, 1, 1, 0, 1], dtype=np.float)
        test_res = bin_data(data, max(data))[0]
        assert_array_equal(test_res, vegan)

        data = np.array([1, 1, 1, 1, 4, 5, 6, 7, 12, 34, 56])
        vegan = np.array([4, 0, 1, 3, 1, 0, 2], dtype=np.float)
        test_res = bin_data(data, max(data))[0]
        assert_array_equal(test_res, vegan)

    def test_bin_data_boundary(self):
        # Test boundary condition
        data = np.array([1, 2])
        vegan = np.array([1, 1], dtype=np.float)
        test_res = bin_data(data, max(data))[0]
        assert_array_equal(test_res, vegan)

        data = np.array([1, 1, 1])
        vegan = np.array([3], dtype=np.float)
        test_res = bin_data(data, max(data))[0]
        assert_array_equal(test_res, vegan)

        data = np.array([1, 2, 3])
        vegan = np.array([1, 1, 1], dtype=np.float)
        test_res = bin_data(data, max(data))[0]
        assert_array_equal(test_res, vegan)

    # def test_lrt(self):

    #     # Test against what the lrtest() R function returns
    #     model1 = 158.0494
    #     model0 = 139.806
    #     R_chisquare = 36.4868
    #     R_p = 1.537e-09

    #     pred_chi, pred_p = lrt(model1, model0, 1)

    #     assert_almost_equal(pred_chi, R_chisquare)
    #     assert_almost_equal(pred_p, R_p)

    # def test_empirical_cdf(self):

    #     #Test against R's ecdf function

    #     # Test Case 1
    #     test_data = [1, 1, 1, 1, 2, 3, 4, 5, 6, 6]
    #     R_res = [.4, .4, .4, .4, .5, .6, .7, .8, 1, 1]
    #     res = empirical_cdf(test_data)
    #     assert_array_equal(R_res, res)

    #     # Test Case 2
    #     test_data = [3, 3, 3, 3]
    #     R_res = [1, 1, 1, 1]
    #     res = empirical_cdf(test_data)
    #     assert_array_equal(R_res, res)
