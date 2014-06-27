from __future__ import division

from numpy.testing import (TestCase, assert_equal, assert_array_equal,
                           assert_almost_equal, assert_array_almost_equal,
                           assert_allclose, assert_, assert_raises)

from macroeco.compare import *
import numpy as np
import scipy.stats as stats
import macroeco.models as mod


class TestNLL(TestCase):

    def test_nll(self):
        # R: sum(dnorm(c(1,2,3,4,5), log=TRUE))
        R_res = 32.09469
        data = np.array([1, 2, 3, 4, 5])
        model = stats.norm(loc=0, scale=1)
        lglk = nll(data, model)
        assert_almost_equal(R_res, lglk, decimal=5)


# TODO: Test LRT


class TestAIC(TestCase):

    def test_aic_basic(self):

        model = stats.norm(loc=0, scale=1)
        data = np.arange(1, 9)
        aic1 = AIC(data, model, corrected=False)
        expected = 222.703016531  # Calculated by hand
        assert_almost_equal(aic1, expected, decimal=6)

        model = stats.gamma(a=2)
        data = [1, 1, 1, 2, 4, 5, 7, 12]
        aic1 = AIC(data, model, corrected=False)
        expected = 51.760607494   # Calculated by hand
        assert_almost_equal(aic1, expected, decimal=6)

        model = stats.gamma(a=2, loc=0)
        aic1 = AIC(data, model, corrected=False)
        expected = 53.760607494   # Calculated by hand
        assert_almost_equal(aic1, expected, decimal=6)

    def test_aic_given_params(self):

        model = stats.norm()
        data = np.arange(1, 9)
        aic1 = AIC(data, model, corrected=False, params=2)
        # statsmodel.tools.eval_measures.aic: aic(L, 8, 2)
        expected = 222.703016531
        assert_almost_equal(aic1, expected)

        model = stats.gamma(2)
        data = [1, 1, 1, 2, 4, 5, 7, 12]
        aic1 = AIC(data, model, corrected=False, params=1)
        # statsmodel.tools.eval_measures.aic: aic(L, 8, 1)
        expected = 51.760607494
        assert_almost_equal(aic1, expected, decimal=6)

        model = stats.gamma(2, 0)
        aic1 = AIC(data, model, corrected=False, params=2)
        # statsmodel.tools.eval_measures.aic: aic(L, 8, 2)
        expected = 53.760607494
        assert_almost_equal(aic1, expected, decimal=6)

    def test_aicc(self):

        model = stats.norm()
        data = np.arange(1, 9)
        aic1 = AIC(data, model, corrected=True, params=2)
        expected = 225.10302  # Calculated by hand
        assert_almost_equal(expected, aic1, decimal=5)


class TestFullModelNLL(TestCase):

    def test_correct_value_for_continuous_models(self):

        # Test that the full model returns what we expect
        data = np.array([3, 4, 5])

        models = [mod.lognorm]
        for model in models:

            params = [model.fit_mle(np.array([td])) for td in data]
            values = [model(*params[i]).logpdf(data[i]) for i in
                            xrange(len(data))]
            pred_nll = -np.sum(values)

            test_nll = full_model_nll(data, model)

            assert_equal(pred_nll, test_nll)

    def test_correct_value_for_discrete_models(self):

        # Test that the full model returns what we expect
        data = np.array([3, 4, 5])

        models = [mod.nbinom]
        for model in models:

            params = [model.fit_mle(np.array([td])) for td in data]
            values = [model(*params[i]).logpmf(data[i]) for i in
                            xrange(len(data))]
            pred_nll = -np.sum(values)

            test_nll = full_model_nll(data, model)

            assert_equal(pred_nll, test_nll)


class TestAICCompare(TestCase):

    def test_aic_delta_and_weights(self):

        data = [1, 1, 1, 2, 3, 4, 7, 23, 78]
        models = [stats.norm(scale=100), stats.norm(scale=99)]
        aic_vals = [AIC(data, tm) for tm in models]
        daic, aicw = AIC_compare(aic_vals)

        pred = np.array([0.47909787, 0.52090213])   # Calculated by hand
        assert_array_almost_equal(aicw, pred)
        assert_array_almost_equal(daic, [daic[0]-daic[1], 0])


class TestRsquared(TestCase):

    def test_r_squared_repeated_data(self):

        # Identical data should lead to an R^2 of 1
        test_data = np.random.randint(5, 100, 100)
        rsq = r_squared(test_data, test_data)
        assert_equal(rsq, 1)

    # TODO: Test known R2 for regression and one-to-one


class TestPrestonBin(TestCase):

    def test_bin_functionality(self):

        # Test against R's vegan prestonfit: prestonfit(data, tiesplit=FALSE)
        # Note that vegan drops the bins with 0 values

        data = np.array([1, 1, 1, 1, 2, 2, 4, 4, 8, 16, 17.1, 89])
        vegan = np.array([4, 2, 2, 1, 1, 1, 0, 1], dtype=np.float)
        test_res = preston_bin(data, max(data))[0]
        assert_array_equal(test_res, vegan)

        data = np.array([1, 1, 1, 1, 4, 5, 6, 7, 12, 34, 56])
        vegan = np.array([4, 0, 1, 3, 1, 0, 2], dtype=np.float)
        test_res = preston_bin(data, max(data))[0]
        assert_array_equal(test_res, vegan)

    def test_bin_data_boundary(self):
        # Test boundary condition
        data = np.array([1, 2])
        vegan = np.array([1, 1], dtype=np.float)
        test_res = preston_bin(data, max(data))[0]
        assert_array_equal(test_res, vegan)

        data = np.array([1, 1, 1])
        vegan = np.array([3], dtype=np.float)
        test_res = preston_bin(data, max(data))[0]
        assert_array_equal(test_res, vegan)

        data = np.array([1, 2, 3])
        vegan = np.array([1, 1, 1], dtype=np.float)
        test_res = preston_bin(data, max(data))[0]
        assert_array_equal(test_res, vegan)

