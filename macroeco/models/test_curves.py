from __future__ import division

from numpy.testing import (TestCase, assert_equal, assert_array_equal,
                           assert_almost_equal, assert_array_almost_equal,
                           assert_allclose, assert_, assert_raises)

import numpy as np
from decimal import Decimal
from macroeco.models import *
import scipy as sp
import scipy.stats as stats


class METE_SAR(TestCase):

    def test_reversible(self):
        S0, N0 = 100, 1e6
        As = np.array([100,50,10])
        Ns = N0 * As / As[0]

        Ss = mete_sar.vals(As, 100, 1e6, approx=True)

        # Start with each smaller base and go directly up to A0
        for A, S, N in zip(As[1:], Ss[1:], Ns[1:]):
            assert_almost_equal(S0,
                mete_sar.vals([A, As[0]], S, N, approx=True)[1])

    def test_vals_down(self):
        pass

    def test_vals_up(self):
        pass

class METE_iterative_SAR(TestCase):

    def test_reversible(self):
        S0, N0 = 100, 1e6
        As = np.array([100,50,10])
        Ns = N0 * As / As[0]

        Ss = mete_sar_iterative.vals(As, 100, 1e6, approx=True)

        assert_array_almost_equal(Ss[::-1],
            mete_sar_iterative.vals(As[::-1], Ss[-1], Ns[-1], approx=True))

    def test_vals_down(self):
        pass

    def test_vals_up(self):
        # ACARI results from Bassett upscaling paper, see SI
        # Note that different approximations are used here and in that analysis
        S0, N0 = 86.6, 2015
        As = [0.04, 0.08, 0.16, 0.32, 0.64, 1.28, 2.56]

        Ss = mete_sar_iterative.vals(As, S0, N0, approx=True)

        assert_array_almost_equal(Ss,
            [86.6, 106.0327113, 127.1223631, 149.7292838,
             173.7360065, 199.0452844, 225.5766732])
