from __future__ import division

from numpy.testing import (TestCase, assert_equal, assert_array_equal,
                           assert_almost_equal, assert_array_almost_equal,
                           assert_allclose, assert_, assert_raises)

import numpy as np
from decimal import Decimal
from macroeco.models import *
import scipy as sp
import scipy.stats as stats
import matplotlib.pyplot as plt

class SAMPLING_SAR(TestCase):

    def test_reversible(self):
        S0, N0 = 20, 1e3
        As = np.array([100, 50, 10])
        Ns = N0 * As / As[0]

        Ss = sampling_sar.vals(As, S0, N0, sad_k=0.5, ssad_k=0.5, approx=True)

        # Start with each smaller base and go directly up to A0
        for A, S, N in zip(As[1:], Ss[1:], Ns[1:]):
            assert_almost_equal(S0,
                sampling_sar.vals([A, As[0]], S, N, sad_k=0.5, ssad_k=0.5,
                                        approx=True)[1], decimal=5)

class SAMPLING_iterative_SAR(TestCase):

    def test_reversible(self):
        S0, N0 = 20, 1e3
        As = 1 / 2**np.arange(0, 4)

        Ns = N0 * As / As[0]

        Ss = sampling_sar_iterative.vals(As, S0, N0, sad_k=0.5, ssad_k=0.5,
                                approx=True)

        assert_array_almost_equal(Ss[::-1],
            sampling_sar_iterative.vals(As[::-1], Ss[-1], Ns[-1], sad_k=0.5,
            ssad_k=0.5, approx=True), decimal=1)

class SAMPLING_EAR(TestCase):

    def test_no_upscale(self):
        # Not allowing for EAR upscaling right now

        assert_raises(NotImplementedError, sampling_ear.vals, [1, 2], 50, 1000, 1, 1)

    def test_make_green_plot(self):

        # Testing whether the sampling EAR shows the same qualitative results
        # as the plot given in Green and Ostling (2003): Endemics-area
        # relationship: the influence of species dominance and spatial
        # aggregation

        ## Uncomment to make plot: Takes a while (~ 2 minutes) to generate plot###

        # # Data from Pasoh
        # S0 = 814
        # N0 = 335356

        # sad_ks = [0, 0.5, 1]
        # areas = 1 / 2**np.arange(0, 20)

        # fig, axes = plt.subplots(2, 1, figsize=(5, 8))
        # axes = axes.ravel()

        # # ssad_k = 10 is approximately random placement
        # for sad_k in sad_ks:

        #     temp_sar = sampling_sar.vals(areas, S0, N0, sad_k=sad_k, ssad_k=10, approx=True)
        #     temp_ear = sampling_ear.vals(areas, S0, N0, sad_k=sad_k, ssad_k=10, approx=True)

        #     axes[0].plot(np.log(areas), temp_sar, '-o')
        #     axes[0].plot(np.log(areas), temp_ear, '-s')

        #     axes[1].plot(np.log(areas), np.log(temp_sar), '-o')
        #     axes[1].plot(np.log(areas), np.log(temp_ear), '-s')

        # axes[0].set_ylabel("Species/Endemics Richness")
        # axes[1].set_ylabel("log(Species/Endemics Richness)")

        # axes[0].set_xlabel("log(Area)")
        # axes[1].set_xlabel("log(Area)")

        pass

class METE_SAR(TestCase):

    def test_reversible(self):
        S0, N0 = 100, 1e6
        As = np.array([100,50,10])
        Ns = N0 * As / As[0]

        Ss = mete_sar.vals(As, S0, N0, approx=True)

        # Start with each smaller base and go directly up to A0
        for A, S, N in zip(As[1:], Ss[1:], Ns[1:]):
            assert_almost_equal(S0,
                mete_sar.vals([A, As[0]], S, N, approx=True)[1], decimal=5)

    def test_reproduce_mcglinn_figure(self):

        # Remake figure from McGlinn et al. 2013 using iterative and
        # non-iterative METE SAR

        ## UNCOMMENT TO MAKE PLOT ##

        # data = [(205096, 301), (7622.5, 174.5), (4326, 138.5), (32320, 124),
        #         (3394, 48), (5469, 40), (8892, 39), (3383, 37), (3879, 30),
        #         (758, 19), (669, 17), (2139, 41), (2584, 36), (5885, 31),
        #         (37182, 24), (7625, 7)]

        # plot_areas = np.array([50, 2, 2, 12.5, 1.7113, 2, 6.5536, 1.44, 1.96, 0.5041,
        #             0.5, 0.845, 1, 4.5, 0.0064, 4]) * 10000

        # # Make plots
        # fig, axes = plt.subplots(4, 4, figsize=(8, 8))
        # axes = axes.ravel()

        # for ax in axes:
        #     ax.tick_params(labelsize=5)

        # for i, (N, S) in enumerate(data):

        #     areas = 1 / 2**np.arange(0, 14)

        #     try:
        #         non_iter_sar = mete_sar.vals(areas, S, N, approx=True)
        #         iter_sar = mete_sar_iterative.vals(areas, S, N, approx=True)
        #     except:
        #         areas = 1 / 2**np.arange(0, 9)
        #         non_iter_sar = mete_sar.vals(areas, S, N, approx=True)
        #         iter_sar = mete_sar_iterative.vals(areas, S, N, approx=True)

        #     axes[i].plot(np.log2(areas*plot_areas[i]), np.log2(non_iter_sar), color="blue",
        #                             label="Non-iterative")
        #     axes[i].plot(np.log2(areas*plot_areas[i]), np.log2(iter_sar), color="red",
        #                             label="Iterative")

        #     if axes[i].is_first_col():
        #         axes[i].set_ylabel("log2(Species Richness", size=9)
        #     if axes[i].is_last_row():
        #         axes[i].set_xlabel("log2(Area)", size=9)

        pass

    def test_vals_up(self):
        pass

class METE_iterative_SAR(TestCase):

    def test_reversible(self):
        S0, N0 = 100, 1e6
        As = 1 / 2**np.arange(0, 4)

        Ns = N0 * As / As[0]

        Ss = mete_sar_iterative.vals(As, S0, N0, approx=True)

        assert_array_almost_equal(Ss[::-1],
            mete_sar_iterative.vals(As[::-1], Ss[-1], Ns[-1], approx=True),
            decimal=1)

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

class METE_EAR(TestCase):

    def test_no_upscale(self):
        # Not allowing for EAR upscaling right now

        assert_raises(NotImplementedError, mete_ear.vals, [1, 2], 50, 1000)

    def test_mete_ear_with_data(self):

        BCI_S = 283
        BCI_N = 208310

        SERP_S = 28
        SERP_N = 60346

        areas_serp = 1 / 2**np.arange(0, 9) * 0.0016 * 10000
        serp_ear = mete_ear.vals(areas_serp, SERP_S, SERP_N, approx=True)

        areas_bci = 1 / 2**np.arange(0, 14) * 50 * 10000
        bci_ear = mete_ear.vals(areas_bci, BCI_S, BCI_N, approx=True)

        ## Uncomment for plots and see Harte (2011) pg 193 for comparison ##

        # fig, axes = plt.subplots(2, 1, figsize=(5, 8))
        # axes = axes.ravel()
        # for ax in axes:
        #     ax.set_xlabel("log(area)")
        #     ax.set_ylabel("log(endemics)")

        # axes[0].plot(np.log(areas_serp), np.log(serp_ear), 'd')
        # axes[1].plot(np.log(areas_bci), np.log(bci_ear), 'd')

        # axes[0].set_ylim((-6, 4))
        # axes[1].set_ylim((-6, 8))




