"""
===============================
Models (:mod:`macroeco.models`)
===============================

This module contains distributions and curves (i.e., standard mathematical
functions) commonly used in analysis of ecological patterns.

Distributions
=============

All of the distributions here are subclasses of either
`~scipy.stats.rv_continuous` and `~scipy.stats.rv_discrete` found in
`scipy.stats`. Several of the distributions here are similar to or based on
existing distributions found in `scipy.stats` but are updated to allow the use
of common ecological parameterizations.

In addition to all of the methods found in `scipy.stats`, methods for fitting
distributions and curves to data and for translating common distribution
arguments into formal parameters (i.e., deriving the ``p`` of the geometric
distribution from the distribution mean) are also provided in these classes.

The following discrete distributions are available.

.. autosummary::
   :toctree: generated/

   geom
   geom_uptrunc
   nbinom
   nbinom_ztrunc
   cnbinom
   logser
   logser_uptrunc
   plnorm
   plnorm_ztrunc
   dgamma

The following continuous distributions are available.

.. autosummary::
   :toctree: generated/

   expon
   expon_uptrunc
   lognorm

Curves
======

Several common curves used in ecologial analysis are included here.

.. autosummary::
   :toctree: generated/

   power_law
   mete_sar
   mete_sar_iterative
   mete_ear
   sampling_sar
   sampling_sar_iterative
   sampling_ear

"""

from _distributions import (geom, geom_uptrunc, nbinom, nbinom_ztrunc,
                            cnbinom, logser, logser_uptrunc, plnorm,
                            plnorm_ztrunc, expon, expon_uptrunc, lognorm,
                            dgamma)

from ._curves import (power_law,
                      mete_sar, mete_ear, mete_sar_iterative,
                      mete_upscale_iterative_alt, sampling_sar,
                      sampling_sar_iterative, sampling_ear)
