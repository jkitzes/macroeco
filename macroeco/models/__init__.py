"""
===============================
Models (:mod:`macroeco.models`)
===============================

This module contains distributions commonly used in analysis of ecological
patterns. At present, all distributions here are univariate.

Most of these distributions are subclasses of `~scipy.stats.rv_continuous` and
`~scipy.stats.rv_discrete` found in `scipy.stats`. Additionally, several of the
distribution classes here are simple wrappers for existing distributions found
in `scipy.stats` that are updated to allow the use of common ecological
parameterizations.

Continouous distributions
=========================

.. autosummary::
   :toctree: generated/

   expon
   expon_uptrunc

Discrete distributions
======================

.. autosummary::
   :toctree: generated/

   geom
   geom_uptrunc
   nbinom
   cnbinom
   logser_uptrunc
   lognorm

"""

from _distributions import (geom, geom_uptrunc, nbinom, cnbinom,
                                 logser_uptrunc, expon, expon_uptrunc,
                                 lognorm)

from ._curves import (power_law,
                      mete_sar, mete_iterative_sar,
                      mete_ear, mete_iterative_ear)
