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

.. DV:
   Our public-facing distributions do not use location and scale parameters, as
   they are not common in quantitative ecology.
"""

from _distributions import (geom, geom_uptrunc, nbinom, cnbinom,
                            expon, expon_uptrunc)
