"""
=================================
Compare (:mod:`macroeco.compare`)
=================================

This module contains functions that compare the goodness of fit of a
distribution/curve to data or the fit of two distributions/curves to each
other.

Comparison Functions
====================

.. autosummary::
   :toctree: generated/

   nll
   lrt
   AIC
   AIC_weights
   sum_of_squares
   r_squared
   bin_data

"""

from .compare import (nll, lrt, AIC, AIC_weights,
                      sum_of_squares, r_squared,
                      bin_data)
