"""
=================================
Compare (:mod:`macroeco.compare`)
=================================

This module contains functions that compare the goodness of fit of a
distribution/curve to data or the fit of two distributions/curves to each
other.

.. autosummary::
   :toctree: generated/

   nll
   lrt
   AIC
   AIC_compare
   sum_of_squares
   r_squared
   preston_bin

"""

from ._compare import (nll, lrt, AIC, AIC_compare,
                       sum_of_squares, r_squared,
                       preston_bin, pueyo_bins)
