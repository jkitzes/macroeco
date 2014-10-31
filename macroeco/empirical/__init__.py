"""
=====================================
Empirical (:mod:`macroeco.empirical`)
=====================================

This module contains functions used in the analysis of ecological patterns in
empirical data sets.

Patch
=====

Patch is the core class of the empirical module. It reads and validates
metadata and data table files, and patch objects are the first argument to all
of the empirical metric functions in this module.

.. autosummary::
   :toctree: generated/

   Patch

Metrics
=======

Each of these functions calculates an empirical ecological metric for a given
patch object.

.. autosummary::
   :toctree: generated/

   sad
   ssad
   sar
   comm_grid
   o_ring

Other
=====

.. autosummary::
   :toctree: generated/

   empirical_cdf

"""

from ._empirical import (Patch,
                         sad, ssad, sar, comm_grid, o_ring,
                         empirical_cdf)
