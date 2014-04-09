"""
=====================================
Empirical (:mod:`macroeco.empirical`)
=====================================

This module contains functions used in the empirical analysis of
macroecological patterns.

Patch
=====

Patch is a class.

.. autosummary::
   :toctree: generated/

   Patch

Metrics
=======

.. autosummary::
   :toctree: generated/

   sad
   ssad
   sar
   ear
   comm_grid

"""

from .empirical import (Patch,
                        sad, ssad, sar, comm_grid,
                        empirical_cdf)
