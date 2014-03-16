"""
==============================================
Distributions (:mod:`macroeco.distributions2`)
==============================================

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

from __future__ import division

from decimal import Decimal
import numpy as np
import numpy.random as nprand

from scipy.stats.distributions import rv_discrete, rv_continuous

class rv_continuous_meco(rv_continuous):
    """
    A modified generic continuous random variable class meant for subclassing.

    This class inherits from the `rv_continuous` class of `scipy.stats` and 
    contains all of its functionality. See the docstring of `rv_continuous` for 
    information on usage and subclassing. In addition, this class adds one new 
    methods.

    Methods
    -------
    translate_args
        takes user-friendly params as input and returns shape params

    fit2
        calls method `fit` with fixed loc=0 and scale=1 (defaults)

    """

    def translate_args(self, *args):
        """
        Translates user-friendly arguments into shape parameters

        See distribution docstring for description of user arguments and shape 
        parameters.

        Parameters
        ----------
        uargs : floats
            User argument(s), usually easily measured and specified

        Returns
        -------
        tuple of floats
            Shape parameter(s) of distribution

        Notes
        -----
        """

        raise NotImplementedError, ("translate_args method not implemented "
                                    "for this distribution")


    def fit2(self, *args):
        """
        Return MLEs for shape parameters from data.

        Parameters
        ----------
        data : array_like
            Data to use in calculating the MLEs.
        args : floats
            Starting value(s) for shape parameters. Some may be held constant 
            (see Notes).

        Returns
        -------
        tuple of floats
            MLEs for shape parameters

        Notes
        -----
        """

        return self.fit(*args, floc=0, fscale=1)[:-2]


