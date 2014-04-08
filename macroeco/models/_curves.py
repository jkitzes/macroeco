from __future__ import division

import numpy as np
import pandas as pd
from scipy import optimize

from ..misc import inherit_docstring_from

_doc_methods = \
"""Methods
    -------
    vals(x, parameters)
        Dependent variable y given independent variable x and curve parameters
    fit_lsq(x, y_obs, params_start=None)
        Least squares fit of parameters given data"""

_doc_parameters = \
"""Parameters
    ----------
    x : iterable
        Independent variable
    y_obs : iterable
        Dependent variable (values observed at x)
    params_start : iterable
        Optional start values for all parameters. Default 1."""


class curve(object):
    """
    Generic function class meant for subclassing
    """

    def __init__(self, name=None, parameters=None):
        """
        Distribution parameters may be given here or to individual methods

        """
        self.name = name
        self.parameters = parameters
        self.n_parameters = len(parameters.split(','))

    def __call__(self, *args, **kwargs):
        raise ValueError, "Choose either the vals or fit_lsq methods"

    def vals(self, x, *args, **kwargs):
        """
        [Docstring]

        """
        self.vals_kwargs = kwargs
        x = np.array(x)
        y = self._vals(x, *args)
        return pd.DataFrame({'x': x, 'y': y})

    def _vals(self, x, *args):
        """
        Return y given x and parameters
        """
        raise NotImplementedError, ("vals not implemented for %s" % self.name)

    def fit_lsq(self, x, y_obs, params_start=None):
        """
        Fit curve by method of least squares.

        Parameters
        ----------
        x : iterable
            Independent variable
        y_obs : iterable
            Dependent variable (values observed at x)
        params_start : iterable
            Optional start values for all parameters. Default 1.

        Returns
        -------
        array
            Best fit values of parameters

        Notes
        -----
        If least squares fit does not converge, ValueError is raised with
        convergence message.

        """

        # Set up variables
        x = np.array(x)
        y_obs = np.array(y_obs)
        if not params_start:
            params_start = np.ones(self.n_parameters)

        # Error checking
        if len(x) != len(y_obs):
            raise ValueError, "x and y_obs must be the same length"
        if len(params_start) != self.n_parameters:
            raise ValueError, "Incorrect number of values in params_start"

        # Calculate fit
        def residuals(params, x, y_obs):
            y_pred = self.vals(x, *params)['y']
            return y_obs - y_pred

        params_fit, _, _, msg, ier = optimize.leastsq(residuals, params_start,
                                             args=(x, y_obs), full_output=True)

        # Check for convergence
        if ier > 4:
            raise ValueError, ("Least squares fit did not converge with "
                               "message %s" % msg)

        return tuple(params_fit)


class power_law_gen(curve):
    """
    A power-law function

    .. math::

       y = c x^z

    or equivalently

    .. math::

       \log(y) = \log(c) + z \log(x)

    Stemming from the log form, ``c`` is often known as the intercept and ``z``
    as the slope of the power law.

    {0}

    {1}
    c, z
        Parameters: Log-log slope and intercept

    """

    def _vals(self, x, c, z):
        return c * x**z

power_law = power_law_gen(name='power_law', parameters='c,z')
power_law.__doc__ = power_law.__doc__.format(_doc_methods, _doc_parameters)


class gen_sar_gen(curve):
    """
    INCOMPLETE NEEDS CONTINUED WORK

    A generic SAR based on a combination of an SAD and SSAD

    .. math::

       y = c x^z

    The generic SAR may be used either for downscaling, when values of A are
    less than A0, or upscaling, when values of A are greater than A0.
    Downscaling creates the traditional SAR known to ecologists, while
    wpscaling is particularly useful for estimating large-scale species
    richness from small-scale plot data.

    A keyword argument iterative is available for the generic SAR (default is
    False). If True, the SAR is calculated at successive A values, with the
    result at each value of A used as the base values of S0 and N0 for the
    subsequent calculation. The iterative SAR form is a generalization of the
    universal SAR proposed by Harte et al [#]_.

    Methods
    -------
    vals(S0, N0, A, SAD_model, SSAD_model)
        Calculate SAR given starting values and two models. See notes.

    Parameters
    ----------
    S0 : float
        Species richness at A0
    N0 : float
        Community abundance at A0
    A : iterable
        Areas at which to calculate SAR (first element is A0)
    SAD_model : object
        Frozen distribution from macroeco.models
    SSAD_model : object
        Frozen distribution from macroeco.models
    iterative : bool
        If true, SAR calculation for subplots are based on variables for next
        larger area instead of initial plot variables. Default False.

    References
    ----------
    .. [#]
       Harte, J., Smith, A. B., & Storch, D. (2009). Biodiversity scales from
       plots to biomes with a universal species-area curve. Ecology Letters,
       12(8), 789-797.

    """

    def _vals(self, x, S0, N0, iterative=False):
        # x is area, y is S

        A0 = x[0]
        y = [S0]

        for A in x[1:]:
            S1, N1 = self._single_step(S0, N0, A/A0)
            y.append(S1)
            if iterative:
                S0, N0, A0 = S1, N1, A

        return np.array(y)

    def _single_step(self, S0, N0, a):
        # if a < 1, solve, if a > 1, guess and check
        if a == 1:
            S1 = S0
            N1 = N0
        elif a < 1:  # "Normal" downscale
            S1 = S0
            N1 = N0
        else:  # Upscale solver
            S1 = S0
            N1 = N0

        return S1, N1

gen_sar = gen_sar_gen(name='gen_sar', parameters='S0,N0')











