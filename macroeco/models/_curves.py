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
            params_start = np.ones(len(self.parameters))

        # Error checking
        if len(x) != len(y):
            raise ValueError, "x and y_obs must be the same length"
        if len(params) != self.n_parameters:
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

        return params_fit


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

class mete_sar_gen(curve):
    """
    The SAR predicted by the Maximum Entropy Theory of Ecology

    .. math::

       y = c x^z

    or equivalently

    .. math::

       \log(y) = \log(c) + z \log(x)

    {0}

    {1}
    S0, N0
        Parameters: Initial species richness and community abundance at largest
        scale
    iterative : bool
        If true, SAR calculation for subplots are based on variables for next
        larger area instead of initial plot variables. Default False.

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
            pass
        else:  # Upscale solver
            pass

        return S1, N1

mete_sar = mete_sar_gen(name='mete_sar', parameters='S0,N0')
mete_sar.__doc__ = mete_sar.__doc__.format(_doc_methods, _doc_parameters)












