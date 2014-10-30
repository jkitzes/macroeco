from __future__ import division

import numpy as np
import pandas as pd
from scipy import optimize
from mpmath import lerchphi

from ..misc import inherit_docstring_from
import _distributions as dist

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
        x = np.array(x)
        return self._vals(x, *args, **kwargs)

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
            y_pred = self.vals(x, *params)
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


class mete_sar_gen(curve):
    """
    A SAR/EAR predicted by the Maximum Entropy Theory of Ecology

    The METE SAR and EAR may be used either for downscaling, when values of A
    are less than A0, or upscaling, when values of A are greater than A0.
    Downscaling creates the traditional SAR known to ecologists, while
    upscaling is useful for estimating large-scale species richness from small-
    scale plot data.

    A keyword argument iterative is available (default is False). If True, the
    SAR is calculated at successive A values, with the result at each value of
    A used as the base values of S and N for the subsequent calculation. The
    iterative form was used in Harte et al [#]_, although note that the
    implementation here uses a different internal equation.

    Methods
    -------
    vals(x, S0, N0, iterative=False)
        Calculate SAR given starting values and two models. See notes.

    Parameters
    ----------
    x : iterable
        Areas at which to calculate SAR (first element is A0)
    S0 : float
        Species richness at A0
    N0 : float
        Community abundance at A0
    iterative : bool (opt)
        If true, SAR calculation for subplots are based on variables for next
        larger area instead of initial plot variables. Default False.
    array_size : int (opt)
        Maximum size of array for SAD pmf's. If N0 is greater than this value,
        calculation proceeds using array_size increments until N0 is reached.
    approx : bool (opt)
        Use non-truncated logseries and geometric distributions. Default False.

    References
    ----------
    .. [#]
       Harte, J., Smith, A. B., & Storch, D. (2009). Biodiversity scales from
       plots to biomes with a universal species-area curve. Ecology Letters,
       12(8), 789-797.

    """

    def __init__(self, name=None, parameters=None, iterative=False, ear=False):
        """
        Provides extra iterative attribute.
        """
        if iterative and ear:
            raise ValueError, "Iterative EAR calculation is not possible"

        self.name = name
        self.parameters = parameters
        self.n_parameters = len(parameters.split(','))
        self.iterative = iterative
        self.ear = ear

    def _vals(self, x, S0, N0, array_size=1e6, approx=False, alt_up=False):
        # x is area, y is S

        A0 = x[0]
        y = [S0]

        for A in x[1:]:
            a = A/A0

            if a == 1:
                S1, N1 = S0, N0
            elif a < 1:
                S1, N1 = self._downscale_step(a, S0, N0, array_size, approx)
            else:
                S1, N1 = self._upscale_step(a, S0, N0, array_size, approx)

            y.append(S1)

            if self.iterative:
                S0, N0, A0 = S1, N1, A

        return np.array(y)

    def _downscale_step(self, a, S0, N0, array_size, approx):
        lower = 1
        upper = array_size + 1
        S = 0
        print S0,N0

        if S0 < 1 or np.isnan(S0):  # Give up if S0 too small
            return np.nan, N0*a

        while lower < N0:
            print lower

            if upper > N0:
                upper = N0 + 1

            n0 = np.arange(lower, upper)

            if approx:
                sad_p = dist.logser.translate_args(N0/S0)
                sad = dist.logser.pmf(n0, sad_p)
            else:
                sad_p, _ = dist.logser_uptrunc.translate_args(N0/S0, N0)
                sad = dist.logser_uptrunc.pmf(n0, sad_p, N0)

            if np.isclose(a, 0.5):
                ssad_p = 1 / (n0 + 1)
            else:
                if approx:
                    ssad_p = dist.geom.translate_args(a*n0)
                else:
                    ssad_p, _ = dist.geom_uptrunc.translate_args(a*n0, N0)

            if self.ear:
                if approx:
                    ssad = dist.geom.pmf(n0, ssad_p)
                else:
                    ssad = dist.geom_uptrunc.pmf(n0, ssad_p, N0)
                S += S0 * np.sum(ssad * sad)
            else:
                if approx:
                    ssad = dist.geom.pmf(0, ssad_p)
                else:
                    ssad = dist.geom_uptrunc.pmf(0, ssad_p, N0)
                S += S0 * np.sum((1 - ssad) * sad)

            lower += array_size
            upper += array_size

        return S, N0*a

    def _upscale_step(self, a, S0, N0, array_size, approx):

        N1 = N0*a
        print a

        def eq(S1, N1, a, S0, array_size, approx):
            return S0-self._downscale_step(1/a, S1, N1, array_size, approx)[0]

        return optimize.brentq(eq,S0,S0*a,args=(N1,a,S0,array_size,approx)), N1

    def fit_lsq(self, df):
        """
        Parameterize generic SAR curve from empirical data set

        Parameters
        ----------
        df : DataFrame
            Result data frame from empirical SAR analysis

        Notes
        -----
        Simply returns S0 and N0 from empirical SAR output, which are two fixed
        parameters of METE SAR and EAR. The first row of the empirical
        dataframe corresponds to area A0. Name ``fit_lsq`` is retained for
        consistency with other curves.

        """
        # Just return S0 and N0 at largest scale, which is first row of df
        return df['n_spp'].values[0], df['n_individs'].values[0]


mete_sar = mete_sar_gen(name='mete_sar', parameters='S0,N0')
mete_sar_iterative = mete_sar_gen(name='mete_iterative_sar',
                                  parameters='S0,N0', iterative=True)

mete_ear = mete_sar_gen(name='mete_ear', parameters='S0,N0', ear=True)

def mete_upscale_iterative_alt(S, N, doublings):
    """
    This function is used to upscale from the anchor area.

    Parameters
    ----------
    S : int or float
        Number of species at anchor scale
    N : int or float
        Number of individuals at anchor scale
    doublings : int
        Number of doublings of A. Result vector will be length doublings + 1.

    Returns
    -------
    result : ndarray
        1D array of number of species at each doubling

    """

    # Arrays to store N and S at all doublings
    n_arr = np.empty(doublings+1)
    s_arr = np.empty(doublings+1)

    # Loop through all scales
    for i in xrange(doublings+1):

        # If this is first step (doubling 0), N and S are initial values
        if i == 0:
            n_arr[i] = N
            s_arr[i] = S

        # If not first step
        else:

            # Get previous S
            SA = s_arr[i-1]

            # N is double previous N
            n_arr[i] = 2 * n_arr[i-1]
            N2A = n_arr[i]

            # Eq 8 from Harte 2009, setup to return S2A given input of x
            # x is exp(-lam_phi, 2A)
            def S2A_calc(x, SA, N2A):
                return ((SA +
                        N2A *
                        (1-x)/(x-x**(N2A+1)) *
                        (1 - (x**N2A)/(N2A+1))) /
                        x**-1)

            # Eq 9 from Harte 2009, setup to equal to zero, used to solve x
            # Note that two summations are replaced by known formulas for sum
            # of geometric and logarithmic series.
            # Note "offset" of 1e-23, which is needed because f(a) and f(b) do
            # not have the same sign in solver otherwise. This introduces no
            # more than a 1e-23 error in the calculation of x, which should not
            # cause a significant problem.
            def x_calc(x, SA, N2A):
                return (S2A_calc(x,SA,N2A) /
                        N2A *
                        x*(x**N2A-1)/(x-1) -
                        (x**N2A * (-lerchphi(x,1,N2A+1))-np.log(1-x)) ) - 1e-23

            # Solve for x
            x = (optimize.brentq(x_calc, 1e-24, 1-1e-16, args=(SA,N2A),
                    xtol=1e-16, maxiter=1000, disp=True) + 1e-23)

            # Given x, calculate S2A
            s_arr[i] = S2A_calc(x,SA,N2A)

    return s_arr
