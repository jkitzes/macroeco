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
        x = np.atleast_1d(x)
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
        x = np.atleast_1d(x)
        y_obs = np.atleast_1d(y_obs)
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

    Examples
    --------

    >>> # Specify a classic power law with z = 0.25
    >>> import macroeco.models as md

    >>> areas = [1, 0.5, 0.25, 0.125]
    >>> c = 20  # Number of species at the base scale
    >>> z = 0.25 # Slope of the power law

    >>> # Get the species richness predictions of the power law
    >>> res = md.power_law.vals(areas, c, z)
    >>> res
    array([ 20.        ,  16.81792831,  14.14213562,  11.89207115])

    >>> # Fit the power law using least squares
    >>> md.power_law.fit_lsq(areas, res)
    (20.0, 0.25000000000000006)

    """

    def _vals(self, x, c, z):
        return c * x**z

power_law = power_law_gen(name='power_law', parameters='c,z')
power_law.__doc__ = power_law.__doc__.format(_doc_methods, _doc_parameters)

class sampling_sar_gen(curve):
    """
    A general sampling SAR/EAR

    As described in Wilber et al. 2015 [#]_, a sampling
    SAR is defined by a species abundance distribution (SAD, :math:`\phi`) and a
    species-level spatial abundance distribution (SSAD, :math:`\Pi`)

    .. math::

        S(A) = S_0 \sum_{n_0=1}^{N_0} \phi(n_0 | \Theta_{\phi}) [1 - \Pi(0 | A / A_0, n_0, \Theta_{\Pi})]

    where :math:`\Theta_{\phi}` and :math:`\Theta_{\Pi}` defines the parameters
    of the SAD and SSAD respectively, :math:`S_0` is the total number of
    species at the base area :math:`A_0`, :math:`N_0` is the total number of
    individual at the base area, and `A` is the area
    at which to calculate species richness.

    A flexible choice for the SAD is a zero-truncated negative binomial with
    parameters :math:`\mu` and :math:`k_{SAD}`. Fisher Logseries
    (:math:`k_{SAD} = 0`), Canonical Lognormal (:math:`0 < k_{SSAD} < 1`), and
    Broken Stick (:math:`k_{SAD} = 1`) SADs can all be exactly or approximately
    represented by this distribution.  A flexible choice for the SSAD is the
    condition (finite) negative binomial distribution (FNBD) with parameters
    :math:`\mu` and :math:`k_{SSAD}`. When :math:`k_{SSAD}` is large a binomial
    distribution is obtained and when :math:`k_{SSAD} = 1` and truncated
    geometric distribution is obtained. See Wilber et al. 2015 for more
    information.  We implement these two distributions in the general sampling
    SAR.

    The general sampling SAR and EAR may be used either for downscaling, when
    values of :math:`A` are less than :math:`A_0`, or upscaling, when values of
    :math:`A` are greater than :math:`A0`. Downscaling creates the traditional
    SAR known to ecologists, while upscaling is useful for estimating large-
    scale species richness from small- scale plot data.

    The parameters required for the sampling SAR are species richness at the
    base scale `A_0` (:math:`S_0`), total community abundance at the base scale
    (:math:`N_0`), the aggregation parameter of the SAD (:math:`k_{SAD}`), and
    the aggregation parameter of the SSAD (:math:`k_{SSAD}`). See examples
    below.

    If the standard SAR is chosen (`sampling_sar`), the SAR is calculated by
    solving the above equation for any given value :math:`A` greater than or
    less than :math:`A_0`.

    If the iterative SAR is chosen (`sampling_sar_iterative`), the SAR is
    calculated by successively halving (if downscaling) or successively
    doubling (if upscaling) the base area, calculating the values S and N at
    this new scale, and then setting these calculated values of S and N as the
    base :math:`S_0` and :math:`N_0` for subsequent calculations. This
    iterative form was used in Harte et al [#]_, although note that the
    implementation here uses a different internal equation. Note that the
    the iterative form of the SAR can only calculate species
    richness at values of `A` that are doublings or halvings of the `A_0`.  Any
    value of `A` can be passed to `sampling_sar_iterative`, but it will return
    the species richness at the closest iterative halving that is less than or equal to
    the given :math:`A` or the closest doubling that is greater than or equal to
    the given :math:`A`. See examples below.

    If the endemics area relationship (`sampling_ear`) is choosen (as
    given in Harte (2011) pg. 46), the number of endemics are calculated at any
    given :math:`A \leq A_0` where :math:`a < 1`. This method requires the
    same parameters as the SAR.

    Methods
    -------
    vals(x, S0, N0, sad_k, ssad_k, approx=True)
        Calculate SAR given starting values and two aggregation parameters.
        See notes.

    Parameters
    ----------
    x : iterable
        Areas at which to calculate SAR or EAR. The first element must be A0,
        the base area.
    S0 : float
        Species richness at A0
    N0 : float
        Community abundance at A0
    sad_k : float
        Aggregation parameter of the SAD. Between [0, infinity)
    ssad_k : float
        Aggregation parameter of the SSAD. Between (0, infinity).
    approx : bool (opt)
        Approximate the truncated logseries. Default True. The approximation is
        much faster and not very different than the exact answer for most
        cases.

    Examples
    --------

    >>> # Use a standard SAR with a Broken Stick SAD and an truncated geometric
    >>> # binomial SSAD. sad_k = 1, ssad_k = 1
    >>> import macroeco.models as md

    >>> # Specify the areas at which to calculate species richness
    >>> areas = [1, 0.5, 0.25, 0.125] # 1 is the base area A0

    >>> # Set community abundance and species richness at A0
    >>> N0, S0 = (10000, 50)

    >>> # Get the standard SAR
    >>> md.sampling_sar.vals(areas, S0, N0, sad_k=1, ssad_k=1, approx=True)
    array([ 50.        ,  48.91333114,  47.3371818 ,  45.00844875])

    >>> # Get the iterative SAR
    >>> md.sampling_sar_iterative.vals(areas, S0, N0, sad_k=1, ssad_k=1, approx=True)
    array([ 50.        ,  48.91333114,  47.13849657,  44.37798735])

    >>> # Get the EAR
    >>> md.sampling_ear.vals(areas, S0, N0, sad_k=1, ssad_k=1, approx=True)
    array([  5.00000000e+01,   1.08666886e+00,   1.23816570e-01, 4.15836438e-02])

    >>> # Upscaling species richness
    >>> areas = [5, 10, 11, 13]
    >>> md.sampling_sar.vals(areas, S0, N0, sad_k=0, ssad_k=1, approx=True)
    array([ 50.        ,  57.19380999,  58.10403751,  59.66376591])

    >>> # Iterative SAR with doubled areas
    >>> areas_up = [1, 2, 4, 8]
    >>> md.sampling_sar_iterative.vals(areas_up, S0, N0, sad_k=0, ssad_k=1, approx=True)
    array([ 50.        ,  57.19380999,  64.72987319,  72.58968402])

    >>> # Iterative SAR with not quite doubled areas
    >>> areas_up = [1, 2.1, 4.1, 8.1]
    >>> md.sampling_sar_iterative.vals(areas_up, S0, N0, sad_k=0, ssad_k=1, approx=True)
    array([ 50.        ,  64.72987319,  72.58968402,  80.75754743])

    >>> # Notice that the iterative method rounds up if the areas are not
    >>> # exact doublings.  This is equivalent to the following
    >>> md.sampling_sar_iterative.vals([1, 4, 8, 16], S0, N0, sad_k=0, ssad_k=1, approx=True)
    array([ 50.        ,  64.72987319,  72.58968402,  80.75754743])

    References
    ----------
    .. [#]
       Wilber, M., Kitzes, J., and Harte, J. (2015). Scale collapse and the
       emergence of the power law species-area relationship. Global Ecology and
       Biogeography. 24(8), 883-895

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

    def _downscale_direct(self, a, S, N, sad_k, ssad_k, approx):

        # Parameterize the appropriate SAD
        if sad_k == 0:
            if approx:
                sad = _logser(N=N, S=S)
            else:
                sad = _logser_uptrunc(N=N, S=S)
        else:
            sad = _nbinom_ztrunc(N=N, S=S, k=sad_k)

        n_vals = np.arange(1, N + 1)
        ssad = _cnbinom(N=n_vals, a=a, k=ssad_k)

        # Sampling SAR formula

        if self.ear:
            down_S = S * np.sum(sad.pmf(n_vals) * ssad.pmf(n_vals))
        else:
            down_S = S * np.sum(sad.pmf(n_vals) * (1 - ssad.pmf(0)))

        return down_S

    def _upscale_direct(self, a, S, N, sad_k, ssad_k, approx):

        # Don't bother trying to upscale the EAR
        if self.ear:
            raise NotImplementedError("Upscaling EAR not implemented")

        up_N = np.round(N * a, decimals=0)
        n_vals = np.arange(1, up_N + 1)

        if sad_k == 0:
            if approx:
                sad = _logser(N=up_N)
            else:
                sad = _logser_uptrunc(N=up_N)
        else:
            sad = _nbinom_ztrunc(N=up_N, k=sad_k)

        ssad = _cnbinom(N=n_vals, a=1 / a, k=ssad_k)

        def up_fxn(up_S):
            # Find the zero of this function to get upscaled richness

            sad.S = up_S
            x1 = up_S * np.sum(sad.pmf(n_vals) * (1 - ssad.pmf(0))) - S
            return x1

        S_calc = optimize.brentq(up_fxn, S, a * S, xtol=1e-5)

        return S_calc

    def _downscale_iterative(self, a, S, N, sad_k, ssad_k, approx, delta=2):
        # delta = 2 says that we move along the iterative SAR by
        # halving and doubling

        S_chain = [S]

        a_down = 1
        N_down = N
        S_down = S

        # Keep iterating while area is bigger than desired area
        while a_down > a:

            N_down = np.round(N_down, decimals=0)
            a_down = a_down * (1 / delta)

            S_down = self._downscale_direct(1 / delta, S_down, N_down, sad_k,
                                            ssad_k, approx)
            S_chain.append(S_down)

            N_down = N_down * (1 / delta)

        return np.array(S_chain)[-1]

    def _upscale_iterative(self, a, S, N, sad_k, ssad_k, approx, delta=2):
        # delta = 2 says that we move along the iterative SAR by
        # halving and doubling

        a_up = 1
        N_up = N
        S_up = S

        S_chain = [S]

        while a_up < a:

            N_up = np.round(N_up, decimals=0)
            a_up = a_up * delta

            S_up = self._upscale_direct(delta, S_up, N_up, sad_k, ssad_k,
                                            approx)
            S_chain.append(S_up)

            # Reset N
            N_up = N_up * (delta)

        # Only return the last value of the chain
        return np.array(S_chain)[-1]

    def _vals(self, x, S0, N0, sad_k, ssad_k, approx=True):
        # x is area

        if self.iterative:
            downscale = self._downscale_iterative
            upscale = self._upscale_iterative
        else:
            downscale = self._downscale_direct
            upscale = self._upscale_direct

        x = np.atleast_1d(x)
        areas = x / x[0]
        S_vals = [S0]

        for a in areas[1:]:

            if a == 1:
                S1 = S0
            elif a < 1:
                S1 = downscale(a, S0, N0, sad_k, ssad_k, approx)
            else:
                S1 = upscale(a, S0, N0, sad_k, ssad_k, approx)

            S_vals.append(S1)

        return np.array(S_vals)

    @inherit_docstring_from(curve)
    def fit_lsq(self, df):

        raise NotImplementedError("Not implemented for sampling sar")


sampling_sar = sampling_sar_gen(name='sampling_sar',
                                parameters='S0,N0,sad_k,ssad_k',
                                iterative=False)

sampling_ear = sampling_sar_gen(name="sampling_ear",
                                parameters='S0,N0,sad_k,ssad_k',
                                iterative=False, ear=True)

sampling_sar_iterative = sampling_sar_gen(name='sampling_sar',
                                parameters='S0,N0,sad_k,ssad_k',
                                iterative=True)


class mete_sar_gen(sampling_sar_gen):
    """
    A SAR/EAR predicted by the Maximum Entropy Theory of Ecology

    The METE SAR/EAR is a special case of the general sampling sar when
    :math:`k_{SAD} = 0` and :math:`k_{SSAD} = 1` described in
    Harte et al. (2009) [#]_.  See the documentation for `sampling_sar` for
    more information.

    The METE SAR and EAR may be used either for downscaling, when values of A
    are less than A0, or upscaling, when values of A are greater than A0.
    Downscaling creates the traditional SAR known to ecologists, while
    upscaling is useful for estimating large-scale species richness from small-
    scale plot data.

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
    approx : bool (opt)
        Approximate the truncated logseries. Default True. The approximation is
        much faster and not very different than the exact answer for most
        cases.

    Examples
    --------
    >>> # For a base area of A = 10, downscale species richness using the
    >>> # METE SAR
    >>> import macroeco.models as md

    >>> # Set the areas at which to downscale species richness
    >>> areas = [10, 8, 5, 3, 0.5]

    >>> # Set the species richness and abundance at the base scale
    >>> S0 = 50
    >>> N0 = 4356

    >>> # Standard METE SAR
    >>> md.mete_sar.vals(areas, S0, N0, approx=True)
    array([ 50.        ,  47.2541949 ,  42.16880014,  37.32787098,  22.94395771])

    >>> # Iterative METE SAR
    >>> md.mete_sar_iterative.vals(areas, S0, N0, approx=True)
    array([ 50.        ,  42.16880014,  42.16880014,  34.89587044,  16.97747262])

    >>> # METE EAR
    >>> md.mete_ear.vals(areas, S0, N0, approx=True)
    Out[12]: array([ 50.        ,  24.4391921 ,   7.83119986,   3.3844614 ,   0.41615551])

    References
    ----------
    .. [#]
       Harte, J., Smith, A. B., & Storch, D. (2009). Biodiversity scales from
       plots to biomes with a universal species-area curve. Ecology Letters,
       12(8), 789-797.

    """

    def _vals(self, x, S0, N0, approx=True):

        sampling_sar = sampling_sar_gen(self.name, self.parameters,
                                                self.iterative, self.ear)

        # sad_k = 0 and ssad_k = 1 for the sampling
        return sampling_sar._vals(x, S0, N0, 0, 1, approx)


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
        parameters of METE SAR and EAR. This simply returns n_spp and
        n_individs from the 1,1 division in
        the dataframe. An error will be thrown if this division is not present
        The ``fit_lsq`` is retained for consistency with other curves.

        """
        tdf = df.set_index('div')
        return tdf.ix['1,1']['n_spp'], tdf.ix['1,1']['n_individs']


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


class _logser:
    """
    Logseries defined in terms of N and S.

    See macroeco.models.logser or full docstring

    Parameters
    ----------
    N : int
        total number of individuals
    S : int
        total number of species

    """

    def __init__(self, N=None, S=None):
        self.N = N
        self.S = S

    def pmf(self, x):

        p = dist.logser.translate_args(self.N / self.S)
        return dist.logser.pmf(x, p) / dist.logser.cdf(self.N, p)

class _logser_uptrunc:
    """
    Logseries defined in terms of N and S.

    See macroeco.models.logser_uptrunc or full docstring

    Parameters
    ----------
    N : int
        total number of individuals
    S : int
        total number of species

    """

    def __init__(self, N=None, S=None, k=None):
        self.N = N
        self.S = S
        self.k = k

    def pmf(self, x):

        p = dist.logser_uptrunc.translate_args(self.N / self.S, self.N)
        return dist.logser_uptrunc.pmf(x, *p)

class _nbinom_ztrunc:
    """
    Zero-truncated negative binomial in terms of N and S

    See macroeco.models.nbinom_ztrunc for full docstring

    Parameters
    ----------
    N : int
        total number of individuals
    S : int
        total number of species
    k : float
        Aggregation parameter. k = 1 in broken stick, the most regular
        realistic SAD.  k between 0 and 1 represents the range of a canoncial
        lognormal.

    """

    def __init__(self, N=None, S=None, k=None):
        self.N = N
        self.S = S
        self.k = k

    def pmf(self, x):

        mu = self.N / self.S
        return dist.nbinom_ztrunc.pmf(x, mu, self.k) / \
                                    dist.nbinom_ztrunc.cdf(self.N, mu, self.k)

class _cnbinom:
    """
    Conditional Negative Binomial distribution defined in terms of N (total
    number of individuals) and area fraction (area relative to base area).

    See macroeco.models.cnbinom for full docstring

    Parameters
    -----------
    N : int
        Number of individuals
    a : float
        Area fraction (can be greater than 1)
    k : float
        k parameter of cnbinom
    """

    def __init__(self, N=None, a=None, k=None):
        self.N = N
        self.a = a
        self.k = k

    def pmf(self, x):

        mu = self.a * self.N

        if self.a == 0.5 and self.k == 1:
            return 1 / (1 + self.N)
        else:
            return dist.cnbinom.pmf(x, mu, self.k, self.N)
