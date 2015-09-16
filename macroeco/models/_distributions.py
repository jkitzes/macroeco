from __future__ import division
import sys

from decimal import Decimal
import numpy as np
import numpy.random as nprand
from scipy.stats.distributions import (rv_discrete, rv_continuous)

import scipy.stats as stats
import scipy.optimize as optim
import scipy.special as special
import scipy.integrate as integrate

from _distributions_docstrings import (docdict, docdict_discrete)
from ..misc import doc_sub, inherit_docstring_from


# Remove header from all methods
_docdict_allmeth = docdict['allmethods'][16:]
_docdict_discrete_allmeth = docdict_discrete['allmethods'][17:]

# **kwds in expect string followed by no space was throwing warning
_docdict_allmeth = _docdict_allmeth.replace(', **kwds','')

# Additional docstrings for custom methods
_docdict_rank_method = \
"""rank(n, %(shapes)s)
    Predicted rank abundance distribution.
"""

_docdict_extra_params = \
"""n : int
    number of values
data : array_like
    values used to fit distribution
"""

# Create docstring helpers
docdict['before_notes'] = ''.join([_docdict_rank_method,
                                   _docdict_allmeth,
                                   docdict['callparams'],
                                   _docdict_extra_params])

docdict_discrete['before_notes'] = ''.join([_docdict_rank_method,
                                            _docdict_discrete_allmeth,
                                            docdict['callparams'],
                                            _docdict_extra_params])

_doc_translate_args = \
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

_doc_rvs_alt = \
"""
Alternative random number generator for discrete distributions.  Uses the
model's cdf function and a uniform random number generator.  Can be faster than
native scipy rvs for some custom models.  Will perform well if the the models
cdf function is also fast.

Parameters
----------
%(shapes)s : array_like
    shape parameters
l : int
    Lower bound of distribution (Either 0 or 1).  Default is 1
b : int
    Upper bound of distribution for computational purposes, even if
    distribution technically has infinite support. Default is 1e5.
size : int
    Number of random variables to draw.  Default is 1.

Returns
-------
array
    Random variables from model

Notes
-----
"""

_doc_fit_mle = \
"""
Return MLEs for shape parameters from data

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

_doc_rank = \
"""
Return predicted rank abundance distribution

Parameters
----------
n : int
    Number of values to return
%(shapes)s : array_like
    shape parameters

Returns
-------
array
    Values of rank abundance distribution

Notes
-----
Describe 0.5 offset. References.

"""
# TODO: Finish doc_rank above

_doc_make_rank = \
"""
obj : discrete distribution object
    Scipy discrete distribution object
crit : float
    A probability between 0 - 1.  Below this value ppf is used, above a this
    value a solver is used.
upper : int
    Upper bound to the solver.  Rank will not return values above
    upper
xtol : float
    Precision of the brentq solver.
"""


class rv_continuous_meco(rv_continuous):
    """
    A modified generic continuous random variable class meant for subclassing.

    This class inherits from the `rv_continuous` class of `scipy.stats` and
    contains all of its functionality. See the docstring of `rv_continuous` for
    information on usage and subclassing. In addition, this class adds two new
    methods.

    Methods
    -------
    translate_args
        Shape parameters given user-friendly parameters (see notes)
    fit_mle
        Shape parameters given data and optional keyword arguments (see notes)
    rank
        Rank abundance distribution

    """

    @doc_sub(_doc_translate_args)
    def translate_args(self, *args):
        """{0}"""
        raise NotImplementedError, ("translate_args method not implemented "
                                    "for this distribution")

    @doc_sub(_doc_fit_mle)
    def fit_mle(self, *args):
        """{0}"""
        return self.fit(*args, floc=0, fscale=1)[:-2]

    @doc_sub(_doc_rank)
    def rank(self, n, *args):
        """{0}"""
        return self.ppf((np.arange(1, n+1) - 0.5) / n, *args)


class rv_discrete_meco(rv_discrete):
    """
    A modified generic discrete random variable class meant for subclassing.

    This class inherits from the `rv_discrete` class of `scipy.stats` and
    contains all of its functionality. See the docstring of `rv_discrete` for
    information on usage and subclassing. In addition, this class adds two new
    methods.

    Methods
    -------
    translate_args
        Shape parameters given user-friendly parameters (see notes)
    fit_mle
        Shape parameters given data and optional keyword arguments (see notes)
    rank
        Rank abundance distribution

    """


    @doc_sub(_doc_translate_args)
    def translate_args(self, *args):
        """{0}"""
        raise NotImplementedError, ("translate_args method not implemented "
                                    "for this distribution")

    @doc_sub(_doc_fit_mle)
    def fit_mle(self, *args):
        """{0}"""
        raise NotImplementedError, ("fit_mle method not implemented "
                                    "for this distribution")

    @doc_sub(_doc_rank)
    def rank(self, n, *args):
        """{0}"""
        return self.ppf((np.arange(1, n+1) - 0.5) / n, *args)

    @doc_sub(_doc_rvs_alt)
    def rvs_alt(self, *args, **kwargs):
        """{0}"""
        l = kwargs.get('l', 1)
        b = kwargs.get('b', 1e5)
        size = kwargs.get('size', 1)

        model_cdf = self.cdf(np.arange(l, b + 1), *args)

        unif_rands = np.random.random(size)
        model_rands = np.array([np.where(tx <= model_cdf)[0][0] + l
                            for tx in unif_rands])

        return model_rands


#
# Discrete
#

class geom_gen(rv_discrete_meco):
    r"""
    A geometric discrete random variable.

    This implementation of the geometric distribution differs from that in
    `scipy.stats`, as the distribution here has support from 0 to inf.

    .. math::
       P(x) = (1-p)^{x} p

    for ``x >= 0``. The ``loc`` parameter is not used.

    Methods
    -------
    translate_args(mu)
        Shape parameter p given distribution mean.
    fit_mle(data)
        ML estimate of shape parameter p given data.
    %(before_notes)s
    mu : float
        distribution mean

    Examples
    --------
    >>> import macroeco.models as md

    >>> # Get the geom_parameters from a mean
    >>> mu = 20
    >>> p = md.geom.translate_args(mu)
    0.047619047619047616

    >>> # Get the pmf
    >>> md.geom.pmf(np.arange(0, 5), p)
    array([ 0.04761905,  0.04535147,  0.04319188,  0.04113512,  0.03917631])

    >>> # Generate a rank abundance distribution
    >>> rad = md.geom.rank(20, p)
    >>> rad
    array([  0.,   1.,   2.,   3.,   5.,   6.,   8.,   9.,  11.,  13.,  15.,
        17.,  20.,  23.,  26.,  30.,  35.,  42.,  53.,  75.])

    >>> # Fit the geom to data
    >>> md.geom.fit_mle(rad)
    (0.048309178743961352,)

    """

    @inherit_docstring_from(rv_discrete_meco)
    def translate_args(self, mu):
        return 1 / (np.array(mu) + 1)

    @inherit_docstring_from(rv_discrete_meco)
    def fit_mle(self, data):
        return self.translate_args(np.mean(data)),

    def _argcheck(self, p):
        return (p <= 1) & (p >= 0)

    def _pmf(self, x, p):
        return (1-p)**x * p

    def _logpmf(self, x, p):
        return x*np.log(1-p) + np.log(p)

    def _cdf(self, x, p):
        x = np.floor(x)
        return (1.0-(1.0-p)**(x+1))

    def _stats(self, p):
        mu = (1.0 - p) / p
        var = (1.0 - p) / p**2
        return mu, var, None, None

geom = geom_gen(name='geom', shapes='p')


class geom_uptrunc_gen(rv_discrete_meco):
    r"""
    An upper-truncated geometric discrete random variable.

    .. math::

       P(x) = \frac{(1-p)^{x} p}{1 - (1-p)^{b+1}}

    for ``x >= 0``. ``geom_uptrunc`` takes two shape parameters: ``p`` and
    ``b``, the upper limit. The ``loc`` parameter is not used.

    Methods
    -------
    translate_args(mu, b)
        Shape parameter p given distribution mean and upper limit.
    fit_mle(data, b=sum(data))
        ML estimate of shape parameter p given data and upper limit.
    %(before_notes)s
    mu : float
        distribution mean
    b : float
        distribution upper limit, defaults to sum of data

    Notes
    -----
    The boundary ``p = 1`` is a special case in which the ratio between
    successive terms of the distribution is 1 (i.e., the pmf is uniform). This
    arises when the mean of the distribution is precisely one-half the upper
    limit.

    This distribution is known as the Pi distribution in the MaxEnt Theory of
    Ecology [#]_, where the ``p`` parameter is equivalent to  ``1 -
    exp(-lambda)``. The special case of a uniform pmf has been described as
    HEAP [#]_.

    References
    ----------
    .. [#]
       Harte, J. (2011). Maximum Entropy and Ecology: A Theory of
       Abundance, Distribution, and Energetics (p. 264). Oxford, United
       Kingdom: Oxford University Press.
    .. [#]
       Harte, J., Conlisk, E., Ostling, A., Green, J. L., & Smith, A. B.
       (2005). A theory of spatial structure in ecological communities at
       multiple spatial scales. Ecological Monographs, 75(2), 179-197.

    Examples
    --------
    >>> import macroeco.models as md

    >>> # Get the geom parameters from a mean and upper limit
    >>> mu = 20; b = 200
    >>> p, b = md.geom_uptrunc.translate_args(mu, b)
    >>> p, b
    (array(0.047592556687674925), 200)


    >>> # Get the pmf
    >>> md.geom_uptrunc.pmf(np.arange(0, 5), p, b)
    array([ 0.04759519,  0.04533002,  0.04317264,  0.04111795,  0.03916104])

    >>> # Generate a rank abundance distribution
    >>> rad = md.geom_uptrunc.rank(20, p, b)
    >>> rad
    array([  0.,   1.,   2.,   3.,   5.,   6.,   8.,   9.,  11.,  13.,  15.,
        17.,  20.,  23.,  26.,  30.,  35.,  42.,  53.,  75.])

    >>> # Fit the geom to data
    >>> md.geom_uptrunc.fit_mle(rad)
    (0.048309175638750035, 394.0)

    """
    # TODO: Should add a warning for b < 5 or 10 or so (p solver gives erratic
    # answers. (This may or may not still be true.)

    @inherit_docstring_from(rv_discrete_meco)
    def translate_args(self, mu, b):
        return _geom_solve_p_from_mu_vect(mu, b), b

    @inherit_docstring_from(rv_discrete_meco)
    def fit_mle(self, data, b=None):
        """%(super)s
        In addition to data, requires ``b``, the upper limit of the
        distribution.
        """
        # Take mean of data as MLE of distribution mean, then calculate p
        mu = np.mean(data)
        if not b:
            b = np.sum(data)
        p = _geom_solve_p_from_mu_vect(mu, b)

        # Just return float, not len 1 array
        if len(np.atleast_1d(p)) == 1:
            return float(p), b
        else:
            return p, b

    def _argcheck(self, p, b):
        # Unlike the traditional geometric, p can be < 0
        return (p <= 1)

    def _pmf(self, x, p, b):
        pmf = (1.0-p)**x * p / (1.0-(1.0-p)**(b+1))
        if len(np.atleast_1d(x)) > 1:
            pmf[x > b] = 0
        elif x > b:
            pmf = 0
        return pmf

    def _cdf(self, x, p, b):
        x = np.floor(x)
        cdf = (1.0-(1.0-p)**(x+1)) / (1.0-(1.0-p)**(b+1))
        if len(np.atleast_1d(x)) > 1:
            cdf[x > b] = 1
        elif x > b:
            cdf = 1
        return cdf

    def _stats(self, p, b):
        mu = (p / (1 - p)) - ((b + 1) / (p**-b - 1))
        return mu, None, None, None

geom_uptrunc = geom_uptrunc_gen(name='geom_uptrunc', shapes='p, b')

def _geom_solve_p_from_mu(mu, b):
    """
    For the geom_uptrunc, given mu and b, return p.
    Ref: Harte 2011, Oxford U Press. Eq. 7.50.
    """

    def p_eq(x, mu, b):
        x, mu, b = Decimal(x), Decimal(mu), Decimal(b)
        return ( (x / (1 - x)) - ((b + 1) / (x**-b - 1)) - mu )

    # x here is the param raised to the k_agg power, or 1 - p
    return 1 - optim.brentq(p_eq, 1e-16, 100, args=(mu, b), disp=True)

_geom_solve_p_from_mu_vect = np.vectorize(_geom_solve_p_from_mu)


class dgamma_gen(rv_discrete_meco):
    r"""
    A discrete gamma random variable.

    .. math::

       P(x) = k * x^{(\alpha - 1)} * e^{(-1 / \theta)*x}

    for ``x >= 1``, ``\theta > 0``. ``k`` is the normalizing constant.
    ``\alpha`` is analogous to the ``k_{agg}`` parameter in the zero-truncated
    negative binomial.

    Methods
    -------
    translate_args(alpha, theta)
        not used, returns alpha and theta.
    fit_mle(data)
        ml estimate of shape parameters alpha and theta given data
    %(before_notes)s
    alpha : float
        distribution parameter
    theta : float
        distribution parameter

    Notes
    -----
    This parameterization of the discrete gamma was taken from [#]_.

    References
    ----------
    .. [#]
       Frank, F. (2011). Measurement scale in maximum entropy models of species
       abundance. Journal of Evolutionary Biology, 24(3), 485-496

    Examples
    --------

    >>> import macroeco.models as md

    >>> # dgamma takes two parameters
    >>> dgamma_dist = md.dgamma(alpha=1, theta=2)

    >>> # When alpha = 1 and theta = 2, should be similar to a to
    >>> # zero_truncated NBD with mu = 2.541.
    >>> dgamma_dist.pmf(np.arange(1, 10))
    array([ 0.39346934,  0.23865122,  0.14474928,  0.08779488,  0.05325028,
        0.03229793,  0.01958968,  0.01188174,  0.00720664])

    >>> md.nbinom_ztrunc.pmf(np.arange(1, 10), mu=2.541, k_agg=1)
    array([ 0.39354585,  0.23866751,  0.1447409 ,  0.08777872,  0.05323377,
        0.03228384,  0.01957867,  0.01187357,  0.00720077])

    >>> # Get the approximate mean and variance
    >>> dgamma_dist.stats()
    (array(2.541494082536799), array(3.917698089032762))

    >>> # Draw random numbers from the discrete gamma distribution
    >>> dgamma_dist.rvs(size=20)
    array([5, 4, 2, 1, 1, 1, 3, 3, 3, 3, 1, 1, 2, 1, 1, 1, 3, 1, 1, 3])


    """

    @inherit_docstring_from(rv_discrete_meco)
    def translate_args(self, alpha, theta):
        return alpha, theta

    @inherit_docstring_from(rv_discrete_meco)
    def fit_mle(self, data, init_vals=(80, 80)):
        """%(super)s
        In addition to data, can take init_vals which allows the user to
        specify initial values for (alpha, theta) during the optimization.

        """

        if len(data) > 1:
            mu = np.mean(data)
            var = np.var(data)
            theta0 = var / mu
            alpha0 = mu / theta0
        else:
            alpha0 = init_vals[0]
            theta0 = init_vals[1]

        def mle(params):
            return -np.sum(np.log(self.pmf(data, params[0], params[1])))

        # Bounded fmin?
        alpha, theta = optim.fmin(mle, x0=[alpha0, theta0], disp=0)

        return alpha, theta

    def _pmf(self, x, alpha, theta):

        b = 1e5
        alpha = np.atleast_1d(alpha)
        theta = np.atleast_1d(theta)
        b = np.atleast_1d(b)
        x = np.atleast_1d(x)

        eq = lambda val, talpha, ttheta: np.exp((talpha - 1) * np.log(val) -
            (val / ttheta))

        # eq = lambda val, talpha, ttheta: val**(talpha - 1) * \
        #                                 np.exp((-1 / ttheta)*val)

        norm = np.sum(eq(np.arange(1, b[0] + 1), alpha[0], theta[0]))

        pmf = eq(x, alpha, theta) / norm
        return pmf

    def _cdf(self, x, alpha, theta):

        alpha = np.atleast_1d(alpha)
        theta = np.atleast_1d(theta)
        x = np.atleast_1d(x)

        max_x = np.max(x)
        pmf_list = self.pmf(np.arange(1, np.int(max_x) + 1), alpha[0],
                        theta[0])
        full_cdf = np.cumsum(pmf_list)

        cdf = np.array([full_cdf[tx - 1] if tx != 0 else 0 for tx in x])

        return cdf

    def _argcheck(self, alpha, theta):

        # TODO: Can theta or alpha be 0 in the discrete version?
        return (theta > 0)

    def _stats(self, alpha, theta):
        # Fixed upper limit
        upper = 10000
        vals = np.arange(1, upper)
        pmf_vals = self.pmf(vals, alpha, theta)
        mom1 = np.sum(vals * pmf_vals)
        mom2 = np.sum(vals**2 * pmf_vals)
        var_est = mom2 - mom1**2
        return mom1, var_est, None, None

dgamma = dgamma_gen(name='dgamma', shapes='alpha, theta')


class nbinom_gen(rv_discrete_meco):
    r"""
    A negative binomial discrete random variable.

    This implementation of the negative binomial distribution differs from that
    in `scipy.stats`, as the distribution here uses the more common ecological
    parameterization.

    .. math::

       P(x) = \frac{\gamma (k + x)}{\gamma(k) x!}
       \left(\frac{k}{k+\mu}\right)^k \left(\frac{\mu}{k+\mu}\right)^x

    for ``x >= 0``. In the traditional parameterization, ``n = k_agg`` (the
    size parameter) and ``p = k_agg / (k_agg + mu)``. the ``loc`` parameter is
    not used.

    Methods
    -------
    translate_args(mu, k_agg)
        not used, returns mu and k_agg.
    fit_mle(data, k_range=(0.1,100,0.1))
        ml estimate of shape parameters mu and k_agg given data
    %(before_notes)s
    mu : float
        distribution mean
    k_agg : float
        clustering parameter

    Examples
    --------
    >>> import macroeco.models as md

    >>> # Define a NBD distribution with mean = 10 and aggregation = 2
    >>> nbd_dist = md.nbinom(mu=10, k_agg=2)

    >>> # Get the pmf for some values
    >>> nbd_dist.pmf(range(1, 10))
    array([ 0.0462963 ,  0.05787037,  0.06430041,  0.0669796 ,  0.0669796 ,
        0.06511905,  0.06201814,  0.05814201,  0.05383519])

    >>> # Get the cdf for some values
    >>> nbd_dist.cdf(range(1, 10))
    array([ 0.07407407,  0.13194444,  0.19624486,  0.26322445,  0.33020405,
        0.3953231 ,  0.45734124,  0.51548325,  0.56931845])

    >>> # Get the logpmf using a different notation
    >>> md.nbinom.logpmf(range(1, 10), 10, 2)
    array([-3.07269331, -2.84954976, -2.74418925, -2.70336725, -2.70336725,
       -2.73153813, -2.78032829, -2.84486682, -2.92182786])

    >>> # Get a random sample
    >>> samp = md.nbinom.rvs(mu=10, k_agg=1, size=10)
    >>> samp
    array([12,  1,  4, 10, 23,  0, 12,  4,  1, 15])


    >>> # Get the rank abundance distribution for n = 20
    >>> rad = md.nbinom.rank(20, 10, 1)
    >>> rad
    array([  0.,   0.,   0.,   1.,   1.,   2.,   3.,   3.,   4.,   5.,   6.,
         7.,   9.,  10.,  12.,  14.,  17.,  20.,  26.,  37.])


    """

    @inherit_docstring_from(rv_discrete_meco)
    def translate_args(self, mu, k_agg):
        return mu, k_agg

    @inherit_docstring_from(rv_discrete_meco)
    def fit_mle(self, data, k_array=np.arange(0.1, 100, 0.1)):
        """%(super)s

        In addition to data, gives an optional keyword argument k_array
        containing the values to search for k_agg. A brute force search is then
        used to find the parameter k_agg.

        """
        # todo: check and mention in docstring biases of mle for k_agg
        data = np.array(data)
        mu = np.mean(data)
        return mu, _solve_k_from_mu(data, k_array, nbinom_nll, mu)

    def _get_p_from_mu(self, mu, k_agg):
        return k_agg / (k_agg + mu)

    def _rvs(self, mu, k_agg):
        p = self._get_p_from_mu(mu, k_agg)
        return nprand.negative_binomial(k_agg, p, self._size)

    def _argcheck(self, mu, k_agg):
        p = self._get_p_from_mu(mu, k_agg)
        return (k_agg >= 0) & (p >= 0) & (p <= 1)

    def _pmf(self, x, mu, k_agg):
        return np.exp(self._logpmf(x, mu, k_agg))

    def _logpmf(self, x, mu, k_agg):
        p = self._get_p_from_mu(mu, k_agg)

        coeff =\
           special.gammaln(k_agg+x)-special.gammaln(x+1)-special.gammaln(k_agg)

        return coeff + k_agg*np.log(p) + x*np.log(1-p)

    def _cdf(self, x, mu, k_agg):
        p = self._get_p_from_mu(mu, k_agg)
        x = np.floor(x)
        return special.betainc(k_agg, x+1, p)

    def _ppf(self, q, mu, k_agg):
        p = self._get_p_from_mu(mu, k_agg)
        vals = np.ceil(special.nbdtrik(q, k_agg, p))
        vals1 = (vals-1).clip(0.0, np.inf)
        temp = self._cdf(vals1, k_agg, p)
        return np.where(temp >= q, vals1, vals)

    def _stats(self, mu, k_agg):
        p = self._get_p_from_mu(mu, k_agg)
        Q = 1.0 / p
        p = Q - 1.0
        mu = k_agg*p
        var = k_agg*p*Q
        g1 = (Q+p)/np.sqrt(k_agg*p*Q)
        g2 = (1.0 + 6*p*Q) / (k_agg*p*Q)
        return mu, var, g1, g2

nbinom = nbinom_gen(name='nbinom', shapes='mu, k_agg')


def nbinom_nll(data, k_agg, mu):
    return -np.sum(nbinom._logpmf(data, mu, k_agg))


class nbinom_ztrunc_gen(rv_discrete_meco):
    r"""
    The zero-truncated negative binomial random variable.

    This distribution is described by Sampford (1955) [#]_.

    .. math::

       P(x) = \frac{(k + x - 1)!}{(k - 1)!x!} \left(\frac{p}
        {1 + p}\right)^{x} \frac{1}{(1 + p)^{k - 1}}

    for ``x >= 1``. ``p`` can be computed directly from the mean of the
    distribution and is calculated internally so that the distribution is
    parameterized by ``\mu`` and ``k_agg`` analogous to ``nbinom``.

    Methods
    -------
    translate_args(mu, k_agg, return_p=False)
        Returns mu and k_agg. Returns p parameter if return_p is True.
    fit_mle(data, k_agg0=0.5)
        ml estimate of shape parameters mu and k_agg given data
    %(before_notes)s
    mu : float
        distribution mean
    k_agg : float
        clustering parameter

    Notes
    -----

    References
    ----------
    .. [#]
       Sampford, M. R. (1955). The truncated negative binomial distribution.
       Biometrika, 42(1), 58-69

    Examples
    --------

    >>> import macroeco.models as md

    >>> # Define a zero_truncated NBD distribution with mean = 10 and aggregation = 2
    >>> nbd_dist = md.nbinom_ztrunc(mu=10, k_agg=2)

    >>> # Get the pmf for some values
    >>> nbd_dist.pmf(range(1, 10))
    array([ 0.04984472,  0.06199534,  0.06854036,  0.07104033,  0.07068624,
        0.06838018,  0.06479938,  0.06044661,  0.05569011])

    >>> # Get the cdf for some values
    >>> nbd_dist.cdf(range(1, 10))
    array([ 0.04984472,  0.11184006,  0.18038041,  0.25142075,  0.32210698,
        0.39048717,  0.45528654,  0.51573315,  0.57142326])

    >>> # Get the logpmf using a different notation
    >>> md.nbinom_ztrunc.logpmf(range(1, 10), 10, 2)
    array([-2.99884273, -2.78069611, -2.68033253, -2.64450747, -2.64950441,
       -2.68267222, -2.73645932, -2.80599478, -2.88795275])

    >>> # Get and fit random sample
    >>> samp = md.nbinom_ztrunc.rvs(mu=10, k_agg=1, size=100)
    >>> md.nbinom_ztrunc.fit_mle(samp)
    (10.210000000000001, 0.93369140625000113)


    >>> # Get the rank abundance distribution for n = 20
    >>> rad = md.nbinom_ztrunc.rank(20, 10, 1)
    >>> rad
    array([  1.,   1.,   2.,   2.,   3.,   4.,   4.,   5.,   6.,   7.,   8.,
         9.,  10.,  11.,  13.,  15.,  17.,  20.,  25.,  36.])

    """

    @inherit_docstring_from(rv_discrete_meco)
    def translate_args(self, mu, k_agg, return_p=False):
        """%(super)s

        The keyword argument return_p computes the p values used to define the
        the truncated negative binomial
        """
        if return_p:
            return nbinom_ztrunc_p(mu, k_agg), k_agg
        else:
            return mu, k_agg

    @inherit_docstring_from(rv_discrete_meco)
    def fit_mle(self, data, k_agg0=0.5):
        """%(super)s

        In addition to data, gives an optional keyword argument k_agg0 that
        specifies the initial value of k_agg used in the optimization.

        """

        mu = np.mean(data)

        def mle(k):

            return -np.sum(np.log(self.pmf(data, mu, k)))

        k = optim.fmin(mle, x0=k_agg0, disp=0)

        return mu, k[0]

    def _pmf(self, x, mu, k_agg):

        x = np.atleast_1d(x)

        norm = np.exp(special.gammaln(k_agg + x) - ((special.gammaln(k_agg) +
                                        special.gammaln(x + 1))))
        p = nbinom_ztrunc_p(mu, k_agg)
        kernel = (p / (1 + p))**x * (1 / ((1 + p)**k_agg - 1))
        pmf = norm * kernel

        pmf[x == 0] = 0

        return pmf

    def _stats(self, mu, k_agg):
        p = nbinom_ztrunc_p(mu, k_agg)
        omega = 1 / (1 + p)
        eta = 1 - omega
        mu = mu

        # From Sampford 1955
        var = (k_agg * eta * (1 + k_agg * eta)) / \
                (omega**2 * (1 - omega**k_agg)) - mu**2
        return mu, var, None, None

nbinom_ztrunc = nbinom_ztrunc_gen(name='nbinom_ztrunc', shapes='mu, k_agg')


def _nbinom_ztrunc_p(mu, k_agg):
        """ Calculates p parameter for truncated negative binomial

        Function given in Sampford 1955, equation 4

        Note that omega = 1 / 1 + p in Sampford
        """

        p_eq = lambda p, mu, k_agg: (k_agg * p) / (1 - (1 + p)**-k_agg) - mu

        # The upper bound needs to be large. p will increase with increasing mu
        # and decreasing k_agg
        p = optim.brentq(p_eq, 1e-10, 1e10, args=(mu, k_agg))
        return p

nbinom_ztrunc_p = np.vectorize(_nbinom_ztrunc_p)


class cnbinom_gen(rv_discrete_meco):
    r"""
    The conditional negative binomial random variable.

    This distribution was described by Zillio and He (2010) [#]_ and Conlisk
    et al. (2007) [#]_

    .. math::

       P(x) = \frac{\binom{x + k - 1}{x}  \binom{b - x + k/a - k -1}{b
                -x}}{\binom{b + k/a - 1}{b}}

    for ``x >= 0``. In this parameterization ``a = E[p(x)] / b`` where ``b`` is
    the upper limit of the distribution.

    Methods
    -------
    translate_args(mu, k_agg, b)
        not used, returns mu, k_agg, and b.
    fit_mle(data, k_array=np.arange(0.1,100,0.1))
        ml estimate of shape parameters mu and k_agg given data
    %(before_notes)s
    mu : float
        distribution mean
    k_agg : float
        clustering parameter (refered to as ``k`` above)
    b : float
        Upper bound of distribution

    References
    ----------
    .. [#]
        Zillio, T. & He, F. (2010). Modeling spatial aggregation of finite
        populations. Ecology, 91(12), 3698-3706
    .. [#]
        Conlisk, E., Bloxham, M., Conlisk, J, Enquist, E., and Harte, J.
        (2007). A new class of models of spatial distribution. Ecological
        Monographs, 77(2), 269-284

    Examples
    --------

    >>> import macroeco.models as md

    >>> # Define a conditional NBD distribution with mean = 10, agg = 2,
    >>> # and an upper bound (b) = 300
    >>> cnbd_dist = md.cnbinom(mu=10, k_agg=2, b=300)

    >>> # Get the pmf for some values
    >>> cnbd_dist.pmf(range(0, 10))
    array([ 0.02662579,  0.04474923,  0.05637649,  0.06309932,  0.06617407,
        0.06658649,  0.06510469,  0.06232244,  0.05869438,  0.05456466])

    >>> # Get the cdf for some values
    >>> cnbd_dist.cdf(range(0, 10))
    array([ 0.02662579,  0.07137502,  0.12775151,  0.19085083,  0.2570249 ,
        0.32361139,  0.38871607,  0.45103851,  0.50973289,  0.56429755])

    >>> # Get the logpmf using a different notation
    >>> md.cnbinom.logpmf(range(1, 10), 10, 2, 300)
    array([-3.10668105, -2.8757031 , -2.76304533, -2.71546655, -2.7092536 ,
       -2.73175874, -2.7754338 , -2.83541131, -2.90836891])

    >>> # Get and fit random sample
    >>> samp = md.cnbinom.rvs(mu=10, k_agg=1, b=300, size=100)
    >>> md.cnbinom.fit_mle(samp)
    (11.640000000000001, 1.2000000000000002, 1164)

    >>> # Be more specific about fit_mle
    >>> md.cnbinom.fit_mle(samp, k_array=np.linspace(1, 1.5, num=1000))
    (11.640000000000001, 1.1966966966966968, 1164)


    >>> # Get the rank abundance distribution for n = 20
    >>> rad = md.cnbinom.rank(20, 10, 1, 300)
    >>> rad
    array([  0.,   0.,   1.,   2.,   2.,   3.,   4.,   5.,   5.,   6.,   7.,
         9.,  10.,  11.,  13.,  15.,  18.,  21.,  26.,  37.])


    """

    @inherit_docstring_from(rv_discrete_meco)
    def translate_args(self, mu, k_agg, b):
        return mu, k_agg, b

    @inherit_docstring_from(rv_discrete_meco)
    def fit_mle(self, data, b=None, k_array=np.arange(0.1, 100, 0.1)):

        data = np.array(data)
        mu = np.mean(data)

        if not b:
            b = np.sum(data)

        return mu, _solve_k_from_mu(data, k_array, _cnbinom_nll, mu, b), b

    def _pmf(self, x, mu, k_agg, b):
        return np.exp(self._logpmf(x, mu, k_agg, b))

    def _logpmf(self, x, mu, k_agg, b):
        a = mu / b
        logpmf = _cnbinom_logpmf(x, b, a, k_agg)
        logpmf[x > b] = -np.inf
        return logpmf

    def _stats(self, mu, k_agg, b):
        mu = mu
        var = ((1 - mu / b) * mu * (k_agg + mu)) / (k_agg + (mu / b))
        return mu, var, None, None

cnbinom = cnbinom_gen(name="cnbinom", shapes="mu, k_agg, b")


def _cnbinom_logpmf(n_i, n, a, k_agg):
    # Logpmf for cnbinom
    return _ln_choose(n_i + k_agg - 1, n_i) + \
        _ln_choose(n - n_i + (k_agg / a) - k_agg - 1, n - n_i) -\
        _ln_choose(n + (k_agg / a) - 1, n)


def _cnbinom_nll(data, k_agg, mu, b):
    # Negative log likelihood for cnbinom
    return -np.sum(cnbinom._logpmf(data, mu, k_agg, b))


def _ln_choose(n, k_agg):
    '''
    log binomial coefficient with extended gamma factorials. n and k_agg may be
    int or array - if both array, must be the same length.

    '''
    gammaln = special.gammaln
    return gammaln(n + 1) - (gammaln(k_agg + 1) + gammaln(n - k_agg + 1))


def _solve_k_from_mu(data, k_array, nll, *args):
    """
    For given args, return k_agg from searching some k_range.

    Parameters
    ----------
    data : array
    k_range : array
    nll : function

    args :

    Returns
    --------
    :float
        Minimum k_agg

    """
    # TODO: See if a root finder like fminbound would work with Decimal used in
    # logpmf method (will this work with arrays?)

    nll_array = np.zeros(len(k_array))

    for i in range(len(k_array)):
        nll_array[i] = nll(data, k_array[i], *args)

    min_nll_idx = np.argmin(nll_array)

    return k_array[min_nll_idx]

class logser_gen(rv_discrete_meco):
    """
    Logseries (logarithmic) random variable.

    .. math::

        P(x) = - p**x / (x*log(1-p))

    Methods
    -------
    translate_args(mu)
        Translates the mean into p
    fit_mle(data)
        ml estimate of shape parameter p
    %(before_notes)s
    p : float
        p parameter of the logseries distribution

    Examples
    --------

    >>> import macroeco.models as md

    >>> # Define a logseries distribution by specifying the necessary parameters
    >>> logser_dist = md.logser(p=0.9)

    >>> # Get the pmf
    >>> logser_dist.pmf(1)
    0.39086503371292664

    >>> # Get the cdf
    >>> logser_dist.cdf(10)
    0.9201603889810761

    >>> # You can also use the following notation
    >>> md.logser.pmf(1, 0.9)
    0.39086503371292664
    >>> md.logser.cdf(10, 0.9)
    0.9201603889810761

    >>> # Get a rank abundance distribution for 30 species
    >>> rad = md.logser.rank(30, 0.9)
    >>> rad
    array([  1.,   1.,   1.,   1.,   1.,   1.,   1.,   1.,   1.,   1.,   1.,
     1.,   2.,   2.,   2.,   2.,   2.,   3.,   3.,   3.,   4.,   4.,
     5.,   5.,   6.,   7.,   8.,  10.,  13.,  21.])

    >>> # Fit the logser to data and estimate the parameters
    >>> md.logser.fit_mle(rad)
    (0.8957380268455059,)

    """

    @inherit_docstring_from(rv_continuous_meco)
    def translate_args(self, mu):
        eq = lambda p, mu: -p/np.log(1-p)/(1-p) - mu
        return optim.brentq(eq, 1e-16, 1-1e-16, args=(mu), disp=True)

    @inherit_docstring_from(rv_continuous_meco)
    def fit_mle(self, data):
        # Use method of moments
        return self.translate_args(np.mean(data)),

    def _rvs(self, p):
        # looks wrong for p>0.5, too few k=1
        # trying to use generic is worse, no k=1 at all
        return stats.logser.rvs(p, size=self._size)
        #return np.random.mtrand.logseries(p, size=self._size)

    def _argcheck(self, p):
        return (p > 0) & (p < 1)

    def _pmf(self, x, p):
        return stats.logser.pmf(x, p)
        # return -np.power(p, x) * 1.0 / x / np.log(1 - p)

    def _cdf(self, x, p):
        return stats.logser.cdf(x, p)

    def _stats(self, p):
        r = np.log(1 - p)
        mu = p / (p - 1.0) / r
        mu2p = -p / r / (p - 1.0)**2
        var = mu2p - mu*mu
        mu3p = -p / r * (1.0+p) / (1.0 - p)**3
        mu3 = mu3p - 3*mu*mu2p + 2*mu**3
        g1 = mu3 / np.power(var, 1.5)

        mu4p = -p / r * (
            1.0 / (p-1)**2 - 6*p / (p - 1)**3 + 6*p*p / (p-1)**4)
        mu4 = mu4p - 4*mu3p*mu + 6*mu2p*mu*mu - 3*mu**4
        g2 = mu4 / var**2 - 3.0
        return mu, var, g1, g2

logser = logser_gen(name="logser", shapes="p")


class logser_uptrunc_gen(rv_discrete_meco):
    r"""
    Upper truncated logseries random variable.

    This distribution was described by Harte (2011) [#]_

    .. math::

        p(x) = \frac{1}{Z} \frac{p^n}{n}

    where ``Z`` is the normalizing factor

    Methods
    -------
    translate_args(mu, b)
        Translates the mean and the upper bound into p and b.
    fit_mle(data)
        ml estimate of shape parameter p
    %(before_notes)s
    p : float
        p parameter of the logseries distribution
    b : float
        Upper bound of the distribution


    Notes
    -----
    Code adapted from Ethan White's macroecology_tools and version 0.1 of
    macroeco

    References
    -----------
    .. [#]
        Harte, J. (2011). Maximum Entropy and Ecology: A Theory of
        Abundance, Distribution, and Energetics. Oxford, United
        Kingdom: Oxford University Press.

    Examples
    --------

    >>> import macroeco.models as md

    >>> # Define a logseries distribution by specifying the necessary parameters
    >>> logser_dist = md.logser_uptrunc(p=0.9, b=1000)

    >>> # Get the pmf
    >>> logser_dist.pmf(1)
    0.39086503371292664

    >>> # Get the cdf
    >>> logser_dist.cdf(10)
    0.9201603889810761

    >>> # You can also use the following notation
    >>> md.logser_uptrunc.pmf(1, 0.9, 1000)
    0.39086503371292664
    >>> md.logser_uptrunc.cdf(10, 0.9, 1000)
    0.9201603889810761

    >>> # Get a rank abundance distribution for 30 species
    >>> rad = md.logser_uptrunc.rank(30, 0.9, 1000)
    >>> rad
    array([  1.,   1.,   1.,   1.,   1.,   1.,   1.,   1.,   1.,   1.,   1.,
         1.,   2.,   2.,   2.,   2.,   2.,   3.,   3.,   3.,   4.,   4.,
         5.,   5.,   6.,   7.,   8.,  10.,  13.,  21.])

    >>> # Fit the logser_uptrunc to data and estimate the parameters
    >>> md.logser_uptrunc.fit_mle(rad)
    (0.8957385644316679, 114.0)

    """

    @inherit_docstring_from(rv_discrete_meco)
    def translate_args(self, mu, b):
        return _trunc_logser_solver((1 / mu) * b, b), b

    @inherit_docstring_from(rv_discrete_meco)
    def fit_mle(self, data, b=None):
        """%(super)s
b : float
    The upper bound of the distribution. If None, fixed at sum(data)
        """

        data = np.array(data)
        length = len(data)

        if not b:
            b = np.sum(data)

        return _trunc_logser_solver(length, b), b

    def _pmf(self, x, p, b):

        x = np.atleast_1d(x)
        p = np.atleast_1d(p)
        b = np.atleast_1d(b)

        if p[0] > 0:
            pmf = stats.logser.pmf(x, p) / stats.logser.cdf(b, p)
        else:
            ivals = np.arange(1, b[0] + 1)
            normalization = np.sum(p[0] ** ivals / ivals)
            pmf = (p[0] ** x / x) / normalization

        return pmf

    def _cdf(self, x, p, b):

        x = np.atleast_1d(x)
        p = np.atleast_1d(p)
        b = np.atleast_1d(b)

        if p[0] < 1:
            return stats.logser.cdf(x, p) / stats.logser.cdf(b, p)
        else:
            cdf_list = [sum(self.pmf(range(1, int(x_i) + 1), p[0], b[0])) for
                                                x_i in x]
            return np.array(cdf_list)

    def _rvs(self, p, b):
        # Code from weecology/macroecotools

        if not self._size:
            self._size = 1

        out = []
        if p < 1:
            for i in range(self._size):
                rand_logser = stats.logser.rvs(p)
                while rand_logser > b:
                    rand_logser = stats.logser.rvs(p)
                out.append(rand_logser)
        else:
            rand_list = stats.uniform.rvs(size = self._size)
            for rand_num in rand_list:
                y = lambda x: self.cdf(x, p, b) - rand_num
                if y(1) > 0: out.append(1)
                else: out.append(int(round(bisect(y, 1, b))))
        return np.array(out)

    def _stats(self, p, b):

        vals = np.arange(1, b + 1)
        full_pmf = self.pmf(vals, p, b)
        mean, var = _mean_var(vals, full_pmf)
        return mean, var, None, None


logser_uptrunc = logser_uptrunc_gen(name="logser_uptrunc", shapes="p, b")


def _trunc_logser_solver(bins, b):
    """
    Given bins (S) and b (N) solve for MLE of truncated logseries
    parameter p

    Parameters
    -----------
    bins : float
        Number of bins. Considered S in an ecological context
    b : float
        Upper truncation of distribution. Considered N in an ecological context

    Returns
    -------
    : float
        MLE estimate of p

    Notes
    ------
    Adapted from Ethan White's macroecology_tools
    """

    if bins == b:
        p = 0

    else:
        BOUNDS = [0, 1]
        DIST_FROM_BOUND = 10 ** -15
        m = np.array(np.arange(1, np.int(b) + 1))
        y = lambda x: np.sum(x ** m / b * bins) - np.sum((x ** m) / m)
        p = optim.bisect(y, BOUNDS[0] + DIST_FROM_BOUND,
                   min((sys.float_info[0] / bins) ** (1 / b), 2),
                   xtol=1.490116e-08, maxiter=1000)
    return p


class plnorm_gen(rv_discrete_meco):
    r"""
    Poisson lognormal random variable.

    Methods
    -------
    translate_args(mean, sigma)
        not implemented
    fit_mle(data)
        ml estimate of shape parameters mu and sigma
    %(before_notes)s
    mu : float
        mu parameter of the poisson lognormal
    sigma : float
        sigma parameter of the poisson lognormal

    Notes
    -----
    The pmf method was adopted directly from the VGAM package in R.
    The VGAM R package was adopted directly from Bulmer (1974) [#]_

    The fit_mle function was adapted from Ethan White's pln_solver function in
    macroeco_distributions (https://github.com/weecology/macroecotools)

    References
    ----------
    .. [#]
        Bulmer, M. G. (1974). On fitting the poisson lognormal distribution to
        species bundance data. Biometrics, 30, 101-110.

    Examples
    --------

    >>> import macroeco.models as md

    >>> # Get the pmf for the poisson lognormal with mu = -1 and sigma = 3
    >>> md.plnorm.pmf(np.arange(1, 11), -1, 3)
    array([ 0.12139284,  0.05769201,  0.03558665,  0.02486353,  0.01868109,
        0.01472104,  0.01199807,  0.01002759,  0.00854552,  0.00739661])

    >>> md.plnorm.pmf([0, 50, 1000], 2.34, 5)
    array([  2.86468926e-01,   1.51922299e-03,   5.25717609e-05])

    >>> # Get the CDF
    >>> md.plnorm.cdf([0, 15, 10000], mu=.1, sigma=2)
    array([ 0.3954088 ,  0.90489995,  0.99999662])

    >>> # Rank abundance distribution
    >>> md.plnorm.rank(10, 1, 1, crit=0.5, upper=40)
    array([  0.,   0.,   1.,   2.,   2.,   4.,   5.,   7.,   8.,  15.])

    >>> # Fit the the plnorm to data
    >>> data = np.array([1,1,1,1,1,2,2,2,3,3,4,4,5,5,6,6,12,45,67])
    >>> md.plnorm.fit_mle(data)
    (1.3195513537335777, 1.1876220131629682)


    """

    @inherit_docstring_from(rv_discrete_meco)
    def translate_args(self, mean, sigma):
        raise NotImplementedError("Translate args not implemented")

    @inherit_docstring_from(rv_discrete_meco)
    def fit_mle(self, data):

        mu0 = np.mean(np.log(np.array(data) + 1))
        sig0 = np.std(np.log(np.array(data) + 1))

        if sig0 == 0:

            sig0 = 1e-5 # can't be zero

        def mle(params):
            return -np.sum(self.logpmf(data, params[0], params[1]))

        # Bounded fmin?
        mu, sigma = optim.fmin_bfgs(mle, x0=[mu0, sig0], disp=0)

        return mu, sigma

    @inherit_docstring_from(rv_discrete_meco)
    @doc_sub(_doc_make_rank)
    def rank(self, n, mu, sigma, crit=.5, upper=10000, xtol=1):
        """%(super)s

Additional Parameters
----------------------
        {0}

        """

        return _make_rank(self, n, mu, sigma, crit=crit, upper=upper,
                                                                    xtol=xtol)

    def _argcheck(self, mu, sigma):
        return True

    def _pmf(self, x, mu, sigma):

        # TODO: Add approx_cut as keyword. Strange parse_args error
        approx_cut = 10
        x = np.array(x)
        pmf = np.empty(len(x), dtype=np.float)
        xbelow = x <= approx_cut
        xabove = x > approx_cut

        # If below, use exact answer
        if np.sum(xbelow) > 0:

            pmf[xbelow] = plognorm_intg_vec(x[xbelow], mu[xbelow],
                                                                sigma[xbelow])

        # If above, use approximation
        if np.sum(xabove) > 0:

            z = (np.log(x[xabove]) - mu[xabove]) / sigma[xabove]

            pmf_above = ((1 + (z**2 + np.log(x[xabove]) - mu[xabove] - 1) /
                    (2 * x[xabove] * sigma[xabove]**2)) * np.exp(-0.5 * z**2) /
                    (np.sqrt(2 * np.pi) * sigma[xabove] * x[xabove]))

            pmf[xabove] = pmf_above

        # If pmf is 0 the likelihood might break
        # TODO: This should be fixed in likelihood function as it might apply
        # to other distributions
        pmf[pmf == 0] = 1e-120

        return pmf

    def _cdf(self, x, mu, sigma):

        mu = np.atleast_1d(mu)
        sigma = np.atleast_1d(sigma)
        x = np.atleast_1d(x)

        max_x = np.max(x)
        pmf_list = self.pmf(np.arange(np.int(max_x) + 1), mu[0], sigma[0])
        full_cdf = np.cumsum(pmf_list)

        cdf = np.array([full_cdf[tx] for tx in x])

        return cdf


plnorm = plnorm_gen(name='plnorm', shapes='mu,sigma')


class plnorm_ztrunc_gen(rv_discrete_meco):
    r"""
    Zero-truncated poisson lognormal random variable.

    Methods
    -------
    translate_args(mean, sigma)
        not implemented
    fit_mle(data)
        ml estimate of shape parameters mu and sigma
    %(before_notes)s
    mu : float
        mu parameter of the poisson lognormal
    sigma : float
        sigma parameter of the poisson lognormal

    Notes
    -----
    The pmf method was adopted directly from the VGAM package in R.
    The VGAM R package was adopted directly from Bulmer (1974) [#]_

    The fit_mle function was adapted from Ethan White's pln_solver function in
    macroeco_distributions (https://github.com/weecology/macroecotools)

    References
    ----------
    .. [#]
        Bulmer, M. G. (1974). On fitting the poisson lognormal distribution to
        species bundance data. Biometrics, 30, 101-110.

    Examples
    --------

    >>> import macroeco.models as md

    >>> # Get the pmf for the zero-truncated poisson lognormal with mu = -1
    >>> # and sigma = 3
    >>> md.plnorm_ztrunc.pmf(np.arange(1, 11), -1, 3)
    array([ 0.27334111,  0.12990549,  0.08013071,  0.05598538,  0.04206434,
        0.03314746,  0.02701614,  0.02257919,  0.019242  ,  0.01665499])

    >>> md.plnorm_ztrunc.pmf([1, 50, 1000], 2.34, 5)
    array([  9.27474648e-02,   2.12916164e-03,   7.36783061e-05])

    >>> # Get the CDF
    >>> md.plnorm_ztrunc.cdf([1, 15, 10000], mu=.1, sigma=2)
    array([ 0.27575055,  0.84270355,  0.99999442])

    >>> # Rank abundance distribution. For speed, set crit = 0 or else the
    >>> # calculation is very slow.
    >>> md.plnorm_ztrunc.rank(20, 1, 1, crit=0, upper=40)
    array([  1.,   1.,   2.,   2.,   2.,   2.,   2.,   3.,   3.,   4.,   5.,
         5.,   6.,   6.,   7.,   7.,   8.,  11.,  14.,  22.])

    >>> # Fit the the plnorm to data
    >>> data = np.array([1,1,1,1,1,2,2,2,3,3,4,4,5,5,6,6,12,45,67])
    >>> md.plnorm_ztrunc.fit_mle(data)
    (0.19056468723097392, 1.7698965441710266)

    """

    @inherit_docstring_from(rv_discrete_meco)
    def translate_args(self, mean, sigma):
        raise NotImplementedError("Translate args not implemented")

    @inherit_docstring_from(rv_discrete_meco)
    def fit_mle(self, data):

        # Copying code...could we make this a generic function with an eval?
        # Or would that slow it down too much?
        mu0 = np.mean(np.log(data))
        sig0 = np.std(np.log(data))

        if sig0 == 0:

            sig0 = 1e-5 # can't be zero

        def mle(params):
            return -np.sum(np.log(self._pmf(data, params[0], params[1])))

        # Bounded fmin?
        mu, sigma = optim.fmin_bfgs(mle, x0=[mu0, sig0], disp=0)

        return mu, sigma

    @inherit_docstring_from(rv_discrete_meco)
    @doc_sub(_doc_make_rank)
    def rank(self, n, mu, sigma, crit=0, upper=10000, xtol=1):
        """%(super)s

Additional Parameters
----------------------
        {0}

        """

        return _make_rank(self, n, mu, sigma, crit=crit, upper=upper,
                                                                    xtol=xtol)
    def _argcheck(self, mu, sigma):
        return True

    def _pmf(self, x, mu, sigma):

        x = np.array(x)
        mu = np.atleast_1d(mu)
        sigma = np.atleast_1d(sigma)

        norm = 1 - plognorm_intg_vec(0, mu[0], sigma[0])
        pmf_vals = plnorm.pmf(x, mu, sigma) / norm
        pmf_vals[x < 1] = 0

        return pmf_vals

    def _cdf(self, x, mu, sigma):

        # Format input
        x = np.array(x)
        mu = np.atleast_1d(mu)
        sigma = np.atleast_1d(sigma)

        # Calculate cdf from plnorm_gen
        norm = 1 - plognorm_intg_vec(0, mu[0], sigma[0])
        cdf_vals = (plnorm.cdf(x, mu, sigma) -
                                        plnorm.cdf(0, mu[0], sigma[0])) / norm

        # Values less than one have zero probability
        cdf_vals = np.atleast_1d(cdf_vals)
        cdf_vals[x < 1] = 0

        return cdf_vals

plnorm_ztrunc = plnorm_ztrunc_gen(name="plnorm_ztrunc",
        shapes='mu, sigma')


def plognorm_intg(x, mu, sigma):
    # Integral for plognorm
    eq = lambda t, x, mu, sigma: np.exp(t * x - np.exp(t) - 0.5 *
                                        ((t - mu) / sigma) ** 2)

    intg = integrate.quad(eq, -np.inf, np.inf, args=(x, mu, sigma))[0]

    norm = np.exp(-0.5 * np.log(2 * np.pi * sigma ** 2) -
                            special.gammaln(x + 1))

    return norm * intg

plognorm_intg_vec = np.vectorize(plognorm_intg)



#
# Continuous
#


class expon_gen(rv_continuous_meco):
    r"""
    An exponential continuous random variable.

    .. math::

       f(x) = \lambda e^{-\lambda x}

    for ``x >= 0``. The ``loc`` and ``scale`` parameters are not used.

    Methods
    -------
    translate_args(mu)
        Shape parameter mu given distribution mean.
    fit_mle(data)
        ML estimate of shape parameter lam given data.
    %(before_notes)s
    mu : float
        distribution mean

    Examples
    --------

    >>> import macroeco.models as md
    >>> import numpy as np

    >>> # Get the rate parameter of the exponential distribution from a mean
    >>> md.expon.translate_args(20)
    0.05

    >>> # Get the pdf
    >>> md.expon.pdf(np.linspace(0.1, 10, num=10), 0.05)
    array([ 0.04975062,  0.04708823,  0.04456831,  0.04218324,  0.03992581,
        0.03778919,  0.0357669 ,  0.03385284,  0.03204121,  0.03032653])

    >>> # Get the cdf
    >>> md.expon.cdf(np.linspace(0.1, 10, num=10), 0.05)
    array([ 0.00498752,  0.05823547,  0.10863386,  0.15633518,  0.20148378,
        0.24421626,  0.28466191,  0.32294313,  0.35917572,  0.39346934])

    >>> # Get the ppf
    >>> md.expon.ppf(0.8, 0.05)
    32.188758248682014

    >>> # Draw a random sample
    >>> samp = md.expon.rvs(0.05, size=100)

    >>> # Fit the model to data
    >>> md.expon.fit_mle(samp)
    0.052277939307395938


    """

    @inherit_docstring_from(rv_continuous_meco)
    def translate_args(self, mu):
        return 1 / mu

    @inherit_docstring_from(rv_continuous_meco)
    def fit_mle(self, data):
        # MLE is method of moments for exponential
        return 1 / (np.sum(data) / len(data))

    def _rvs(self, lam):
        return nprand.exponential(1/lam, self._size)

    def _pdf(self, x, lam):
        return lam * np.exp(-lam*x)

    def _cdf(self, x, lam):
        return 1 - np.exp(-lam*x)

    def _entropy(self, lam):
        return 1 - np.ln(lam)

    def _stats(self, lam):
        return lam**-1, lam**-2, 2, 6

expon = expon_gen(a=0.0, name='expon', shapes='lam')


class expon_uptrunc_gen(rv_continuous_meco):
    r"""
    An upper-truncated exponential continuous random variable.

    .. math::

       f(x) = \frac{\lambda e^{-\lambda x}}{1 - e^{-\lambda b}}

    for ``b >= x >= 0``. The ``loc`` and ``scale`` parameters are not used.

    Methods
    -------
    translate_args(mu, b)
        Shape parameter lam given distribution mean and upper limit.
    fit_mle(data, b=sum(data))
        ML estimate of shape parameter lam given data and upper limit.
    %(before_notes)s
    mu : float
        distribution mean
    b : float
        distribution upper limit, defaults to sum of data

    Examples
    ---------
    >>> import macroeco.models as md
    >>> import numpy as np

    >>> # Get the rate parameter of the exponential distribution from a mean
    >>> md.expon_uptrunc.translate_args(20, 100)
    (array(0.04801007549722518), 100)

    >>> # Get the pdf
    >>> md.expon_uptrunc.pdf(np.linspace(0.1, 10, num=10), 0.05, 100)
    array([ 0.05008812,  0.04740766,  0.04487064,  0.0424694 ,  0.04019665,
        0.03804554,  0.03600953,  0.03408249,  0.03225857,  0.03053226])

    >>> # Get the cdf
    >>> md.expon_uptrunc.cdf(np.linspace(0.1, 10, num=10), 0.05, 100)
    array([ 0.00502135,  0.05863052,  0.10937079,  0.15739571,  0.20285058,
        0.24587294,  0.28659296,  0.32513386,  0.36161225,  0.3961385 ])

    >>> # Get the ppf
    >>> md.expon_uptrunc.ppf(0.8, 0.05, 100)
    31.656858541834165

    >>> # Draw a random sample
    >>> samp = md.expon_uptrunc.rvs(0.05, 100, size=100)

    >>> # Fit the model to data
    >>> md.expon_uptrunc.fit_mle(samp)
    (0.06080396315704938, 1644.6296393823973)

    """

    @inherit_docstring_from(rv_continuous_meco)
    def translate_args(self, mu, b):
        return _expon_solve_lam_from_mu_vect(mu, b), b

    @inherit_docstring_from(rv_continuous_meco)
    def fit_mle(self, data, b=None):
        """%(super)s

        Additional Parameters
        ----------------------
        b : float
            The upper limit of the distribution
        """
        # Take mean of data as MLE of distribution mean, then calculate p
        mu = np.mean(data)
        if not b:
            b = np.sum(data)
        lam = _expon_solve_lam_from_mu_vect(mu, b)

        # Just return float, not len 1 array
        if len(np.atleast_1d(lam)) == 1:
            return float(lam), b
        else:
            return lam, b

    def _argcheck(self, lam, b):
        return True

    def _pdf(self, x, lam, b):
        return (lam * np.exp(-lam*x)) / (1 - np.exp(-lam*b))

    def _cdf(self, x, lam, b):
        return (1 - np.exp(-lam*x)) / (1 - np.exp(-lam*b))

expon_uptrunc = expon_uptrunc_gen(a=0.0, name='expon_uptrunc', shapes='lam, b')

def _expon_solve_lam_from_mu(mu, b):
    """
    For the expon_uptrunc, given mu and b, return lam.
    Similar to geom_uptrunc
    """

    def lam_eq(lam, mu, b):
        # Small offset added to denominator to avoid 0/0 erors
        lam, mu, b = Decimal(lam), Decimal(mu), Decimal(b)
        return ( (1 - (lam*b + 1) * np.exp(-lam*b)) /
                 (lam - lam * np.exp(-lam*b) + Decimal(1e-32)) - mu )

    return optim.brentq(lam_eq, -100, 100, args=(mu, b), disp=True)

_expon_solve_lam_from_mu_vect = np.vectorize(_expon_solve_lam_from_mu)


class lognorm_gen(rv_continuous_meco):
    r"""
    A lognormal random variable.

    .. math::

        f(x) = \frac{1}{\sigma x \sqrt{2 \pi}} e^{(\log{x} - \mu)^2 / 2
        \sigma^2}

    Methods
    -------
    translate_args(mean, sigma)
        Shape parameters mu and sigma given mean and sigma
    fit_mle(data, b=sum(data))
        ML estimate of shape parameters mu and sigma
    %(before_notes)s
    mu : float
        mu parameter of lognormal distribution. Mean log(x)
    sigma : float
        sigma parameter of lognormal distribution. sd of log(x)

    Examples
    --------

    >>> import macroeco.models as md
    >>> import numpy as np

    >>> # Given mean = 20 and sigma = 2, get the parameters for the lognormal
    >>> md.lognorm.translate_args(20, 2)
    (0.99573227355399085, 2)

    >>> # Define a lognormal distribution
    >>> lgnorm = md.lognorm(mu=0.99573227355399085, sigma=2)

    >>> # Get the PDF
    >>> lgnorm.pdf(np.linspace(0.1, 100, num=10))
    array([  5.12048368e-01,   1.38409265e-02,   5.13026600e-03,
         2.71233007e-03,   1.68257934e-03,   1.14517879e-03,
         8.28497146e-04,   6.26026253e-04,   4.88734158e-04,
         3.91396174e-04])

    >>> # Get the stats for the lognormal
    >>> lgnorm.stats()
    (array(19.99935453933589), array(21437.87621568711))

    >>> # Similarly you could use
    >>> md.lognorm.stats(mu=0.99573227355399085, sigma=2)
    (array(20.0), array(21439.260013257688))

    >>> # Draw some random numbers from the lognormal
    >>> samp = md.lognorm.rvs(mu=1.5, sigma=1.3, size=100)

    >>> # Fit model to data
    >>> md.lognorm.fit_mle(samp)
    (1.2717334369626212, 1.3032723732257057)

    >>> # Get the rank abundance distribution for a lognormal for 10 species
    >>> md.lognorm.rank(10, 1.5, 1.3)
    array([  0.52818445,   1.16490157,   1.86481775,   2.71579138,
         3.806234  ,   5.27701054,   7.39583206,  10.77077745,
        17.24226102,  38.02750505])

    """

    @inherit_docstring_from(rv_continuous_meco)
    def translate_args(self, mean, sigma):
        return np.log(mean) - (sigma ** 2 / 2), sigma

    @inherit_docstring_from(rv_continuous_meco)
    def fit_mle(self, data, fix_mean=False):
        """%(super)s

        Additional Parameters
        ----------------------
        fix_mean : bool
            Default False.  If True, fixes mean before optimizing sigma

        """

        if not fix_mean:
            sigma, _, scale = stats.lognorm.fit(data, floc=0)
            return np.log(scale), sigma

        else:
            mean = np.mean(data)

            # MLE fxn to be optmimized
            mle = lambda sigma, x, mean: -1 *\
                                    np.sum(self._pdf_w_mean(x, mean, sigma))

            sigma = optim.fmin(mle, np.array([np.std(np.log(data), ddof=1)]),
                                            args=(data, mean), disp=0)[0]

            return self.translate_args(mean, sigma)

    def _pdf_w_mean(self, x, mean, sigma):
        """
        Calculates the pdf of a lognormal distribution with parameters mean
        and sigma

        Parameters
        ----------
        mean : float or ndarray
            Mean of the lognormal distribution
        sigma : float or ndarray
            Sigma parameter of the lognormal distribution

        Returns
        -------
        : float or ndarray
            pdf of x
        """

        # Lognorm pmf with mean for optimization
        mu, sigma = self.translate_args(mean, sigma)
        return self.logpdf(x, mu, sigma)

    def _argcheck(self, mu, sigma):
        return True

    def _rvs(self, mu, sigma):
        return stats.lognorm.rvs(sigma, scale=np.exp(mu), size=self._size)

    def _pdf(self, x, mu, sigma):
        return stats.lognorm.pdf(x, sigma, scale=np.exp(mu))

    def _cdf(self, x, mu, sigma):
        return stats.lognorm.cdf(x, sigma, scale=np.exp(mu))

    def _stats(self, mu, sigma):
        return stats.lognorm.stats(sigma, scale=np.exp(mu), moments="mvsk")

lognorm = lognorm_gen(name="lognorm", shapes="mu, sigma")


@doc_sub(_doc_make_rank)
def _make_rank(dist_obj, n, mu, sigma, crit=0.5, upper=10000, xtol=1):
    """
    Make rank distribution using both ppf and brute force.

    Setting crit = 1 is equivalent to just using the ppf

    Parameters
    ----------
    {0}

    """
    qs = (np.arange(1, n + 1) - 0.5) / n
    rank = np.empty(len(qs))

    brute_ppf = lambda val, prob: prob - dist_obj.cdf(val, mu, sigma)

    qs_less = qs <= crit
    ind = np.sum(qs_less)

    # Use ppf if qs are below crit
    rank[qs_less] = dist_obj.ppf(qs[qs_less], mu, sigma)

    # Use brute force if they are above
    for i, tq in enumerate(qs[~qs_less]):

        j = ind + i
        try:
            # TODO: Use an adaptable lower bound to increase speed
            rank[j] = np.abs(np.ceil(optim.brentq(brute_ppf, -1, upper,
                                        args=(tq,), xtol=xtol)))

        except ValueError:

            # If it is above the upper bound set all remaining values
            # to the previous value
            rank[j:] = np.repeat(rank[j - 1], len(rank[j:]))
            break

    return rank


def _mean_var(vals, pmf):
    """
    Calculates the mean and variance from vals and pmf

    Parameters
    ----------
    vals : ndarray
        Value range for a distribution
    pmf : ndarray
        pmf values corresponding with vals

    Returns
    -------
    : tuple
        (mean, variance)

    """

    mean = np.sum(vals * pmf)
    var = np.sum(vals ** 2 * pmf) - mean ** 2
    return mean, var


