from __future__ import division
import sys

from decimal import Decimal
import numpy as np
import numpy.random as nprand

from scipy.stats.distributions import (rv_discrete, rv_continuous, docdict,
                                       docdict_discrete)
import scipy.stats as stats
import scipy.optimize as optim
import scipy.special as special
import scipy.integrate as integrate

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

class nbinom_gen(rv_discrete_meco):
    r"""
    A negative binomial discrete random variable.

    This implementation of the negative binomial distribution differs from that
    in `scipy.stats`, as the distribution here uses the more common ecological
    parameterization.

    .. math::

       p(x) = \frac{\gamma (k + x)}{\gamma(k) x!}
       \left(\frac{k}{k+\mu}\right)^k \left(\frac{\mu}{k+\mu}\right)^x

    for ``x >= 0``. in the traditional parameterization, ``n = k_agg`` (the
    size parameter) and ``p = k_agg / (k_agg + mu)``. the ``loc`` parameter is
    not used.

    Methods
    -------
    translate_args(mu, k_agg)
        not used, returns mu and k_agg.
    fit_mle(data, k_range=(0.1,100,0.1))
        ml estimate of shape parameters mu and k_agg given data, with k_agg evaluated at (min, max, step) values given by k_range.
    %(before_notes)s
    mu : float
        distribution mean
    k_agg : float
        clustering parameter

    """

    @inherit_docstring_from(rv_discrete_meco)
    def translate_args(self, mu, k_agg):
        return mu, k_agg

    @inherit_docstring_from(rv_discrete_meco)
    def fit_mle(self, data, k_range=(0.1, 100, 0.1)):
        """%(super)s

        In addition to data, gives an optional keyword argument k_range
        contains a tuple of the start, stop, and step values to search for
        k_agg. default is ``k_range=(0.1,100,0.1)``. a brute force search is
        then used to find the parameter k_agg.

        """
        # todo: check and mention in docstring biases of mle for k_agg
        mu = np.mean(data)
        return mu, _solve_k_from_mu(data, k_range, nbinom_nll, mu)

    def _get_p_from_mu(self, mu, k_agg):
        return k_agg / (k_agg + mu)

    def _rvs(self, mu, k_agg):
        p = self._get_p_from_mu(mu, k_agg)
        return nprand.negative_binomial(k_agg, p, self._size)

    def _argcheck(self, mu, k_agg):
        p = self._get_p_from_mu(mu, k_agg)
        return (k_agg >= 0) & (p >= 0) & (p <= 1)

    def _pmf(self, x, mu, k_agg):
        p = self._get_p_from_mu(mu, k_agg)
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


class cnbinom_gen(rv_discrete_meco):
    r"""
    The conditional negative binomial random variable.

    This distribution was described by Zillio and He (2010) [#]_ and Conlisk
    et al. (2007) [#]_

    .. math::

       p(x) = \frac{\binom{x + k - 1}{x}  \binom{b - x + k/a - k -1}{b
                -x}}{\binom{b + k/a - 1}{b}}

    for ``x >= 0``. In this parameterization ``a = E[p(x)] / b`` where ``b`` is
    the upper limit of the distribution.

    Methods
    -------
    translate_args(mu, k_agg, b)
        not used, returns mu, k_agg, and b.
    fit_mle(data, k_range=(0.1,100,0.1))
        ml estimate of shape parameters mu and k_agg given data, with k_agg evaluated at (min, max, step) values given by k_range.
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
    """

    @inherit_docstring_from(rv_discrete_meco)
    def translate_args(self, mu, k_agg, b):
        return mu, k_agg, b

    @inherit_docstring_from(rv_discrete_meco)
    def fit_mle(self, data, b=None, k_range=(0.1, 100, 0.1)):

        data = np.array(data)
        mu = np.mean(data)

        if not b:
            b = np.sum(data)

        return mu, _solve_k_from_mu(data, k_range, _cnbinom_nll, mu, b), b

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


def _solve_k_from_mu(data, k_range, nll, *args):
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

    k_array = np.arange(*k_range)
    nll_array = np.zeros(len(k_array))

    for i in range(len(k_array)):
        nll_array[i] = nll(data, k_array[i], *args)

    min_nll_idx = np.argmin(nll_array)

    return k_array[min_nll_idx]


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


    """

    @inherit_docstring_from(rv_discrete_meco)
    def translate_args(self, mu, b):
        return _trunc_logser_solver((1 / mu) * b, b), b

    @inherit_docstring_from(rv_discrete_meco)
    def fit_mle(self, data, b=None):

        data = np.array(data)
        length = len(data)

        if not b:
            b = np.sum(data)

        return _trunc_logser_solver(length, b), b

    def _pmf(self, x, p, b):
        x = np.array(x)

        if p[0] > 0:
            pmf = stats.logser.pmf(x, p) / stats.logser.cdf(b, p)
        else:
            ivals = np.arange(1, b[0] + 1)
            normalization = sum(p[0] ** ivals / ivals)
            pmf = (p[0] ** x / x) / normalization

        return pmf

    def _cdf(self, x, p, b):
        x = np.array(x)
        if p[0] < 1:
            return stats.logser.cdf(x, p) / stats.logser.cdf(b, p)
        else:
            cdf_list = [sum(self.pmf(range(1, int(x_i) + 1), p[0], b[0])) for
                                                x_i in x]
            return np.array(cdf_list)

    def _rvs(self, p, b):
        # Code from weecology/macroecotools

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
        mean, var = mean_var(vals, full_pmf)
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

    """

    @inherit_docstring_from(rv_discrete_meco)
    def translate_args(self, mean, sigma):
        raise NotImplementedError("Translate args not implemented")

    @inherit_docstring_from(rv_discrete_meco)
    def fit_mle(self, data):

        mu0 = np.mean(np.log(data))
        sig0 = np.std(np.log(data))

        def mle(params):
            return -np.sum(self.logpmf(data, params[0], params[1]))

        # Bounded fmin?
        mu, sigma = optim.fmin(mle, x0=[mu0, sig0], disp=0)

        return mu, sigma

    @inherit_docstring_from(rv_discrete_meco)
    def rank(self, n, mu, sigma, **kwds):
        """%(super)s

        Uses approximation of rank distribution. The keyword ``upper`` defines
        the upper bound used in the approximation.  Default is 100000.
        """
        upper = kwds.get('upper', 100000)
        return make_rank(self.pmf(np.arange(upper + 1), mu, sigma), n,
            min_supp=0)

    def _argcheck(self, mu, sigma):
        return True

    def _pmf(self, x, mu, sigma, approx_cut=10):

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

        return pmf

    def _cdf(self, x, mu, sigma, approx_cut=10):

        mu = np.atleast_1d(mu)
        sigma = np.atleast_1d(sigma)
        x = np.atleast_1d(x)

        max_x = np.max(x)
        pmf_list = self.pmf(np.arange(np.int(max_x) + 1), mu[0], sigma[0])
        full_cdf = np.cumsum(pmf_list)

        cdf = np.array([full_cdf[tx] for tx in x])

        return cdf

    def _stats(self, mu, sigma, upper=100000):

        vals = np.arange(0, upper + 1)
        full_pmf = self.pmf(vals, mu, sigma)

        mean, var = mean_var(vals, full_pmf)

        return mean, var, None, None

plnorm = plnorm_gen(name="plnorm", shapes="mu, sigma")


class plnorm_lowtrunc_gen(rv_discrete_meco):
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

        def mle(params):
            return -np.sum(np.log(self._pmf(data, params[0], params[1])))

        # Bounded fmin?
        mu, sigma = optim.fmin(mle, x0=[mu0, sig0], disp=0)

        return mu, sigma

    @inherit_docstring_from(rv_discrete_meco)
    def rank(self, n, mu, sigma, upper=100000):
        """%(super)s
        Uses approximation of rank distribution. Increasing ``upper`` will
        give a closer approximation.
        """

        return make_rank(self.pmf(np.arange(upper + 1), mu, sigma), n,
            min_supp=1)

    def _argcheck(self, mu, sigma):
        return True

    def _pmf(self, x, mu, sigma):

        mu = np.atleast_1d(mu)
        sigma = np.atleast_1d(sigma)

        norm = 1 - plognorm_intg_vec(0, mu[0], sigma[0])
        return plnorm.pmf(x, mu, sigma) / norm

    def _cdf(self, x, mu, sigma):

        mu = np.atleast_1d(mu)
        sigma = np.atleast_1d(sigma)

        norm = 1 - plognorm_intg_vec(0, mu[0], sigma[0])
        return (plnorm.cdf(x, mu, sigma) -
                                        plnorm.cdf(0, mu[0], sigma[0])) / norm

    def _stats(self, mu, sigma, upper=100000):

        vals = np.arange(1, upper + 1)
        full_pmf = self.pmf(vals, mu, sigma)
        mean, var = mean_var(vals, full_pmf)

        return mean, var, None, None


plnorm_lowtrunc = plnorm_lowtrunc_gen(name="plnorm_lowtrunc",
        shapes="mu, sigma")

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

    """

    @inherit_docstring_from(rv_continuous_meco)
    def translate_args(self, mu):
        return 1 / mu

    @inherit_docstring_from(rv_continuous_meco)
    def fit_mle(self, data):
        expon = expon_gen(a=0.0)
        return 1 / expon.fit(data, floc=0)[2],

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

       f(x) = \frac{\lambda e^{-\lambda x}}{1 - e^{-\lambda x}}

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

    """

    # Internally, class works by creating a new expon_gen object with the
    # appropriate upper limit and calling its methods.

    # TODO: Do all of these broadcast correctly, or should we call _pdf, etc.?

    @inherit_docstring_from(rv_continuous_meco)
    def translate_args(self, mu, b):
        raise NotImplementedError, "Translation of mu to lam not implemented"

    @inherit_docstring_from(rv_continuous_meco)
    def fit_mle(self, data, b=None):
        """%(super)s

        Additional Parameters
        ----------------------
        b : float
            The upper limit of the distribution
        """
        if not b:
            b = np.sum(data)
        expon = expon_gen(a=0.0, b=b)
        return 1/expon.fit(data, floc=0)[2], b

    def _rvs(self, lam, b):
        expon = expon_gen(a=0.0, b=b)
        return expon.rvs(lam)

    def _pdf(self, x, lam, b):
        expon = expon_gen(a=0.0, b=b)
        return expon.pdf(x, lam)

    def _cdf(self, x, lam, b):
        expon = expon_gen(a=0.0, b=b)
        return expon.cdf(x, lam)

    def _entropy(self, lam, b):
        expon = expon_gen(a=0.0, b=b)
        return expon.entropy(lam)

    def _stats(self, lam, b):
        expon = expon_gen(a=0.0, b=b)
        return expon.stats(lam)

expon_uptrunc = expon_uptrunc_gen(a=0.0, name='expon_uptrunc', shapes='lam, b')


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
            sigma = optim.fmin(mle, np.array([np.std(np.log(data), ddof=1)]),
                                            args=(data, mean), disp=0)[0]
            return self.translate_args(mean, sigma)

    def _argcheck(self, mu, sigma):
        return True

    def _rvs(self, mu, sigma):
        return stats.lognorm.rvs(sigma, scale=np.exp(mu))

    def _pdf(self, x, mu, sigma):
        return stats.lognorm.pdf(x, sigma, scale=np.exp(mu))

    def _cdf(self, x, mu, sigma):
        return stats.lognorm.cdf(x, sigma, scale=np.exp(mu))

    def _stats(self, mu, sigma):
        return stats.lognorm.stats(sigma, scale=np.exp(mu))

lognorm = lognorm_gen(name="lognorm", shapes="mu, sigma")


def tpdf(x, mean, sigma):
    # Lognorm pmf with mean for optimization
    mu, sigma = lognorm.translate_args(mean, sigma)
    return lognorm.logpdf(x, mu, sigma)


def mle(sigma, x, mean):
    # MLE function for lognormal
    return -1 * np.sum(tpdf(x, mean, sigma))



def mean_var(vals, pmf):
    # Calculates the mean and variance from vals and pmf

    mean = np.sum(vals * pmf)
    var = np.sum(vals ** 2 * pmf) - mean ** 2
    return mean, var


def make_rank(pmf, n, min_supp=1):
    '''
    Convert any pmf into a rank curve for S species using cumulative
    distribution function.

    Parameters
    ----------
    pmf : ndarray
        Probability of observing a species from 1 to length pmf individs.
    n : int
        Total number of samples
    min_supp : int
        The minimum support of the distribution. Often either 1 or 0.

    Returns
    -------
    ndarray
        1D array of predicted ranks

    Notes
    -----
    Function actually implements (philosophically) a step quantile function.
    Use if ppf in rv_discrete_meco is too slow

    '''

    pmf = pmf / np.sum(pmf)  # Ensure distribution is normalized

    points = np.arange(1 / (2 * n), 1, 1 / n)
    counts = np.zeros(n)

    if min_supp == 1:
        pmf = np.array([0] + list(pmf)) # Add 0 to start of pmf
    cum_pmf = np.cumsum(pmf)

    for cutoff in cum_pmf:
        greater_thans = (points >= cutoff)
        counts[greater_thans] += 1

        if not greater_thans.any():  # If no greater thans, done with samples
            break

    return counts
