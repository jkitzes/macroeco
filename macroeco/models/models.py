from __future__ import division

from decimal import Decimal
import numpy as np
import numpy.random as nprand

from scipy.stats.distributions import (rv_discrete, rv_continuous, docdict,
                                       docdict_discrete)
import scipy.stats.distributions as spdist
import scipy.optimize as optim
import scipy.special as special

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
        return x*np.log(1-p) + log(p)

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
        pmf[x > b] = 0
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

    # x here is the param raised to the k power, or 1 - p
    return 1 - optim.brentq(p_eq, 1e-9, 20, args=(mu, b), disp=True)

_geom_solve_p_from_mu_vect = np.vectorize(_geom_solve_p_from_mu)


class nbinom_gen(spdist.nbinom_gen):
    r"""
    A negative binomial discrete random variable.

    This implementation of the negative binomial distribution differs from that
    in `scipy.stats`, as the distribution here uses the more common ecological
    parameterization.

    .. math::

       P(x) =
       \frac{\Gamma (k + x)}{\Gamma(k) x!} \left(\frac{k}{k+\mu}\right)^k
       \left(\frac{\mu}{k+\mu}\right)^x

    for ``x >= 0``. In the traditional parameterization, ``n = k`` (the size
    parameter) and ``p = k / (k + mu)``. The ``loc`` parameter is not used.

    Methods
    -------
    translate_args(mu, k)
        Not used, returns mu and k.
    fit_mle(data, k_range=(0.1,100,0.1))
        ML estimate of shape parameters mu and k given data, with k evaluated
        at (min, max, step) values given by k_range.
    %(before_notes)s
    mu : float
        distribution mean
    k : float
        clustering parameter

    """

    @inherit_docstring_from(rv_discrete_meco)
    def translate_args(self, mu, k):
        return mu, k

    @inherit_docstring_from(rv_discrete_meco)
    def fit_mle(self, data, k_range=(0.1,100,0.1)):
        """%(super)s
        In addition to data, gives an optional keyword argument
        k_range contains a tuple of the start, stop, and step values to search
        for k. Default is ``k_range=(0.1,100,0.1)``. A brute force search is
        then used to find the parameter k.

        """
        # TODO: Check and mention in docstring biases of MLE for k
        mu = np.mean(data)
        return mu, _nbinom_solve_k_from_mu(data, mu, k_range)

    def _get_p_from_mu(self, mu, k):
        return k / (k + mu)

    def _rvs(self, mu, k):
        p = self._get_p_from_mu(mu, k)
        return nprand.negative_binomial(k, p, self._size)

    def _argcheck(self, mu, k):
        p = self._get_p_from_mu(mu, k)
        return (k >= 0) & (p >= 0) & (p <= 1)

    def _pmf(self, x, mu, k):
        p = self._get_p_from_mu(mu, k)
        return np.exp(self._logpmf(x, mu, k))

    def _logpmf(self, x, mu, k):
        p = self._get_p_from_mu(mu, k)
        coeff = special.gammaln(k+x)-special.gammaln(x+1)-special.gammaln(k)
        return coeff + k*np.log(p) + x*np.log(1-p)

    def _cdf(self, x, mu, k):
        p = self._get_p_from_mu(mu, k)
        x = np.floor(x)
        return special.betainc(k, x+1, p)

    def _ppf(self, q, mu, k):
        p = self._get_p_from_mu(mu, k)
        vals = np.ceil(special.nbdtrik(q, k, p))
        vals1 = (vals-1).clip(0.0, np.inf)
        temp = self._cdf(vals1, k, p)
        return np.where(temp >= q, vals1, vals)

    def _stats(self, mu, k):
        p = self._get_p_from_mu(mu, k)
        Q = 1.0 / p
        P = Q - 1.0
        mu = k*P
        var = k*P*Q
        g1 = (Q+P)/np.sqrt(k*P*Q)
        g2 = (1.0 + 6*P*Q) / (k*P*Q)
        return mu, var, g1, g2

nbinom = nbinom_gen(name='nbinom', shapes='mu, k')

def _nbinom_solve_k_from_mu(data, mu, k_range):
    """
    For the nbinom, given mu, return k from searching some k_range.
    """
    # TODO: See if a root finder like fminbound would work with Decimal used in
    # logpmf method (will this work with arrays?)

    def nll(data, mu, k):
        return -np.sum(nbinom._logpmf(data, mu, k))

    k_array = np.arange(*k_range)
    nll_array = np.zeros(len(k_array))

    for i in range(len(k_array)):
        nll_array[i] = nll(data, mu, k_array[i])

    min_nll_idx = np.argmin(nll_array)

    return k_array[min_nll_idx]

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
        In addition to data, requires ``b``, the upper limit of the
        distribution.
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
