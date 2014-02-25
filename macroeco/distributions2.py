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

from scipy.misc.doccer import inherit_docstring_from
from scipy.stats.distributions import (rv_discrete, rv_continuous, docdict, 
                                       docdict_discrete, docheaders)
import scipy.stats.distributions as spdist
import scipy.optimize as optim
import scipy.special as spec


_doc_param_note = \
"""There are many available methods of `%(name)s`, each of which require one or 
more of the parameters listed below.
"""

_doc_custom_methods = \
"""fit2(data, %(shapes)s)
    MLE estimates of shapes given initial guesses (use instead of `fit`)."""

_doc_discrete_custom_methods = \
"""translate_args(uargs)
    Get shape parameters from user-friendly args.
fit2(data, %(shapes)s)
    MLE estimates of shapes given initial guesses."""

# Remove header from all methods
_docdict_allmeth_sh = docdict['allmethods'][16:]
_docdict_discrete_allmeth_sh = docdict_discrete['allmethods'][17:]

# **kwds in expect string followed by no space was throwing warning
_docdict_allmeth_sh = _docdict_allmeth_sh.replace(', **kwds','')

docdict['before_notes'] = ''.join([_doc_param_note,
                                   docheaders['methods'],
                                   _doc_custom_methods,
                                   _docdict_allmeth_sh,
                                   docdict['callparams']])

docdict_discrete['before_notes'] = ''.join([_doc_param_note, 
                                            docheaders['methods'],
                                            _doc_discrete_custom_methods,
                                            _docdict_discrete_allmeth_sh,
                                            docdict_discrete['callparams']])



class rv_continuous_meco(rv_continuous):
    """
    A modified generic continuous random variable class meant for subclassing.

    This class inherits from the `rv_continuous` class of `scipy.stats` and 
    contains all of its functionality. See the docstring of `rv_continuous` for 
    information on usage and subclassing. In addition, this class adds one new 
    methods.

    Methods
    -------
    fit2
        calls method `fit` with fixed loc=0 and scale=1 (defaults)

    """

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
        takes user-friendly params as input and returns shape params
    fit
        estimates distribution params from data

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
            Subset of shape parameters that are not fit. See Notes.

        Returns
        -------
        tuple of floats
            MLEs for shape parameters
    
        Notes
        -----
        """

        raise NotImplementedError, ("fit method not implemented for this "
                                    "distribution")


#
# Discrete
#

class geom_gen(rv_discrete_meco):
    r"""
    A geometric discrete random variable.

    This implementation of the geometric distribution differs from that in 
    `scipy.stats`, as the distribution here has support from 0 to inf.

    .. math::
       \mathrm{pmf(x)} = (1-p)^{x} p

    for ``x >= 0``. The location parameter ``loc`` is not used.

    %(before_notes)s
    uargs : float
        distribution mean

    """

    @inherit_docstring_from(rv_discrete_meco)
    def translate_args(self, mu):
        return 1 / (np.array(mu) + 1)

    @inherit_docstring_from(rv_discrete_meco)
    def fit2(self, data):
        """%(super)s
        Requires one argument containing data to fit.
        """
        return self.translate_args(np.mean(data)), 

    def _argcheck(self, p):
        return (p <= 1) & (p >= 0)

    def _pmf(self, x, p):
        return (1-p)**x * p

    def _logpmf(self, x, p):
        return k*np.log(1-p) + log(p)

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

       \mathrm{pmf(x)} = \frac{(1-p)^{x} p}{1 - (1-p)^{b+1}}

    for ``x >= 0``.
    
    `geom_uptrunc` takes two shape parameters: ``p`` and ``b``, the upper 
    limit. The location parameter ``loc`` is not used.

    %(before_notes)s

    uargs : float
        distribution mean, upper limit

    Notes
    -----
    The boundary ``p = 1`` is a special case in which the ratio between 
    successive terms of the distribution is 1 (i.e., the pmf is uniform). This 
    arises when the mean of the distribution is precisely one-half the upper 
    limit.

    This distribution is known as the Pi distribution in the MaxEnt Theory of 
    Ecology [#]_, where the ``p`` parameter is known as ``exp(-lambda)``.

    References
    ----------
    .. [#]
       Harte, J. (2011). Maximum Entropy and Ecology: A Theory of
       Abundance, Distribution, and Energetics (p. 264). Oxford, United
       Kingdom: Oxford University Press.

    ..
       DEV: There is a difficult implicit equation needed to determine the p 
       parameter from the mu and b arguments. We've employed the brentq solver 
       here but note that it fails regularly for certain shape combinations.

    """

    # TODO: Should add a warning for b < 5 or 10 or so (p solver gives erratic 
    # answers.

    @inherit_docstring_from(rv_discrete_meco)
    def translate_args(self, mu, b):
        return _geom_solve_p_from_mu_vect(mu, b), b

    @inherit_docstring_from(rv_discrete_meco)
    def fit2(self, data, b):
        """%(super)s
        Requires two arguments consisting of data to fit and ``b``, the upper 
        limit of the distribution (held constant).
        """
        # Take mean of data as MLE of distribution mean, then calculate p
        mu = np.mean(data)
        return _geom_solve_p_from_mu_vect(mu, b), b

    def _argcheck(self, p, b):
        # Unlike the traditional geometric, p can be > 0
        return (p >= 0)

    def _pmf(self, x, p, b):
        pmf = (1.0-p)**x * p / (1.0-(1.0-p)**(b+1))
        pmf[x > b] = 0
        return pmf

    def _cdf(self, x, p, b):
        k = np.floor(x)
        cdf = (1.0-(1.0-p)**(x+1)) / (1.0-(1.0-p)**(b+1))
        cdf[x > b] = 1
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

    def p_eq(p, mu, b):
        p, mu, b = Decimal(p), Decimal(mu), Decimal(b)
        return ( (p / (1 - p)) - ((b + 1) / (p**-b - 1)) - mu )

    return optim.brentq(p_eq, 1e-9, 20, args=(mu, b), disp=True)

_geom_solve_p_from_mu_vect = np.vectorize(_geom_solve_p_from_mu)


class nbinom_gen(spdist.nbinom_gen):
    r"""
    A negative binomial discrete random variable.

    This implementation of the geometric distribution differs from that in 
    `scipy.stats`, as the distribution here uses the more common ecological 
    parameterization.

    .. math::
        
       \mathrm{pmf(x)} =
       \frac{\Gamma (k + x)}{\Gamma(k) x!} \left(\frac{k}{k+\mu}\right)^k
       \left(\frac{\mu}{k+\mu}\right)^x

    for ``x >= 0``. In the traditional parameterization, ``n = k`` (the size 
    parameter) and ``p = k / (k + mu)``. The location parameter ``loc`` is not 
    used.

    %(before_notes)s

    """

    @inherit_docstring_from(rv_discrete_meco)
    def fit2(self, x, k_range=(0.1,100,0.1)):
        """%(super)s
        Requires one argument containing data to fit. A keyword argument 
        k_range contains a tuple of the start, stop, and step values to search 
        for k. Default is ``k_range=(0.1,100,0.1)``.

        This method recognizes that the MLE of the mu parameter is simply equal 
        to the mean of the data. A brute force search is then used to find the 
        parameter k.

        """
        assert len(x) > 20, "nbinom fit is not stable with <20 data points"
        mu = np.mean(x)
        return mu, _nbinom_solve_k_from_mu(x, mu, k_range)

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
        coeff = spec.gammaln(k+x) - spec.gammaln(x+1) - spec.gammaln(k)
        return coeff + k*np.log(p) + x*np.log(1-p)

    def _cdf(self, x, mu, k):
        p = self._get_p_from_mu(mu, k)
        x = np.floor(x)
        return spec.betainc(k, x+1, p)

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

def _nbinom_solve_k_from_mu(x, mu, k_range):
    """
    For the nbinom, given mu, return k from searching some k_range.
    """

    # TODO: See if a root finder like fminbound would work with Decimal used in 
    # logpmf method (will this work with arrays?)

    def nll(x, mu, k):
        return -np.sum(nbinom._logpmf(x, mu, k))

    k_array = np.arange(*k_range)
    nll_array = np.zeros(len(k_array))

    for i in range(len(k_array)):
        nll_array[i] = nll(x, mu, k_array[i])

    min_nll_idx = np.argmin(nll_array)

    return k_array[min_nll_idx]

#
# Continuous
#

class expon_gen(rv_continuous_meco):
    r"""
    An exponential continuous random variable.

    .. math::
        
       \mathrm{pdf(x)} = \lambda e^{-\lambda x}

    for ``x >= 0``.

    %(before_notes)s

    """

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
        
       \mathrm{pdf(x)} = \frac{\lambda e^{-\lambda x}}{1 - e^{-\lambda x}}

    for ``b >= x >= 0``.

    %(before_notes)s

    """

    # Internally, class works by creating a new expon_gen object with the
    # appropriate upper limit and calling its methods.

    # TODO: Do all of these broadcast correctly, or should we call _pdf, etc.?

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

    @inherit_docstring_from(rv_discrete_meco)
    def fit2(self, data, lam, b):
        expon = expon_gen(a=0.0, b=b)
        return expon.fit(data, lam, floc=0, fscale=1), b

expon_uptrunc = expon_uptrunc_gen(a=0.0, name='expon_uptrunc', shapes='lam, b')




