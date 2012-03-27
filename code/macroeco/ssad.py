#!/usr/bin/python

'''
Calculate pmf and likelihood of species-level spatial-abundance distributions.

All distributions have an argument summary, which if False returns the entire 
pmf for the inputted values of n, and if true returns the summed negative 
log-likelihood of the inputted values (useful for likelihood ratio tests or 
AIC).

Distributions
-------------
- `bin` -- binomial
- `pois` -- poisson
- `nbd` -- negative binomial (Zillio and He 2010)
- `fnbd` -- finite negative binomial (Zillio and He 2010)
- `cnbd` -- conditioned negative binomial (Conlisk et al 2007b)
- `geo` -- geometric (wrapper of nbd)
- `fgeo` -- finite negative binomial with k = 1 (wrapper of fnbd)
- `tgeo` -- truncated geometric (Harte et al 2008)

Supporting functions
--------------------
- `_ln_choose` -- log of binomial coefficient with gamma factorials
- `_ln_F` -- log of F(k, n) from Conlisk et al 2007 Am Nat

References
----------
Conlisk E, Bloxham M, Conlisk J, Enquist B, Harte J (2007a) A new class of 
models of spatial distribution. Ecological Monographs 77:269-284.

Conlisk E, Conlisk J, Harte J (2007b) The impossibility of estimating a 
negative binomial clustering parameter from presence-absence data: a comment on 
He and Gaston. The American Naturalist 170:651-654.

Harte J, Conlisk E, Ostling A, Green JL, Smith AB (2005) A theory of spatial 
structure in ecological communities at multiple spatial scales. Ecological 
Monographs 75:179-197.

Zillio T, He F (2010) Modeling spatial aggregation of finite populations. 
Ecology 91:3698-3706.
'''

from __future__ import division
import numpy as np
import scipy
import scipy.special
import scipy.stats
import scipy.optimize


__author__ = "Justin Kitzes"
__copyright__ = "Copyright 2012, Regents of the University of California"
__credits__ = ["John Harte"]
__license__ = None
__version__ = "0.1"
__maintainer__ = "Justin Kitzes"
__email__ = "jkitzes@berkeley.edu"
__status__ = "Development"


#
# DISTRIBUTIONS
#

# TODO: Check that summary is bool (decorator?) to avoid accidental confusion
# of leaving in a k when not needed and having summary = k argument.

def bin(n, N, a, summary = False):
    '''
    Binomial pmf (ie, random placement model).

    Parameters
    ----------
    n : ndarray or int
        Values of n for which to calculate pmf or likelihood
    N : ndarray or int
        Total number of individuals in landscape
    a : ndarray or int
        Ratio of cell size to area of whole landscape

    Returns
    -------
    : ndarray or float
        If summary is False, returns array with pmf. If summary is True, 
        returns the summed log likelihood of all values in n.
    '''
   
    if summary: return -sum(np.log(scipy.stats.binom.pmf(n, N, a)))
    else:       return scipy.stats.binom.pmf(n, N, a)


def pois(n, N, a, summary = False):
    '''
    Poisson pmf.

    Parameters
    ----------
    n : ndarray or int
        Values of n for which to calculate pmf or likelihood
    N : ndarray or int
        Total number of individuals in landscape
    a : ndarray or int
        Ratio of cell size to area of whole landscape

    Returns
    -------
    : ndarray or float
        If summary is False, returns array with pmf. If summary is True, 
        returns the summed log likelihood of all values in n.
    '''
    
    mu = N * a

    if summary: return -sum(np.log(scipy.stats.poisson.pmf(n, mu)))
    else:       return scipy.stats.poisson.pmf(n, mu)


def nbd(n, N, a, k, summary = False):
    '''
    Negative binomial pmf.

    Parameters
    ----------
    n : ndarray or int
        Values of n for which to calculate pmf or likelihood
    N : ndarray or int
        Total number of individuals in landscape
    a : ndarray or int
        Ratio of cell size to area of whole landscape
    k : int
        Aggregation parameter

    Returns
    -------
    : ndarray or float
        If summary is False, returns array with pmf. If summary is True, 
        returns the summed log likelihood of all values in n.
    '''

    mu = N * a
    p = 1 / (mu / k + 1)  # See Bolker book Chapt 4
    pmf = scipy.stats.nbinom.pmf(n, k, p)
    
    if summary:
        if (pmf == 0).any():
            return np.nan
        else:
            return -sum(np.log(pmf))
    else:
        return pmf


def fnbd(n, N, a, k, summary = False):
    '''
    Finite negative binomial pmf (Zillio and He 2010).

    Parameters
    ----------
    n : ndarray or int
        Values of n for which to calculate pmf or likelihood
    N : ndarray or int
        Total number of individuals in landscape
    a : ndarray or int
        Ratio of cell size to area of whole landscape
    k : ndarray or int
        Aggregation parameter

    Returns
    -------
    : ndarray or float
        If summary is False, returns array with pmf. If summary is True, 
        returns the summed log likelihood of all values in n.

    Notes
    -----
    If any arg other than n is iterable, all args must be iterable of same 
    length.

    The fnbd with k = 1 is not a true geometric distribution - calculate a pmf 
    z and run z[1:]/z[0:-1], noting that the ratio is not constant.
    '''
    
# TODO: Fix to work if n and N are one value
#    if not (n <= N).all():
#        raise Exception, "All values of n must be <= N."
#    elif (a <= 0) or (a >= 1):
#        raise Exception, "a must be between 0 and 1"

    ln_L = lambda n_i,N,a,k: _ln_choose(n_i+k-1,n_i) + \
        _ln_choose(N-n_i+(k/a)-k-1,N-n_i) - _ln_choose(N +(k/a)-1,N)

    pmf = ln_L(n, N, a, k)  # Already log

    if summary:
        return -sum(pmf)
    else:
        return np.exp(pmf)


def cnbd(n, N, a, k, summary = False):
    '''
    Conditioned negative binomial pmf (Conlisk et al 2007b)
    
    Notation modified to match Zillio and He 2010 - especially note that all 
    a's in Conlisk et al 2007 are k's in Zillio and He notation as per Theorem 
    1.5. Identical parameters and output to fnbd (see docstring there), but 
    with different calculation method.
    '''

    M = 1 / a  # Number of cells

    if not (n <= N).all():
        raise Exception, "All values of n must be <= N."
    elif (a <= 0) or (a >= 1):
        raise Exception, "a must be between 0 and 1"

    ln_F_a_n_i = np.zeros(np.size(n))  # Ln first term in num, Theorem 1.3
    for i in xrange(0,np.size(n)):
        ln_F_a_n_i[i] = _ln_F(k, n[i])

    ln_F_second = np.zeros(np.size(n))  # Ln second term in num, Theorem 1.3
    for j in xrange(0,np.size(n)):
        ln_F_a_n_i[j] = _ln_F((M - 1) * k, N - n[j])

    ln_F_Ma_N = _ln_F(M * k, N)  # Ln denom, Theorem 1.3

    if summary: return -sum(ln_F_a_n_i) - ln_F_Ma_N
    else:       return np.exp(ln_F_a_n_i + ln_F_second - ln_F_Ma_N)


def geo(n, N, a, summary = False):
    '''
    Geometric pmf (Zillio and He 2010). Wrapper for nbd with k = 1, see 
    docstring there.
    '''

    return nbd(n, N, a, 1, summary)


def fgeo(n, N, a, summary = False):
    '''
    Finite geometric pmf (Zillio and He 2010). Wrapper for fnbd with k = 1, see 
    docstring there.

    Notes
    -----
    This is not a true geometric distribution - calculate a pmf z and run 
    z[1:]/z[0:-1], noting that the ratio is not constant.
    '''

    return fnbd(n, N, a, 1, summary)


def tgeo(n, N, a, summary = False):
    '''
    Truncated geometric pmf (Harte et al 2008).

    Parameters
    ----------
    n : ndarray or int
        Values of n for which to calculate pmf or likelihood
    N : ndarray or int
        Total number of individuals in landscape
    a : ndarray or int
        Ratio of cell size to area of whole landscape
    k : int
        Aggregation parameter

    Returns
    -------
    : ndarray (1D)
        If summary is False, returns array with pmf. If summary is True, 
        returns the summed log likelihood of all values in n.

    Notes
    -----
    N and a must be the same length.

    This is a true geometric distribution in which p = exp(-lambda). Confirmed 
    with z[1:]/z[0:-1].

    Plotting eq for various realistic values of p, N, and a confirms that this 
    is a smooth, well-behaved function that should be amenable to using a root 
    finder. As written, function works (ie, root is in the interval 0,2) for 
    all tested N >= 10 and a <= 0.9. If root is outside, will throw error.

    eq is from Harte 2011, eq2 from Harte et al 2008 - both give the same 
    answer for p.
    '''

# TODO: Fix to work if n and N are one value
#    if not (n <= N).all():
#        raise Exception, "All values of n must be <= N."
#    elif (a <= 0) or (a >= 1):
#        raise Exception, "a must be between 0 and 1"

    tol = 1e-16
    N = np.round(N)  # N must be integer, brentq seems to fail if not

    eq = lambda p,N,a: (p / (1-p)) - (N+1)*(p**(N+1)) / (1 - (p**(N+1))) - N*a
    eq2 = lambda p,N,a: (1/(1-p**(N+1))) * ( (p/(1-p)) - p**(N+1) * \
            (N + (1/(1-p))) ) - N*a

    if type(N) != int and N.shape[0] > 1:  # If N is array

        if type(n) == int or n.shape[0] == 1:  # If n is an int
            n_array = np.repeat(n, N.shape[0])
        else:
            n_array = n

        pmf = np.zeros(N.shape[0])
        for i in xrange(0, N.shape[0]):
            p_i = scipy.optimize.brentq(eq, tol, 1 - tol, args = (N[i], a[i]), 
                                        disp = True)
            Z_i = (1 - p_i**(N[i]+1)) / (1 - p_i)
            pmf[i] = p_i**n_array[i] / Z_i

    else:
        p = scipy.optimize.brentq(eq, tol, 1 - tol, args = (N, a), disp = True)
        Z = (1 - p**(N+1)) / (1 - p)
        pmf = p**n / Z

    if summary: return -sum(np.log(pmf))
    else:       return pmf


#
# SUPPORTING FUNCTIONS
#

def _ln_choose(n, k):
    '''
    Log binomial coefficient with extended gamma factorials. n and k may be int 
    or array - if both array, must be the same length.
    '''
    gammaln = scipy.special.gammaln
    return gammaln(n + 1) - (gammaln(k + 1) + gammaln(n - k + 1))


def _ln_F(k, n):
    '''
    Log of F(k, n) as appearing in Conlisk et al 2007 Am Nat, also F(a, n) in 
    Conlisk et al 2007 Ecol Mono. This form of F from Ecol Mono paper.
    '''
    gammaln = scipy.special.gammaln

    return gammaln(k + n) - (gammaln(k) + gammaln(n + 1))
