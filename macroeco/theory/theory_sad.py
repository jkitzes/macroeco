#!/usr/bin/python

'''
Calculate pmf and likelihood of species abundance distributions.

Distributions
-------------
- `lgsr_pmf` -- Fisher's log series (Fisher et al. 1943)
- `neg_binom_pmf` -- Negative Binomial
- `geo_pmf` -- Geometric Distribution
- `plognorm_pmf` -- Poisson lognormal (Bulmer 1974)
- `trun_plognorm_pmf` -- Truncated poisson log-normal (Bulmer 1974)
- `mete_logser_pmf` -- METE log series (Harte 2011)
- `mete_logser_approx_pmf` -- METE log series using approximation (Harte 2011)

Misc Functions
--------------
- `make_rank_abund` -- convert any SAD pmf into a rank abundance curve
- `plognorm_MLE` -- Get MLE estimates for poisson lognormal
- `macroeco_pmf`
- `get_sad_cdf`
- `nll`
- `distr_parameters`
- `fit_neg_binom`

References
----------
Bulmer, M. G. 1974. On fitting the poisson lognormal distribution to species
abundance data. Biometrics, 30:101-110.

Fisher, R. A., Corbet, A. S., and C. B. Williams. 1943. The relation between
The number of species and the number of individuals in a random sample
of an animal population. Journal of Animal Ecology, 12:42-58.

Harte, J. 2011. Maximum Entropy and Ecology: A Theory of Abundance,
Distribution, and Energetics. Oxford University Press.

Hubbell, S. P. 2001. The unified theory of biodiversity and biogeography. 
Monographs in Population Biology, 32,1:375.

Magurran, A. E. 1988. Ecological Diversity and Its Measuremnt. Princeton
University Press.

May, R. M. 1975. Patterns of species abundance and diversity. In Ecology and
Evolution of Communities (eds M. L. Cody and J. M. Diamond), Harvard University
Press.

Motomura, L. 1932. A statistical treatment of associations. Japan Journal of
Zoology, 44:379-383.

'''
from __future__ import division
import numpy as np
import scipy.stats as stats
import scipy.optimize 
import scipy.special
import math as m
import scipy.integrate as integrate
import sys
#NOTE: Assertion statements needed!

class RootError(Exception):
    '''Error if no root or multiple roots exist for the equation generated
    for specified values of S and N'''

    def __init__(self, value=None):
        Exception.__init__(self)
        self.value = value
    def __str__(self):
        return '%s' % self.value

def lgsr_pmf(n, S, N, param_out=False):
    '''
    Fisher's log series pmf (Fisher et al. 1943, Hubbel 2001)

    Parameters
    ----------
    n : int, float or array-like object
        Abundances at which to calculate the pmf
    S : int
        Total number of species in landscape
    N : int
        Total number of indviduals in landscape
    abundances : array like object
        Array-like object containing sad

    Returns
    -------
    : ndarray (1D)
        Returns array with pmf for values for the given values n. If 
        param_out = True, returns the array as well as the parameter estimates.


    Notes
    -----
    Multiplying the pmf by S yields the predicted number of species
    with a given abundance.
    '''
    assert S < N, "S must be less than N"
    assert S > 1, "S must be greater than 1"
    assert N > 0, "N must be greater than 0"
    if type(n) is int or type(n) is float:
        n = np.array([n])
    else:
        n = np.array(n)

    start = -2
    stop = 1 - 1e-10
    
    eq = lambda x, S, N: (((N/x) - N) * (-(np.log(1 - x)))) - S
    
    x = scipy.optimize.brentq(eq, start, stop, args=(S,N), disp=True)
    pmf = stats.logser.pmf(n, x)

    
    if param_out == True:
        return (pmf, {'x' : x})
    else:
        return pmf

def neg_binom_pmf(n, S, N, k, param_out=False):
    '''
    Negative binomial distribution 
    
    Parameters
    ----------
    n : int, float or array-like object
        Abundances at which to calculate the pmf
    S : int
        Total number of species in landscape
    N : int
        Total number of individuals in landscape
    k : int
        Aggregation parameter

    Returns
    -------
     : ndarray (1D)
        Returns array with pmf for values for the given values n. If 
        param_out = True, returns the array as well as the parameter estimates.

    '''
    assert S < N, "S must be less than N"
    assert S > 1, "S must be greater than 1"
    assert N > 0, "N must be greater than 0"
    if type(n) is int or type(n) is float:
        n = np.array([n])
    else:
        n = np.array(n)

    mu = float(N) / S
    p = 1 / (mu / k + 1)  # See Bolker book Chapt 4
    pmf = stats.nbinom.pmf(n, k, p)

    if param_out == True:
        return (pmf, {'k' : k, 'p' : p})
    else:
        return pmf

def geo_pmf(n, S, N, param_out=False):
    '''
    Geometric distribution. Using neg_binom_pmf with k=1 as a wrapper
    
    Parameters
    ----------
    n : int, float or array-like object
        Abundances at which to calculate the pmf
    S : int
        Total number of species in landscape
    N : int
        Total number of individuals in landscape
        
    Returns
    -------
     : ndarray (1D)
        Returns array with pmf for values for the given values n. If 
        param_out = True, returns the array as well as the parameter estimates.

    '''
    if param_out == True:
        pmf, params =  neg_binom_pmf(n, S, N, 1, param_out=param_out)
        p = params['p']
        return (pmf, {'p' : p})
    else:
        return neg_binom_pmf(n, S, N, 1, param_out=param_out)

def geo_series_rank_abund(S, N, k):
    '''
    geo_series_rank_abund(S, N, k)

    Given S, N and k, this function returns the predicted rank-abundance
    distribution for a geometric series distribution (Motomura 1932 and 
    Magurran 1988).

    Parameters
    ----------
    S : int
        Total number of species in landscape
    N : int
        Total number of individuals in landscape
    k : float
        The fraction of resources that each species acquires. Range is (0, 1].

    Returns
    -------
    : np.ndarray
        An array containing the the predicted SAD with element 0 containing the
        species with the most individuals and element S - 1 containing the
        species with the least individuals.
    
    '''

    assert S < N, "S must be less than N"
    assert S > 1, "S must be greater than 1"
    assert N > 0, "N must be greater than 0"
    assert k > 0 and k <= 1, "k must be between on the interval (0, 1]"

    C = (1 - (1 - k )** S) ** - 1
    pred_rank_abund = N * C * k * (1 - k) ** (np.arange(1, S + 1) - 1)
    return pred_rank_abund

def broken_stick_pmf(n, S, N, param_out=False):
    '''
    broken_stick_pmf(n, S, N, param_out=False)
    
    McArthur's broken stick species abundance distribution (May 1975)

    Parameters
    ----------
    n : int, float or array-like object
        Abundances at which to calculate the pmf
    S : int
        Total number of species in landscape
    N : int
        Total number of individuals in landscape

    Returns
    -------
     : ndarray (1D) or tuple
        Returns array with pmf for values for the given values n. If 
        param_out = True, returns a tuple with the pmf as an array as well as
        a dictionary of the parameter estimates.
    
    '''

    assert S < N, "S must be less than N"
    assert S > 1, "S must be greater than 1"
    assert N > 0, "N must be greater than 0"
    if type(n) is int or type(n) is float:
        n = np.array([n])
    else:
        n = np.array(n)

    eq = lambda x: ((S - 1) / N) * (1 - (x / N)) ** (S - 2)
    pmf = eq(n)

    if param_out == True:
        return (pmf, {'S' : S})
    else:
        return pmf

def broken_stick_rank_abund(S, N):
    '''
    broken_stick_rank_abund(S, N)

    This distribution returns the predicted rank-abundance distribution for
    McArthur's broken-stick distribution (May 1975).

    Parameters
    ----------
    S : int
        Total number of species in landscape
    N : int
        Total number of individuals in landscape

    Returns
    -------
    : np.ndarray
        An array containing the the predicted SAD with element 0 containing the
        species with the most individuals and element S - 1 containing the
        species with the least individuals.

    '''

    assert S < N, "S must be less than N"
    assert S > 1, "S must be greater than 1"
    assert N > 0, "N must be greater than 0"

    pred_rank_abund = np.empty(S)
    for i in xrange(S):
        n = np.arange(i + 1, S + 1) 
        pred_rank_abund[i] = (N / S) * sum(1 / n)
    return pred_rank_abund


def plognorm_pmf(ab, mu, sigma, param_out=False):
    '''
    Poisson log-normal pmf (Bulmer 1974)

    Parameters
    ----------
    ab : int, float or array-like object
        Abundances at which to calculate the pmf 
    mu : float
        the mu parameter of the poisson log normal
    sigma : float
        the sigma parameter of the poisson log normal
        
     Returns
    -------
     : ndarray (1D)
        Returns array with pmf for values for the given values ab. If 
        param_out = True, returns the array as well as the parameter estimates.

    Notes
    -----
    This fuction was adopted directly from the VGAM package in R by Mark
    Wilber. The VGAM R package was adopted directly from Bulmer (1974).

    '''
    if type(ab) is int or type(ab) is float:
        ab = np.array([ab])
    else:
        ab = np.array(ab)
    n_unq = np.unique(ab)
    
    if sigma**2 <= 0 or mu <= 0: #Parameters could be negative
        pmf = np.repeat(1e-120, len(n_unq))
    else:
        eq = lambda t, x: np.exp(t * x - np.exp(t) - 0.5*((t - mu) / sigma)**2)
        pmf = np.empty(len(n_unq), dtype=np.float)
        for i, n in enumerate(n_unq):
            if n <= 170:
                integral = integrate.quad(eq, -np.inf, np.inf, args=(n))[0]
                norm = np.exp((-0.5 * m.log(2 * m.pi * sigma**2) -\
                                                m.lgamma(n + 1)))
                pmf[i] = norm * integral
            else:
                z = (m.log(n) - mu) / sigma
                pmf[i] = (1 + (z**2 + m.log(n) - mu - 1) / (2 * n * sigma**2))\
                         * np.exp(-0.5 * z**2) / (m.sqrt(2 * m.pi) * sigma * n)   

    #Only calculated unique abundances to save computational time.
    #Get full pmf again
    pmf_full = np.empty(len(ab))
    for i, n in enumerate(ab):
        index = np.where(n_unq == n)[0][0]
        pmf_full[i] = pmf[index]
    pmf = pmf_full
     
    if param_out == True:
        #NOTE: How many parameters should I be returning, 1 or 2?
        return (pmf, {'mu' : mu, 'sigma' : sigma})
    else:
        return pmf

def trun_plognorm_pmf(ab, mu, sigma, param_out=False):
    '''
    Truncated Poisson log-normal (Bulmer 1974)

    Parameters
    ----------
    ab : int, float or array-like object
        Abundances at which to calculate the pmf
    mu : float
        the mu parameter of the poisson log normal
    sigma : float
        the sigma parameter of the poisson log normal
        
    Returns
    -------
      : ndarray (1D)
        Returns array with pmf for values for the given values ab. If 
        param_out = True, returns the array as well as the parameter estimates.

    Notes:  This function was adopted from both Bulmer (1974) and Ethan White's
    code from weecology.  Truncating the plognormal changes the mean of the 
    distribution.  Need to make this change!'''
    
    if param_out == True:
        untr_pmf, mn, sd = plognorm_pmf(ab, mu, sigma, param_out=param_out)
        pmf0 = plognorm_pmf(0, mu, sigma)
        tr_pmf = (untr_pmf / (1 - pmf0))#Truncating based on Bulmer equation A1
        return (tr_pmf, {'mu' : mn, 'sigma' : sd})
    else:
        untr_pmf = plognorm_pmf(ab, mu, sigma)
        pmf0 = plognorm_pmf(0, mu, sigma)
        tr_pmf = (untr_pmf / (1 - pmf0))
        return tr_pmf

def lognorm_pmf(ab, mu, sigma, param_out=False):
    '''
    Lognormal pmf

    Parameters
    ----------
    ab : int, float or array-like object
        Abundances at which to calculate the pmf 
    mu : float
        the mu parameter of the poisson log normal
    sigma : float
        the sigma parameter of the poisson log normal

    Returns
    -------
     : ndarray (1D)
        Returns array with pmf for values for the given values ab. If 
        param_out = True, returns the array as well as the parameter estimates.

    Notes
    -----
    scipy.stats.lognorm is coded very poorly and the docstring is not 
    helpful so we coded our own lognormal distribution

    '''
    if type(ab) is int or type(ab) is float:
        ab = np.array([ab])
    else:
        ab = np.array(ab)
    pmf = stats.norm.pdf(np.log(ab), loc=mu, scale=sigma) / ab
    if param_out == True:
        return (pmf, {'mu' : mu, 'sigma' : sigma})
    else:
        return pmf

def trun_lognormal(ab, mu, sigma, param_out=False):
    '''
    Zero-Truncated Lognormal pmf

    Parameters
    ----------
    ab : int, float or array-like object
        Abundances at which to calculate the pmf 
    mu : float
        the mu parameter of the poisson log normal
    sigma : float
        the sigma parameter of the poisson log normal

    Returns
    -------
     : ndarray (1D)
        Returns array with pmf for values for the given values ab. If 
        param_out = True, returns the array as well as the parameter estimates.

    Notes
    -----
    Truncating shifts the mean. Need to take this into account
    '''

    if param_out == True:
        untr_pmf, mn, sd = lognorm_pmf(ab, mu, sigma, param_out=param_out)
        pmf0 = lognorm_pmf(0, mu, sigma)
        tr_pmf = (untr_pmf / (1 - pmf0))
        return (tr_pmf, {'mu' : mn, 'sigma' : sd})
    else:
        untr_pmf = lognorm_pmf(ab, mu, sigma)
        pmf0 = lognorm_pmf(0, mu, sigma)
        tr_pmf = (untr_pmf / (1 - pmf0))
        return tr_pmf

def plognorm_MLE(ab, trun=True):
    '''
    Maximum likelihood Estimates for Poisson log normal

    Parameter
    ---------
    ab : A list, np.array, or tuple of the abundance of each
                 species in the plot.  len(abundance) should equal 
                 the total number of species in the plot and sum(abundances)
                 should equal the total number of individuals in the plot
    trun : bool
        If true, calculates the MLE's for the truncated poisson lognormal.
        If false, calulates the MLE's for the untruncated poisson lognormal.

    Returns
    -------
    mu, var : float
        The ML estimates of mu and var

    Notes
    -----
    This function was adapted from Ethan White's pln_solver function in 
    weecology. 

    '''
    assert type(ab) == list or type(ab) == tuple \
           or type(ab) == np.ndarray, "Invalid parameter type"

    assert len(ab) >= 1, "len(ab) must be greater than or equal to 1"
    ab = np.array(ab)
    mu0 = np.mean(np.log(ab))
    var0 = np.var(np.log(ab), ddof=1)
    def pln_func(x):
        if trun == True:
            return -sum(np.log(trun_plognorm_pmf(ab, x[0], x[1])))
        else:
            return -sum(np.log(plognorm_pmf(ab, x[0], x[1])))
    mu, var = scipy.optimize.fmin(pln_func, x0 = [mu0, var0], disp=0)
    return mu, var


def mete_lgsr_pmf(n, S, N, param_out=False):
    '''
    mete_logsr_pmf(n, S, N, param_out=False)

    Truncated log series pmf (Harte 2011)

    Parameters:
    -----------
    n : int, float or array-like object
        Abundances at which to calculate the pmf
    S : int
        Total number of species in landscape
    N : int
        Total number of individuals in landscape
        
    Returns
    -------
     : ndarray (1D)
        Returns array with pmf for values for the given values n. If 
        param_out = True, returns the array as well as the parameter estimates.

    Notes
    -----
    This function uses the truncated log series as described in Harte 2011
    eq (7.32).  The equation used in this function to solve for the Lagrange 
    multiplier is equation (7.27) as described in Harte 2011. 
    
    Also note that realistic values of x where x = e**-(beta) (see Harte 2011) 
    arein the range (1/e, 1) (exclusive). Therefore, the start and stop 
    parameters for the brentq procedure are close to these values. However, x 
    can occasionally be greater than one so the maximum stop value of the 
    brentq optimizer is 2.
    '''
    assert S < N, "S must be less than N"
    assert S > 1, "S must be greater than 1"
    assert N > 0, "N must be greater than 0"
    if type(n) is int or type(n) is float:
        n = np.array([n])
    else:
        n = np.array(n)
    start = 0.3
    stop = 2
    
       
    k = np.linspace(1, N, num=N)
    eq = lambda x: sum(x ** k / float(N) * S) -  sum((x ** k) / k)
    x = scipy.optimize.brentq(eq, start, min((sys.float_info[0] / S)**\
                              (1/float(N)), stop), disp=True)
    #Taken from Ethan White's trun_logser_pmf
    norm = sum(x ** k / k)        
    pmf = (x ** n / n) / norm
    
    if param_out == True:
        return (pmf, {'beta' : -np.log(x)}) 
    else:
        return pmf

def mete_lgsr_approx_pmf(n, S, N, param_out=False, root=2):
    '''
    mete_lgsr_approx_pmf(n, S, N, param_out=False, root=2)

    Truncated log series using approximation (7.30) and (7.32) in Harte 2011

    Parameters
    ----------
    n : int, float or array-like object
        Abundances at which to calculate the pmf
    S : int
        Total number of species in landscape
    N : int
        Total number of individuals in landscape
    root: int (optional)
        1 or 2.  Specifies which root to use for pmf calculations        
    Returns
    -------
     : ndarray (1D)
        Returns array with pmf for values for the given values n. If 
        param_out = True, returns the array as well as the parameter estimates.
        
     Notes:
    ------
    This function uses the truncated log series as described in Harte 2011
    eq (7.32).  The equation used in this function to solve for the Lagrange 
    multiplier is equation (7.30) as described in Harte 2011.     
       
    Also note that realistic values of x where x = e^-(beta) (see Harte 2011) 
    are in the range (1/e, 1) (exclusive). Therefore, the start and stop 
    parameters for the brentq optimizer have been chosen near these values.
    '''

    assert S < N, "S must be less than N"
    assert S > 1, "S must be greater than 1"
    assert N > 0, "N must be greater than 0"
    if type(n) is int or type(n) is float:
        n = np.array([n])
    else:
        n = np.array(n)
    start = 0.3
    stop = 1 - 1e-10
    
    eq = lambda x: ((-m.log(x))*(m.log(-1/(m.log(x))))) - (float(S)/N)
    try:
        x = scipy.optimize.brentq(eq, start, stop, disp=True)
    except ValueError:
        values = np.linspace(start, stop, num=1000)
        solns = []
        for i in values:
            solns.append(eq(i))#Find maximum
        ymax = np.max(np.array(solns))
        xmax = values[solns.index(ymax)]
        if ymax > 0:
            print "Warning: More than one solution."
            if root == 1:
                x = scipy.optimize.brentq(eq, start, xmax, disp=True)
            if root == 2:
                x = scipy.optimize.brentq(eq, xmax, stop, disp=True)
                       
        if ymax < 0:
            raise RootError('No solution to constraint equation with given '\
                            + 'values of S and N') 
    g = -1/m.log(x)
    pmf = (1/m.log(g)) * ((x**n)/n) 
    
    if param_out == True:
        return (pmf, {'beta' : -np.log(x)})
    else:
        return pmf

def make_rank_abund(sad_pmf, S):
    '''
    Convert any SAD pmf into a rank abundance curve for S species using 
    cumulative distribution function.
 
    Parameters
    ----------
    sad_pmf : ndarray
        Probability of observing a species from 1 to length sad_pmf individs
    S : int
        Total number of species in landscape

    Returns
    -------
    : ndarray (1D)
        If summary is False, returns array with pmf. If summary is True,
        returns the summed log likelihood of all values in n.

    Notes
    -----
    Function actually implements (philosophically) a step quantile function. 
    '''

    S_points = np.arange(1/(2*S), 1, 1/S)
    S_abunds = np.zeros(S)

    sad_pmf_w_zero = np.array([0] + list(sad_pmf)) # Add 0 to start of pmf
    cum_sad_pmf_w_zero = np.cumsum(sad_pmf_w_zero)
    
    for cutoff in cum_sad_pmf_w_zero:
        greater_thans = (S_points >= cutoff)
        S_abunds[greater_thans] += 1

        if not greater_thans.any():  # If no greater thans, done with all S
            break
    
    return S_abunds

def get_sad_cdf(S, N, distribution, sad=[]):
    '''
    Gets predicted cdf for a given distribution based on the sad

    CDF for METE logseries distribution

    Parameters
    ----------
    S : int
        Total number of species in landscape
    N : int
        Total number of individuals in landscape
    distribution : string
        The predicted distribution:
        'mete' - METE
        'mete_approx' - METE with approximation
        'plognorm' - Poisson lognormal
        'trun_plognorm' - Truncated poisson lognormal
        'neg_binom' - Negative binomial
        'geo' - Geometric distribution
        'lgsr' - Fisher's log series

        
    sad : array like object
        SAD can be empty if distribution is not plognorm, trun_plognorm
        or neg binom

    Returns:
    : 1D structured array
        Field names in the structured array are 'n' (number of individuals) and 'cdf'


    '''
    assert distribution == 'mete' or\
           distribution == 'mete_approx' or\
           distribution == 'neg_binom' or\
           distribution == 'geo' or\
           distribution == 'lgsr' or\
           distribution == 'trun_plognorm' or\
           distribution == 'plognorm', "Do not recognize %s distribution" % (distribution)

    #These three distributions require the full sad.
    if distribution == 'plognorm' or\
       distribution == 'trun_plognorm' or\
       distribution == 'neg_binom':
        assert len(sad) == S, "Length of SAD must equal S"
        assert np.sum(sad) == N, "Sum of SAD must equal N"
    sad = np.array(sad)

    
    if distribution == 'mete':
        pmf = mete_lgsr_pmf(S, N, abundances=sad, pmf_ret=True)
    if distribution == 'mete_approx':
        pmf = mete_lgsr_approx_pmf(S, N, abundances=sad, pmf_ret=True)
    if distribution == 'neg_binom':
        mlek = fit_neg_binom_pmf(sad)
        pmf = neg_binom_pmf(S, N, mlek, abundances=sad, pmf_ret=True)
    if distribution == 'geo':
        pmf = geo_pmf(S, N, abundances=sad, pmf_ret=True)
    if distribution == 'trun_plognorm':
        mu, var = plognorm_MLE(sad, trun=True)
        pmf = trun_plognorm_pmf(mu, var, sad, pmf_ret=True)
    if distribution == 'plognorm':
        mu, var = plognorm_MLE(sad, trun=False)
        pmf = plognorm_pmf(mu, var, sad, pmf_ret=True)
    if distribution == 'lgsr':
        pmf = lgsr_pmf(S, N, abundances=sad, pmf_ret=True)

    cdf = np.cumsum(pmf)
    cdf_struct = np.empty(len(cdf), dtype=[('cdf', np.float), ('n', np.int)])
    cdf_struct['cdf'] = cdf
    cdf_struct['n'] = np.arange(1, N + 1)
    return cdf_struct


def nll(sad, distribution):
    '''
    Calculates the negative log-likelihood for different distributions

    Parameters
    ----------
    sad : ndarray
        An array-like object with species abundances

    distribution : string
        Specifies the distribution for which to get the negative log-likelihood
        'mete' - mete distribution
        'mete_approx' - mete distribution with approximation
        'neg_binom' - negative binomial distribution
        'geo' - geometric distribution
        'plognorm' - poisson log-normal distribution
        'trun_plognorm' - truncated poisson log-normal distribution
        'lgsr' - Fisher's log series 

    Returns
    -------
    : float
        The negative log-likelihood

    '''

    assert distribution == 'mete' or\
           distribution == 'mete_approx' or\
           distribution == 'neg_binom' or\
           distribution == 'geo'or\
           distribution == 'lgsr' or\
           distribution == 'trun_plognorm' or\
           distribution == 'plognorm', "Do not recognize %s distribution" % (distribution)

    if distribution == 'mete':
        pmf = mete_lgsr_pmf(sad, len(sad), sum(sad))
    if distribution == 'mete_approx':
        pmf = mete_lgsr_approx_pmf(len(sad), sum(sad), abundances=sad)
    if distribution == 'neg_binom':
        mlek = fit_neg_binom_pmf(sad)
        pmf = neg_binom_pmf(len(sad), sum(sad), mlek, abundances=sad)
    if distribution == 'geo':
        pmf = geo_pmf(len(sad), sum(sad), abundances=sad)
    if distribution == 'trun_plognorm':
        mu, var = plognorm_MLE(sad, trun=True)
        pmf = trun_plognorm_pmf(mu, var, sad)
    if distribution == 'plognorm':
        mu, var = plognorm_MLE(sad, trun=False)
        pmf = plognorm_pmf(mu, var, sad)
    if distribution == 'lgsr':
        pmf = lgsr_pmf(len(sad), sum(sad), abundances=sad)

    return -sum(np.log(pmf))

def macroeco_pmf(S, N, distribution, sad=[]):
    '''
    This function returns a the full pmf described by a distribution

    Parameters
    ----------
    S : int
        Total number of species in landscape
    N : int
        Total number of individuals in landscape
    distribution : string
        The predicted distribution:
        'mete' - METE
        'mete_approx' - METE with approximation
        'plognorm' - Poisson lognormal
        'trun_plognorm' - Truncated poisson lognormal
        'neg_binom' - Negative binomial
        'geo' - Geometric
        'lgsr' - Fisher's log series

        
    sad : array like object
        SAD can be empty if distribution is not plognorm, trun_plognorm, or neg_binom

    Returns
    -------
    :ndarray
        The full pmf/pdf with support [1,N]



    '''
    assert distribution == 'mete' or\
           distribution == 'mete_approx' or\
           distribution == 'neg_binom' or\
           distribution == 'geo' or\
           distribution == 'lgsr' or\
           distribution == 'trun_plognorm' or\
           distribution == 'plognorm', "Do not recognize %s distribution" % (distribution)

    #These three distributions require the full sad.
    if distribution == 'plognorm' or\
       distribution == 'trun_plognorm' or\
       distribution == 'neg_binom':
        assert len(sad) == S, "Length of SAD must equal S"
        assert np.sum(sad) == N, "Sum of SAD must equal N"
    sad = np.array(sad)

    
    if distribution == 'mete':
        pmf = mete_lgsr_pmf(S, N, abundances=sad, pmf_ret=True)
    if distribution == 'mete_approx':
        pmf = mete_lgsr_approx_pmf(S, N, abundances=sad, pmf_ret=True)
    if distribution == 'neg_binom':
        mlek = fit_neg_binom_pmf(sad)
        pmf = neg_binom_pmf(S, N, mlek, abundances=sad, pmf_ret=True)
    if distribution == 'geo':
        pmf = geo_pmf(S, N, abundances=sad, pmf_ret=True)
    if distribution == 'trun_plognorm':
        mu, var = plognorm_MLE(sad, trun=True)
        pmf = trun_plognorm_pmf(mu, var, sad, pmf_ret=True)
    if distribution == 'plognorm':
        mu, var = plognorm_MLE(sad, trun=False)
        pmf = plognorm_pmf(mu, var, sad, pmf_ret=True)
    if distribution == 'lgsr':
        pmf = lgsr_pmf(S, N, abundances=sad, pmf_ret=True)

    return pmf

def distr_parameters(S, N, distribution, sad=[]):
    '''
    This function returns the best fit parameters for the given distribution

    Parameters
    ----------
    S : int
        Total number of species in landscape
    N : int
        Total number of individuals in landscape
    distribution : string
        The predicted distribution:
        'mete' - METE
        'mete_approx' - METE with approximation
        'plognorm' - Poisson lognormal
        'trun_plognorm' - Truncated poisson lognormal
        'neg_binom' - Negative binomial
        'geo' - Geometric
        'lgsr' - Fisher's log series

        
    sad : array like object
        SAD can be empty if distribution is not plognorm, trun_plognorm, or neg_binom

    Returns
    -------
    :dictionary 
        Dictionary of parameters for given distribution


    '''
    assert distribution == 'mete' or\
           distribution == 'mete_approx' or\
           distribution == 'neg_binom' or\
           distribution == 'geo' or\
           distribution == 'lgsr' or\
           distribution == 'trun_plognorm' or\
           distribution == 'plognorm', "Do not recognize %s distribution" % (distribution)

    #These three distributions require the full sad.
    if distribution == 'plognorm' or\
       distribution == 'trun_plognorm' or\
       distribution == 'neg_binom':
        assert len(sad) == S, "Length of SAD must equal S"
        assert np.sum(sad) == N, "Sum of SAD must equal N"
    sad = np.array(sad)
    
    parameters = {}
    if distribution == 'mete':
        beta = -np.log(mete_lgsr_pmf(S, N, abundances=sad, \
                            pmf_ret=True, param_out=True)[1])
        parameters['beta'] = beta
    if distribution == 'mete_approx':
        beta = -np.log(mete_lgsr_approx_pmf(S, N, abundances=sad, \
                            pmf_ret=True, param_out=True)[1])
        parameters['beta'] = beta
    if distribution == 'neg_binom':
        mlek = fit_neg_binom_pmf(sad)
        k, p = neg_binom_pmf(S, N, mlek, abundances=sad, pmf_ret=True, param_out=True)[1:]
        parameters['k'] = k
        parameters['p'] = p
    if distribution == 'geo':
        p = geo_pmf(S, N, abundances=sad, pmf_ret=True, param_out=True)[1]
        parameters['p'] = p
    if distribution == 'trun_plognorm':
        mu, var = plognorm_MLE(sad, trun=True)
        parameters['mu'] = mu
        parameters['var'] = var
    if distribution == 'plognorm':
        mu, var = plognorm_MLE(sad, trun=False)
        parameters['mu'] = mu
        parameters['var'] = var
    if distribution == 'lgsr':
        x = lgsr_pmf(S, N, abundances=sad, pmf_ret=True, param_out=True)[1]
        parameters['x'] = x

    return parameters

def fit_neg_binom_pmf(sad, guess_for_k=1):
    '''
    Function fits a negative binomial to the sad and returns the MLE for k 

    Parameters
    ----------
    sad : ndarray
        Array like object containing a species abundance distribution
    
    guess_for_k : float
        A default parameter for the approximate k given the data

    Returns
    -------
    : float
        the mle for k

    Notes
    -----
    scipy.stats.nbinom.fit does not exist. That is why the MLE estimator is
    hardcoded.

    '''

    #NOTE: Need to check for convergence
    def nll_nb(k, sad):
        return -sum(np.log(neg_binom_pmf(sad, len(sad), np.sum(sad), k)))
    mlek = scipy.optimize.fmin(nll_nb, np.array([guess_for_k]), args=(sad,),\
                                                                    disp=0)[0]
    return mlek






        






