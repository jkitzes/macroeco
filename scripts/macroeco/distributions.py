#!/usr/bin/python

'''
Macroecological distributions.

Distributions
-------------
- `lgsr` -- Fisher's log series (Fisher et al. 1943)
- `trun_neg_binom` -- Truncated negative binomial
- `geo_series' -- Geometric series distribution (Motomura 1932)
- `broken_stick` -- McArthur's broken stick distribution (May 1975)
- `plognorm` -- Poisson log-normal (Bulmer 1974)
- `trun_plognorm` -- Truncated poisson log-normal (Bulmer 1974)
- `lognorm_pmf` -- Lognormal distribution
- `canonical_lognormal_pmf` -- Preston's canonical lognormal parameterized by
  May (1975)
- `sugihara` -- Sugihara's sequential breakage model (Sugihara 1980)
- `mete_lgsr` -- METE log series (Harte 2011)
- `mete_lgsr_approx` -- METE log series using approximation (Harte 2011)
- `binm` - Binomial distribution (Random Placement Model)
- `pois` - Poisson distribution
- `nbd` - Negative binomial distribution
- `fnbd` - Finite negative binomial (Zillio and He 2010)
- `cnbd` - Conditional negative binomial distribution (Conlisk et al 2007b)
- `geo` - Geometric distribution
- `fgeo` - Finite geometric distribution (Zillio and He 2010)
- `tgeo` - Truncated geometric distribution (Harte et al. 2008)
- `mete_sar` - METE sar functions (Harte 2011)
- `SAR` - General non-METE sar functions

Misc Functions
--------------
- `make_rank_abund` -- convert any SAD pmf into a rank abundance curve
- `_ln_choose`
- `_downscale_sar_`
- `_upscale_sar_`
- `_generate_areas_`

References
----------
Bulmer, M. G. 1974. On fitting the poisson lognormal distribution to species
abundance data. Biometrics, 30:101-110.

Conlisk E, Bloxham M, Conlisk J, Enquist B, Harte J (2007a) A new class of 
models of spatial distribution. Ecological Monographs 77:269-284.

Conlisk E, Conlisk J, Harte J (2007b) The impossibility of estimating a 
negative binomial clustering parameter from presence-absence data: a comment on 
He and Gaston. The American Naturalist 170:651-654.

Fisher, R. A., Corbet, A. S., and C. B. Williams. 1943. The relation between
The number of species and the number of individuals in a random sample
of an animal population. Journal of Animal Ecology, 12:42-58.

Harte, J. 2011. Maximum Entropy and Ecology: A Theory of Abundance,
Distribution, and Energetics. Oxford University Press.

Harte J, Conlisk E, Ostling A, Green JL, Smith AB (2005) A theory of spatial 
structure in ecological communities at multiple spatial scales. Ecological 
Monographs 75:179-197.

Hubbell, S. P. 2001. The unified theory of biodiversity and biogeography. 
Monographs in Population Biology, 32,1:375.

Magurran, A. E. 1988. Ecological Diversity and Its Measuremnt. Princeton
University Press.

May, R. M. 1975. Patterns of species abundance and diversity. In Ecology and
Evolution of Communities (eds M. L. Cody and J. M. Diamond), Harvard University
Press.

Motomura, L. 1932. A statistical treatment of associations. Japan Journal of
Zoology, 44:379-383.

Sugihara, G. 1980. Minimal community structure: an explanation of species
abundance patterns. The American Naturalist, 116:770-787.

Zillio T, He F (2010) Modeling spatial aggregation of finite populations. 
Ecology 91:3698-3706.

'''
from __future__ import division
import numpy as np
import scipy.stats as stats
import scipy.optimize 
import scipy.special
import math as m
import scipy.integrate as integrate
import sys

__author__ = "Justin Kitzes and Mark Wilber"
__copyright__ = "Copyright 2012, Regents of the University of California"
__credits__ = ["John Harte"]
__license__ = None
__version__ = "0.1"
__maintainer__ = "Justin Kitzes and Mark Wilber"
__email__ = "jkitzes@berkeley.edu"
__status__ = "Development"

#TODO: Add truncated log-normal?
#TODO: Which distributions should have param_out?
#TODO: Check test_distributions and account for changes 

class RootError(Exception):
    '''Error if no root or multiple roots exist for the equation generated
    for specified values of S and N'''

    def __init__(self, value=None):
        Exception.__init__(self)
        self.value = value
    def __str__(self):
        return '%s' % self.value

class DownscaleError(Exception):
    '''Catch downscale errors'''
    def __init__(self, value=None):
        Exception.__init__(self)
        self.value = value
    def __str__(self):
        return '%s' % self.value

class Distribution(object):

    def __init__(self, **kwargs):
        '''
        Generic constructor

        **kwargs : keyword parameters for distribution

        '''
        self.params = kwargs
    
    #Better option than lower bound?
    def cdf(self, n, lower_bound=1):
        '''
        Cumulative distribution function.  Determined by summing pmf

        Parameters
        ----------
        n : int, float or array-like object
            Values at which to calculate the cdf

        Returns
        -------
        : ndarray (1D)
            Returns array with cdf values for the given values of n.

        '''
        try:
            len(n); n = np.array(n)
        except:
            n = np.array([n])
        max_n = np.max(n)
        cdf = np.cumsum(self.pmf(np.arange(lower_bound, max_n + 1)))
        return np.array([cdf[x - lower_bound] for x in n])

    def rad(self):
        '''
        Rank abundance distribution. Calculated using pmf

        Parameters
        ----------
        None

        Returns
        -------
        : ndarray (1D)
            Returns and array of length S with the expected abundances given
            pmf

        '''
        N = self.params.get('N', None)
        S = self.params.get('S', None)
        assert N != None, "N parameter not given"
        assert S != None, "S parameter not given"
        full_pmf = self.pmf(np.arange(1, N + 1))
        return make_rank_abund(full_pmf, S)

    def fit(self, data, sad=True):
        '''
        Generic fit to sad or ssads

        Parameters
        ----------
        data : array-like object
            Data used to calculate fit
        sad : bool
            If True, calculate S and N parameter, else only calculate N.

        Return
        ------
        None
        '''
        #Might not be the best check here
        if sad:
            self.params['N'] = np.sum(data)
            self.params['S'] = len(data)
        else:
            self.params['N'] = sum(data)

    def nll(self, n):
        '''
        Calcuates the negative log-likelihood for a given distribution

        Parameters
        ----------
        n : array-like object
            Values for which to calculate nll

        Returns
        : float
            Negative log-likelihood for given values

        '''
        return -sum(np.log(self.pmf(n)))

class lgsr(Distribution):
    '''
    Fisher's log series distribution (Fisher et al. 1943, Hubbel 2001).

    Methods
    -------
    pmf(n, param_out); S and N parameters passed into __init__
        Probability mass function
    cdf(n); S and N parameters passed into __init__
        Cumulative distribution function
    rad(); S and N parameters passed into __init__
        Rank abundance distribution

    Parameters
    ----------
    S : int
        Total number of species in landscape
    N : int
        Total number of indviduals in landscape

    Notes
    -----
    S and N are passed into the constructor (__init__) as keyword arguments.
    Example: lgsr(S=34, N=345)

    '''
    def pmf(self, n, param_out=False):
        '''
        Probability mass function of Fisher log series generated with S and N.

        Parameters
        ----------
        n : int, float or array-like object
            Abundances at which to calculate the pmf
        
        Passed into __init__
        ---------------------
        S : int
            Total number of species in landscape
        N : int
            Total number of indviduals in landscape

        Returns
        -------
        : ndarray (1D)
            Returns array with pmf for values for the given values n. If 
            param_out = True, returns the array as well as the parameter 
            estimates.

        Notes
        -----
        Multiplying the pmf by S yields the predicted number of species
        with a given abundance.

        '''
        S = self.params.get('S', None)
        N = self.params.get('N', None) 
        assert S != None, "S parameter not given"
        assert N != None, "N parameter not given"
        assert S < N, "S must be less than N"
        assert S > 1, "S must be greater than 1"
        assert N > 0, "N must be greater than 0"
        try:
            len(n); n = np.array(n)
        except:
            n = np.array([n])
        assert np.max(n) <= N, "Maximum n cannot be greater than N"

        start = -2
        stop = 1 - 1e-10
    
        eq = lambda x, S, N: (((N/x) - N) * (-(np.log(1 - x)))) - S
    
        x = scipy.optimize.brentq(eq, start, stop, args=(S,N), disp=True)
        pmf = stats.logser.pmf(n, x)

    
        if param_out == True:
            return (pmf, {'x' : x})
        else:
            return pmf

class trun_neg_binom(Distribution):
    '''
    Zero-truncated negative binomial distribution .  
    Because this distribution has infinite support the fit the cdf fron 1:N 
    will be less than one.
    Methods
    -------
    pmf(n, param_out=False); S, N, and k parameters passed into __init__
        Probability mass function
    cdf(n);  S, N, and k parameters passed into __init__
        Cumulative distribution function
    rad(); S, N, and k parameters passed into __init__
        Rank abundance distribution
    fit(sad, guess_for_k=1)
        Maximum likelihood fitting for neg_binom

    Parameters
    ----------
    S : int
        Total number of species in landscape
    N : int
        Total number of individuals in landscape
    k : int
        Aggregation parameter

    Notes
    -----
    S, N, and k are passed into the constructor (__init__) as keyword 
    arguments. Example: trun_neg_binom(S=34, N=345, k=1)



    '''
    
    def pmf(self, n, param_out=False):
        '''
        Probablity mass function of the negative binomial generated with S and
        N

        Parameters
        ----------
        n : int, float or array-like object
            Abundances at which to calculate the pmf

        Passed into __init__
        ---------------------
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
            param_out = True, returns the array as well as the parameter 
            estimates.

        '''
        S = self.params.get('S', None)
        N = self.params.get('N', None)
        k = self.params.get('k', None)
        assert S != None, "S parameter not given"
        assert N != None, "N parameter not given"
        assert k != None, "k parameter not given"
        assert S < N, "S must be less than N"
        assert S > 1, "S must be greater than 1"
        assert N > 0, "N must be greater than 0"
        try:
            len(n); n = np.array(n)
        except:
            n = np.array([n])
        assert np.max(n) <= N, "Maximum n must be less than or equal to N"

        mu = float(N) / S
        p = 1 / (mu / k + 1)  # See Bolker book Chapt 4
        pmf = stats.nbinom.pmf(n, k, p)
        pmf = pmf / (1 - stats.nbinom.pmf(0, k, p)) #Truncation

        if param_out == True:
            return (pmf, {'k' : k, 'p' : p})
        else:
            return pmf   
        
    def fit(self, sad, guess_for_k=1):
        '''
        Fits the truncated negative binomial to the given sad

        Parameters
        ----------
        sad : array-like object
            Observed species abundances 
        guess_for_k : float
            Default parameter for the approximate k given the data

        Returns
        -------
        None

        '''
        def nll_nb(k, sad):
            self.params['k'] = k
            self.params['N'] = np.sum(sad)
            self.params['S'] = len(sad)
            return -sum(np.log(self.pmf(sad)))
        mlek = scipy.optimize.fmin(nll_nb, np.array([guess_for_k]), args=\
                                   (sad,), disp=0)[0]
        self.params['k'] = mlek
        super(trun_neg_binom, self).fit(sad)

class geo_series(Distribution):
    '''
    Geometric series distribution (Motomura 1932 and Magurran 1988).

    Methods
    -------
    rad(); S, N, and k parameters passed into __init__
        Rank abundance distribution
    fit(data)
        Fit k to data

    Parameters
    ----------
    S : int
        Total number of species in landscape
    N : int
        Total number of individuals in landscape
    k : float
        The fraction of resources that each species acquires. Range is 
        (0, 1].

    Notes
    -----
    S, N, and k are passed into the constructor (__init__) as keyword 
    arguments. Example: geo_series(S=34, N=345, k=.4)

    '''

    #TODO:  Need to derive the pmf and cdf for the Geometric series
    def rad(self):
        '''
        Given S, N and k, this function returns the predicted rank-abundance
        distribution for a geometric series distribution (Motomura 1932 and 
        Magurran 1988).

        Parameters
        ----------
        None

        Returns
        -------
        : ndarray (1D)
            An array containing the the predicted SAD with element 0 containing
            thespecies with the most individuals and element S - 1 containing 
            the species with the least individuals.
        '''
        S = self.params.get('S', None)
        N = self.params.get('N', None)
        k = self.params.get('k', None)
        assert S != None, "S parameter not given"
        assert N != None, "N parameter not given"
        assert k != None, "k parameter not given"
        assert S < N, "S must be less than N"
        assert S > 1, "S must be greater than 1"
        assert N > 0, "N must be greater than 0"
        assert k > 0 and k <= 1, "k must be between on the interval (0, 1]"

        C = (1 - (1 - k )** S) ** - 1
        pred_rank_abund = N * C * k * (1 - k) ** (np.arange(1, S + 1) - 1)
        return pred_rank_abund
    
    def fit(self, data):
        '''
        Fit k parameter of geometric series (May 1975)

        Parameters
        ----------
        data : array-like object
            The observed SAD to be fit to a geometric series

        Returns
        -------
        None

        '''
        data = np.array(data)
        S = len(data)
        N = np.sum(data)
        Nmin = np.min(data)
        #Equation from May (1975)
        eq = lambda x: ((x / (1 - x)) * ((1 - x) ** S / (1 - (1 - x) ** S))) -\
                       (Nmin / N)
        k = scipy.optimize.brentq(eq, 1e-10, 1 - 1e-10, disp=True)
        self.params['k'] = k
        super(geo_series, self).fit(data)

class broken_stick(Distribution):
    '''
    McArthur's broken stick species abundance distribution (May 1975)

    Methods
    -------
    pmf(n, param_out=False); S and N parameters passed into __init__
        Probability mass function
    cdf(n); S and N parameters passed into __init__
        Cumulative distribution function
    rad(); S and N parameters passed into __init__
        Rank abundance distribution

    Parameters
    ----------
    S : int
        Total number of species in landscape
    N : int
        Total number of indviduals in landscape

    Notes
    -----
    S and N are passed into the constructor (__init__) as keyword arguments.
    Example: broken_stick(S=22, N=567)

    '''
    #TODO:  PMF is not quite summing to one 
    def pmf(self, n, param_out=False):
        '''
        Probability mass function for MacArthur's broken stick

        Parameters
        ----------
        n : int, float or array-like object
            Abundances at which to calculate the pmf

        Passed into __init__
        ---------------------
        S : int
            Total number of species in landscape
        N : int
            Total number of indviduals in landscape

        Returns
        -------
        : ndarray (1D) or tuple
            Returns array with pmf for values for the given values n. If 
            param_out = True, returns a tuple with the pmf as an array as well 
            as a dictionary of the parameter estimates.
        '''

        S = self.params.get('S', None)
        N = self.params.get('N', None) 
        assert S != None, "S parameter not given"
        assert N != None, "N parameter not given"
        assert S < N, "S must be less than N"
        assert S > 1, "S must be greater than 1"
        assert N > 0, "N must be greater than 0"
        try:
            len(n); n = np.array(n)
        except:
            n = np.array([n])
        #assert np.max(n) <= N, "Maximum n must be less than or equal to N"

        eq = lambda x: ((S - 1) / N) * (1 - (x / N)) ** (S - 2)
        pmf = eq(n)

        if param_out == True:
            return (pmf, {'S' : S})
        else:
            return pmf

    def rad(self):
        '''

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
            An array containing the the predicted SAD with element 0 containing
            the species with the most individuals and element S - 1 containing 
            the species with the least individuals.

        '''
        S = self.params.get('S', None)
        N = self.params.get('N', None) 
        assert S != None, "S parameter not given"
        assert N != None, "N parameter not given"
        assert S < N, "S must be less than N"
        assert S > 1, "S must be greater than 1"
        assert N > 0, "N must be greater than 0"

        pred_rank_abund = np.empty(S)
        for i in xrange(S):
            n = np.arange(i + 1, S + 1) 
            pred_rank_abund[i] = (N / S) * sum(1 / n)
        return pred_rank_abund
            
class plognorm(Distribution):
    '''
    Poisson log-normal distribution (Bulmer 1974)

    Methods
    -------
    pmf(ab, param_out=False); mu and sigma parameters passed into __init__
        Probability mass function
    cdf(ab); mu and sigma parameters passed into __init__
        Cumulative distribution function
    rad(); mu, sigma, S, and N parameters passed into __init__
        Rank abundance distribution
    fit(sad)
        Maximum likelihood estimator for poisson lognormal

    Parameters
    ----------
    mu : float
        the mu parameter of the poisson log normal
    sigma : float
        the sigma parameter of the poisson log normal
    S : int
        Total number of species in landscape
    N : int
        Total number of individuals in landscape

    Notes
    -----
    mu, sigma, S and N are passed into the constructor (__init__) as keyword 
    arguments. Example: plognorm(S=22, N=567, mu=.23, sigma=.89).  S and N
    are needed for plognorm().rad
    '''
    
    def pmf(self, ab, param_out=False):
        '''
        Probability mass function
        
        Parameters
        ----------
        ab : int, float or array-like object
            Abundances at which to calculate the pmf

        Passed into __init__
        ---------------------
        mu : float
            the mu parameter of the poisson log normal
        sigma : float
            the sigma parameter of the poisson log normal

         Returns
        --------
        : ndarray (1D)
            Returns array with pmf for values for the given values ab. If 
            param_out = True, returns the array as well as the parameter 
            estimates.

        Notes
        -----
        This fuction was adopted directly from the VGAM package in R by Mark
        Wilber. The VGAM R package was adopted directly from Bulmer (1974).

        '''
        mu = self.params.get('mu', None)
        sigma = self.params.get('sigma',  None)
        assert mu != None; "mu paramater not given"
        assert sigma != None; "sigma parameter not given"
        try:
            len(ab); ab = np.array(ab)
        except:
            ab = np.array([ab])
        n_unq = np.unique(ab)
    
        if sigma**2 <= 0 or mu <= 0: #Parameters could be negative
            pmf = np.repeat(1e-120, len(n_unq))
        else:
            #TODO: Throwing overflow warning but not affecting result
            eq = lambda t, x: np.exp(t * x - np.exp(t) - 0.5*((t - mu) / \
                                                                    sigma)**2)
            pmf = np.empty(len(n_unq), dtype=np.float)
            for i, n in enumerate(n_unq):
                if n <= 170:
                    integral = integrate.quad(eq, -np.inf, np.inf, args=(n))[0]
                    norm = np.exp((-0.5 * m.log(2 * m.pi * sigma**2) -\
                                                m.lgamma(n + 1)))
                    pmf[i] = norm * integral
                else:
                    z = (m.log(n) - mu) / sigma
                    pmf[i] = (1 + (z**2 + m.log(n) - mu - 1) / (2 * n * 
                             sigma**2)) * np.exp(-0.5 * z**2) / (m.sqrt(2 *\
                             m.pi) * sigma * n)   

        #Only calculated unique abundances to save computational time.
        #Get full pmf again
        pmf_full = np.empty(len(ab))
        for i, n in enumerate(n_unq):
            pmf_full[np.where(ab == n)[0]] = pmf[i]
        pmf = pmf_full

        if param_out == True:
            return (pmf, {'mu' : mu, 'sigma' : sigma})
        else:
            return pmf      
    
    def fit(self, data):
        '''
        Maximum likelihood Estimates for Poisson log normal

        Parameter
        ---------
        data : array-like object
            The observed abundances which will be fit to a poisson lognormal

        Returns
        -------
        None

        Notes
        -----
        This function was adapted from Ethan White's pln_solver function in 
        weecology.

        '''
        try:
            len(data); data = np.array(data)
        except:
            data = np.array([data])
        assert len(data) >= 1, "Data length must be greater than or equal to 1"
        mu0 = np.mean(np.log(data))
        sigma0 = np.std(np.log(data), ddof=1)
        def pln_func(x):
            self.params['mu'] = x[0]
            self.params['sigma'] = x[1]
            return -sum(np.log(self.pmf(data)))
        mu, sigma = scipy.optimize.fmin(pln_func, x0 = [mu0, sigma0], disp=0)
        self.params['mu'] = mu
        self.params['sigma'] = sigma
        super(plognorm, self).fit(data)
        

class trun_plognorm(Distribution):
    '''
    Truncated Poisson lognormal (Bulmer 1974)

    Methods
    -------
    pmf(ab, param_out=False); mu and sigma parameters passed into __init__
        Probability mass function
    cdf(ab); mu and sigma parameters passed into __init__
        Cumulative distribution function
    rad(); mu, sigma, S, and N parameters passed into __init__
        Rank abundance distribution
    fit(sad)
        Maximum likelihood estimator for poisson lognormal

    Parameters
    ----------
    mu : float
        the mu parameter of the poisson log normal
    sigma : float
        the sigma parameter of the poisson log normal
    S : int
        Total number of species in landscape
    N : int
        Total number of individuals in landscape
    
    Notes
    -----
    mu, sigma, S and N are passed into the constructor (__init__) as keyword 
    arguments. Example: trun_plognorm(S=22, N=567, mu=.23, sigma=.89).  S and N
    are needed for trun_plognorm().rad


    '''

    def pmf(self, ab, param_out=False):
        '''
        Truncated Poisson log-normal (Bulmer 1974)

        Parameters
        ----------
        ab : int, float or array-like object
            Abundances at which to calculate the pmf

        Passed into __init__
        ---------------------
        mu : float
            the mu parameter of the poisson log normal
        sigma : float
            the sigma parameter of the poisson log normal
        
        Returns
        -------
        : ndarray (1D)
            Returns array with pmf for values for the given values ab. If 
            param_out = True, returns the array as well as a dict with
            parameter estimates.

        '''
        mu = self.params.get('mu', None)
        sigma = self.params.get('sigma',  None)
        assert mu != None; "mu paramater not given"
        assert sigma != None; "sigma parameter not given"
        if param_out == True:
            plgn = plognorm(mu=mu, sigma=sigma)
            untr_pmf, params = plgn.pmf(ab, param_out=param_out)
            pmf0 = plgn.pmf(0)
            tr_pmf = (untr_pmf / (1 - pmf0))#Truncating based on Bulmer eq. A1
            return (tr_pmf, {'mu' : params['mu'], 'sigma' : params['sigma']})
        else:
            plgn = plognorm(mu=mu, sigma=sigma)
            untr_pmf = plgn.pmf(ab)
            pmf0 = plgn.pmf(0)
            tr_pmf = (untr_pmf / (1 - pmf0))
            return tr_pmf
            
    def fit(self, sad):
        '''
        Maximum likelihood estimates for truncated poisson log normal

        Parameter
        ---------
        sad : array-like object
            The observed abundances

        Returns
        -------
        : dict
            The maximum likelihood estimates of mu and sigma

        Notes
        -----
        This function was adapted from Ethan White's pln_solver function in 
        weecology.
        '''
        assert type(sad) == list or type(sad) == tuple \
           or type(sad) == np.ndarray, "Invalid parameter type"
        assert len(sad) >= 1, "len(sad) must be greater than or equal to 1"
        sad = np.array(sad)
        mu0 = np.mean(np.log(sad))
        sigma0 = np.std(np.log(sad), ddof=1)
        def pln_func(x):
            self.params['mu'] = x[0]
            self.params['sigma'] = x[1]
            return -sum(np.log(self.pmf(sad)))
        mu, sigma = scipy.optimize.fmin(pln_func, x0 = [mu0, sigma0], disp=0)
        self.params['mu'] = mu
        self.params['sigma'] = sigma
        super(trun_plognorm, self).fit(sad)

class lognorm(Distribution):
    '''
    Lognormal distribution

    Methods
    -------
    pmf(ab, param_out=False); mu and sigma parameters passed into __init__
        Probability mass function
    cdf(ab); mu and sigma parameters passed into __init__
        Cumulative distribution function
    rad(); mu, sigma, S, and N parameters passed into __init__
        Rank abundance distribution
    fit(sad)
        Maximum likelihood estimator for lognormal

    Parameters
    ----------
    mu : float
        the mu parameter of the log normal
    sigma : float
        the sigma parameter of the poisson log normal
    S : int
        Total number of species in landscape
    N : int
        Total number of individuals in landscape
    
    Notes
    -----
    mu, sigma, S and N are passed into the constructor (__init__) as keyword 
    arguments. Example: lognorm(S=22, N=567, mu=.23, sigma=.89).  S and N
    are needed for lognorm().rad

    Log-normal distributions are continuous and there for the cdf should be an
    integral. However, the cdf of an integral from a to a is 0 and the
    probability of there being a single individual within a species given an
    SAD is certainly not 0.  Therefore, we consider the lognormal "discrete"
    and calcuate the cdf by summing.  Note that this is one of the many 
    problems that Blackburn and Gaston have with using the lognormal for SADs.


    '''
    
    def pmf(self, ab, param_out=False):
        '''
        Lognormal pmf

        Parameters
        ----------
        ab : int, float or array-like object
            Abundances at which to calculate the pmf

        Pass into __init__
        ------------------
        mu : float
            the mu parameter of the poisson log normal
        sigma : float
            the sigma parameter of the poisson log normal

        Returns
        -------
        : ndarray (1D)
            Returns array with pmf for values for the given values ab. If 
            param_out = True, returns the array as well as the parameter 
            estimates.

        Notes
        -----
        scipy.stats.lognorm is coded very unclearly and the docstring is not 
        helpful so we coded our own lognormal distribution
        '''

        mu = self.params.get('mu', None)
        sigma = self.params.get('sigma',  None)
        assert mu != None; "mu paramater not given"
        assert sigma != None; "sigma parameter not given"
        try:
            len(ab); ab = np.array(ab)
        except:
            ab = np.array([ab])
        pmf = stats.norm.pdf(np.log(ab), loc=mu, scale=sigma) / ab
        if param_out == True:
            return (pmf, {'mu' : mu, 'sigma' : sigma})
        else:
            return pmf
        
    def fit(self, sad):
        '''
        Maximum likelihood Estimates for log normal

        Parameter
        ---------
        sad : array-like object
            The observed abundances

        Returns
        -------
        None

        '''
        assert len(sad) >= 1, "len(sad) must be greater than or equal to 1"
        sad = np.array(sad)
        mu0 = np.mean(np.log(sad))
        sigma0 = np.std(np.log(sad), ddof=1)
        def ln_func(x):
            self.params['mu'] = x[0]
            self.params['sigma'] = x[1]
            return -sum(np.log(self.pmf(sad)))
        mu, sigma = scipy.optimize.fmin(ln_func, x0 = [mu0, sigma0], disp=0)
        self.params['mu'] = mu
        self.params['sigma'] = sigma
        super(lognorm, self).fit(sad)

class sugihara(Distribution):
    '''
    Sugihara Rank Abundance Distribution (Sugihara 1980)

    Methods
    -------
    rad(sample_size=10000); N and S parameters are passed in __init__

    Notes
    -----
    S and N are passed into the constructor (__init__) as keyword 
    arguments. Example: sugihara(S=22, N=567)

    '''
    #TODO: Back-derive the pmf?

    def rad(self, sample_size=10000):
        '''
        Parameters
        ----------
        None

        Passed into __init__
        ---------------------
        S : int
            Total number of species in landscape
        N : int
            Total number of indviduals in landscape

        Returns
        -------
        : 1D np.ndarray
            The predicted rank abundance distribution of the species in the
            landscape

        Notes
        -----
        As S gets large, we will start to under sample some of the possible 
        breakage sequences and this model begins to fail. Adapted from 
        Sugihara (1980)

        '''
        S = self.params.get('S', None)
        N = self.params.get('N', None) 
        assert S != None, "S parameter not given"
        assert N != None, "N parameter not given"
        total = []
        for i in xrange(sample_size):
            U = np.random.triangular(0.5, 0.75, 1, size=S - 1)
            p = []
            #Could this be refactored to look sexier and perform better?
            for i in xrange(S):
                if i == 0:
                    p.append(1)
                else:
                    index = np.random.random_integers(0, len(p) - 1)
                    p_new1 = p[index] * U[i - 1]
                    p_new2 = p[index] * (1 - U[i - 1])
                    p[index] = p_new1
                    p.append(p_new2)
                    p.sort(reverse=True)
            total.append(p)
      
        total_array = np.array(total)
        means = []
        for i in xrange(S):
            means.append(np.mean(total_array[:,i]))
        means = np.array(means)
        return N * means

class mete_lgsr(Distribution):
    '''
    METE truncated log series (Harte 2011)

    Methods
    -------
    pmf(n, param_out); S and N parameters passed into __init__
        Probability mass function
    cdf(n); S and N parameters passed into __init__
        Cumulative distribution function
    rad(); S and N parameters passed into __init__
        Rank abundance distribution

    Parameters
    ----------
    S : int
        Total number of species in landscape
    N : int
        Total number of indviduals in landscape

    Notes
    -----
    S and N are passed into the constructor (__init__) as keyword arguments. 
    Example: mete_lgsr(S=22, N=567)
    
    '''
    
    def pmf(self, n, param_out=False):
        '''
        Truncated log series pmf (Harte 2011)

        Parameters
        ----------
        n : int, float or array-like object
            Abundances at which to calculate the pmf
        
        Passed into __init__
        ---------------------
        S : int
            Total number of species in landscape
        N : int
            Total number of indviduals in landscape

        Returns
        -------
        : ndarray (1D)
            Returns array with pmf for values for the given values n. If 
            param_out = True, returns the array as well as the parameter 
            estimates.

        Notes
        -----
        This function uses the truncated log series as described in Harte 2011
        eq (7.32).  The equation used in this function to solve for the 
        Lagrange multiplier is equation (7.27) as described in Harte 2011. 
    
        Also note that realistic values of x where x = e**-(beta) 
        (see Harte 2011) are in the range (1/e, 1) (exclusive). Therefore, the
        start and stop parameters for the brentq procedure are close to these 
        values. However, x can occasionally be greater than one so the maximum
        stop value of the brentq optimizer is 2.

        '''
        S = self.params.get('S', None)
        N = self.params.get('N', None) 
        assert S != None, "S parameter not given"
        assert N != None, "N parameter not given"
        assert S < N, "S must be less than N"
        assert S > 1, "S must be greater than 1"
        assert N > 0, "N must be greater than 0"
        try:
            len(n); n = np.array(n)
        except:
            n = np.array([n])
        assert np.max(n) <= N, "Maximum n must be less than or equal to N"
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
    
class mete_lgsr_approx(Distribution):
    '''
    METE log series distribution with approximation (Harte 2011)

    Methods
    -------
    pmf(n, param_out); S and N parameters passed into __init__
        Probability mass function
    cdf(n); S and N parameters passed into __init__
        Cumulative distribution function
    rad(); S and N parameters passed into __init__
        Rank abundance distribution

    Parameters
    ----------
    S : int
        Total number of species in landscape
    N : int
        Total number of indviduals in landscape

    Notes
    -----
    S and N are passed into the constructor (__init__) as keyword arguments. 
    Example: mete_lgsr_approx(S=22, N=567)

    '''
    
    def pmf(self, n, param_out=False, root=2):
        '''
        Truncated log series using approximation (7.30) and (7.32) in Harte 
        2011

        Parameters
        ----------
        n : int, float or array-like object
            Abundances at which to calculate the pmf

        root: int (optional)
            1 or 2.  Specifies which root to use for pmf calculations   

        Passed into __init__
        ---------------------
        S : int
            Total number of species in landscape
        N : int
            Total number of indviduals in landscape

        Returns
        -------
        : ndarray (1D)
            Returns array with pmf for values for the given values n. If 
            param_out = True, returns the array as well as the parameter 
            estimates.
        
        Notes:
        ------
        This function uses the truncated log series as described in Harte 2011
        eq (7.32).  The equation used in this function to solve for the 
        Lagrange multiplier is equation (7.30) as described in Harte 2011.     
       
        Also note that realistic values of x where x = e^-(beta) (see Harte 
        2011) are in the range (1/e, 1) (exclusive). Therefore, the start and 
        stop parameters for the brentq optimizer have been chosen near these 
        values.
        
        '''
        S = self.params.get('S', None)
        N = self.params.get('N', None) 
        assert S != None, "S parameter not given"
        assert N != None, "N parameter not given"
        assert S < N, "S must be less than N"
        assert S > 1, "S must be greater than 1"
        assert N > 0, "N must be greater than 0"
        try:
            len(n); n = np.array(n)
        except:
            n = np.array([n])
        assert np.max(n) <= N, "Maximum n must be less than or equal to N"

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
                print "Warning: More than one root."
                if root == 1:
                    x = scipy.optimize.brentq(eq, start, xmax, disp=True)
                if root == 2:
                    x = scipy.optimize.brentq(eq, xmax, stop, disp=True)
                       
            if ymax < 0:
                raise RootError('No solution to constraint equation with ' +\
                                ' given values of S and N') 
        g = -1/m.log(x)
        pmf = (1/m.log(g)) * ((x**n)/n) 
    
        if param_out == True:
            return (pmf, {'beta' : -np.log(x)})
        else:
            return pmf

class binm(Distribution):
    ''' Binomial distribution (ie, random placement model)

    Methods
    -------
    pmf(n); N and a parameters are passed in the __init__ method
        Probability mass function
    cdf(n); N an a parameters are passed in the __init__method
        Cumulative distribution function
    
    Parameters
    ----------
    N : ndarray or int
        Total number of individuals in landscape
    a : ndarray or int
        Ratio of cell size to area of whole landscape

    Notes
    -----
    N and a are passed into the constructor (__init__) as keyword arguments. 
    Example: binm(N=567, a=0.2)

    '''

    def pmf(self, n):
        '''
        Binomial pmf (ie, random placement model).

        Parameters
        ----------
        n : ndarray or int
            Values of n for which to calculate pmf or likelihood

        Pass to __init__
        -----------------
        N : ndarray or int
            Total number of individuals in landscape
        a : ndarray or int
            Ratio of cell size to area of whole landscape

        Returns
        -------
        : ndarray
            Returns array with pmf.
    
        '''
        N = self.params.get('N', None)
        a = self.params.get('a', None)
        assert N != None, "N parameter not given"
        assert a != None, "a parameter not given"
        assert a < 1, "a must be less than 1"
        try:
            len(n); n = np.array(n)
        except:
            n = np.array([n])
        assert np.max(n) <= N, "n maximum must be less than or equal to N"
        return stats.binom.pmf(n, N, a)

    def cdf(self, n):
        '''
        Cumulative distribution function

        Parameters
        ----------
        n : ndarray or int
            Values of n for which to calculate the cdf

        Pass to __init__
        -----------------
        N : ndarray or int
            Total number of individuals in landscape
        a : ndarray or int
            Ratio of cell size to area of whole landscape

        Returns
        -------
        : ndarray
            Returns array with cdf.
        
        '''
        N = self.params.get('N', None)
        a = self.params.get('a', None)
        assert N != None, "N parameter not given"
        assert a != None, "a parameter not given"
        assert a < 1, "a must be less than 1"
        try:
            len(n); n = np.array(n)
        except:
            n = np.array([n])
        assert np.max(n) <= N, "n maximum must be less than or equal to N"
        return scipy.stats.binom.cdf(n, N, a)

class pois(Distribution):
    '''
    Poisson distribution

    Methods
    -------
    pmf(n); N and a parameters are passed in the __init__ method
        Probability mass function
    cdf(n); N an a parameters are passed in the __init__method
        Cumulative distribution function
    
    Parameters
    ----------
    N : ndarray or int
        Total number of individuals in landscape
    a : ndarray or int
        Ratio of cell size to area of whole landscape

    Notes
    -----
    N and a are passed into the constructor (__init__) as keyword arguments. 
    Example: pois(N=567, a=0.2)


    '''
    
    def pmf(self, n):
        '''
        Poisson pmf.

        Parameters
        ----------
        n : ndarray or int
            Values of n for which to calculate pmf or likelihood

        Pass to __init__
        -----------------
        N : ndarray or int
            Total number of individuals in landscape
        a : ndarray or int
            Ratio of cell size to area of whole landscape

        Returns
        -------
        : ndarray
            Returns array with pmf. 

        '''
        N = self.params.get('N', None)
        a = self.params.get('a', None)
        assert N != None, "N parameter not given"
        assert a != None, "a parameter not given"
        assert a < 1, "a must be less than 1"
        try:
            len(n); n = np.array(n)
        except:
            n = np.array([n])
        assert np.max(n) <= N, "n maximum must be less than or equal to N" 
        mu = N * a
        return scipy.stats.poisson.pmf(n, mu)
    
    def cdf(self, n):
        '''
        Cumulative distribution function

        Parameters
        ----------
        n : ndarray or int
            Values of n for which to calculate the cdf

        Pass to __init__
        -----------------
        N : ndarray or int
            Total number of individuals in landscape
        a : ndarray or int
            Ratio of cell size to area of whole landscape

        Returns
        -------
        : ndarray
            Returns array with cdf.

        '''

        N = self.params.get('N', None)
        a = self.params.get('a', None)
        assert N != None, "N parameter not given"
        assert a != None, "a parameter not given"
        assert a < 1, "a must be less than 1"
        try:
            len(n); n = np.array(n)
        except:
            n = np.array([n])
        assert np.max(n) <= N, "n maximum must be less than or equal to N"
        mu = N * a
        return scipy.stats.poisson.cdf(n, mu)

class nbd(Distribution):
    '''
    Negative binomial distribution

    Methods
    -------
    pmf(n); N, a, k parameters are passed in the __init__ method
        Probability mass function
    cdf(n); N, a, and k parameters are passed in the __init__ method
        Cumulative distribution function
    fit(data): a parameter is passed in the __init__ method

    Parameters
    ----------
    N : ndarray or int
        Total number of individuals in landscape
    a : ndarray or int
        Ratio of cell size to area of whole landscape
    k : int
        Aggregation parameter

    Notes
    -----
    N, a, and k are passed into the constructor (__init__) as keyword arguments. 
    Example: nbd(N=567, a=0.2, k=1)


    '''

    def pmf(self, n):

        '''
        Negative binomial pmf.

        Parameters
        ----------
        n : ndarray or int
            Values of n for which to calculate pmf

        Pass in __init__
        -----------------
        N : ndarray or int
            Total number of individuals in landscape
        a : ndarray or int
            Ratio of cell size to area of whole landscape
        k : int
            Aggregation parameter

        Returns
        -------
        : ndarray
            Returns array with pmf.
        '''
        N = self.params.get('N', None)
        a = self.params.get('a', None)
        k = self.params.get('k', None)
        assert N != None, "N parameter not given"
        assert a != None, "a parameter not given"
        assert k != None, "k parameter not given"
        assert a < 1, "a must be less than 1"
        try:
            len(n); n = np.array(n)
        except:
            n = np.array([n])
        assert np.max(n) <= N, "n maximum must be less than or equal to N"
        mu = N * a
        p = 1 / (mu / k + 1)  # See Bolker book Chapt 4
        return scipy.stats.nbinom.pmf(n, k, p)
    
    def cdf(self, n):
        '''
        Cumuluative distribution function

        Parameters
        ----------
        n : ndarray or int
            Values of n for which to calculate the cdf

        Pass in __init__
        ----------------
        N : ndarray or int
            Total number of individuals in landscape
        a : ndarray or int
            Ratio of cell size to area of whole landscape
        k : int
            Aggregation parameter

        Returns
        -------
        : ndarray or float
            Returns array with cdf.

        '''
        N = self.params.get('N', None)
        a = self.params.get('a', None)
        k = self.params.get('k', None)
        assert N != None, "N parameter not given"
        assert a != None, "a parameter not given"
        assert k != None, "k parameter not given"
        assert a < 1, "a must be less than 1"
        try:
            len(n); n = np.array(n)
        except:
            n = np.array([n])
        assert np.max(n) <= N, "n maximum must be less than or equal to N"
        mu = N * a
        p = 1 / (mu / k + 1)  # See Bolker book Chapt 4
        return scipy.stats.nbinom.cdf(n, k, p)

    def fit(self, data, guess_for_k=1):
        '''
        Fits negative binomial to data

        Parameters
        ----------
        data : array-like object
            Individuals per cell
        guess_for_k : float
            Default parameter for the approximate k given the data

        Pass to __init__
        -----------------
        a : ndarray or int
            Ratio of cell size to area of whole landscape


        Returns
        -------
        : dict
            The maximum likelihood estimator (MLE) for k
        
        '''
        a = self.params.get('a', None)
        assert a != None, "a parameter not given"
        def nll_nb(k, data):
            self.params['k'] = k
            self.params['N'] = sum(data)
            return -sum(np.log(self.pmf(data)))
        mlek = scipy.optimize.fmin(nll_nb, np.array([guess_for_k]), args=\
                                   (data,), disp=0)[0]
        self.params['k'] = mlek
        super(nbd, self).fit(data, sad=False)

class fnbd(Distribution):
    '''
    Finite negative binomial (Zillio and He 2010).

    Methods
    -------
    pmf(n); N, a, k parameters are passed in the __init__ method
        Probability mass function
    cdf(n); N, a, and k parameters are passed in the __init__ method
        Cumulative distribution function
    fit(data): a parameter is passed in the __init__ method

    Parameters
    ----------
    N : ndarray or int
        Total number of individuals in landscape
    a : ndarray or int
        Ratio of cell size to area of whole landscape
    k : int
        Aggregation parameter

    Notes
    -----
    N, a, and k are passed into the constructor (__init__) as keyword arguments. 
    Example: fnbd(N=567, a=0.2, k=1)

    
    '''
    
    def pmf(self, n):
        '''
        Finite negative binomial pmf (Zillio and He 2010).

        Parameters
        ----------
        n : ndarray or int
            Values of n for which to calculate pmf

        Pass in __init__
        -----------------
        N : ndarray or int
            Total number of individuals in landscape
        a : ndarray or int
            Ratio of cell size to area of whole landscape
        k : int
            Aggregation parameter

        Returns
        -------
        : ndarray
            Returns array with pmf.

        Notes
        -----
        If any arg other than n is iterable, all args must be iterable of same 
        length.

        The fnbd with k = 1 is not a true geometric distribution - calculate a 
        pmf z and run z[1:]/z[0:-1], noting that the ratio is not constant.

        '''

        # TODO: Fix to work if n and N are one value
        #    if not (n <= N).all():
        #        raise Exception, "All values of n must be <= N."
        #    elif (a <= 0) or (a >= 1):
        #        raise Exception, "a must be between 0 and 1"
        N = self.params.get('N', None)
        a = self.params.get('a', None)
        k = self.params.get('k', None)
        assert N != None, "N parameter not given"
        assert a != None, "a parameter not given"
        assert k != None, "k parameter not given"
        try:
            len(n); n = np.array(n)
        except:
            n = np.array([n])
        assert np.max(n) <= N, "n maximum must be less than or equal to N"
        ln_L = lambda n_i,N,a,k: _ln_choose(n_i+k-1,n_i) + \
            _ln_choose(N-n_i+(k/a)-k-1,N-n_i) - _ln_choose(N +(k/a)-1,N)

        pmf = ln_L(n, N, a, k)  # Already log

        return np.exp(pmf)
    
    def cdf(self, n):
        '''
        See Distribution.cdf() docstring
        '''
        return super(fnbd, self).cdf(n, lower_bound=0)
            
    def fit(self, data, upper_bnd=5):
        '''
        Fits finite negative binomial to data

        Parameters
        ----------
        data : array-like object
            Individuals per cell
        upper_bnd : float
            Default parameter for upper_bnd on k estimate

        Pass to __init__
        -----------------
        a : ndarray or int
            Ratio of cell size to area of whole landscape


        Returns
        -------
        : dict
            The maximum likelihood estimator (MLE) for k 
        
        '''

        a = self.params.get('a', None)
        assert a != None, "a parameter not given"
        def nll_nb(k, data):
            self.params['k'] = k
            self.params['N'] = sum(data)
            return -sum(np.log(self.pmf(data)))
        mlek = scipy.optimize.brute(nll_nb, ((1e-10,upper_bnd),), args=(data,))
        self.params['k'] = mlek
        super(fnbd, self).fit(data, sad=False)

class geo(Distribution):
    '''
    Geometric distribution.  Uses nbd object as wrapper with k = 1

    Methods
    -------
    pmf(n); N and a parameters are passed in the __init__ method
        Probability mass function
    cdf(n); N and a parameters are passed in the __init__ method
        Cumulative distribution function

    Parameters
    ----------
    N : ndarray or int
        Total number of individuals in landscape
    a : ndarray or int
        Ratio of cell size to area of whole landscape

    Notes
    -----
    N and aare passed into the constructor (__init__) as keyword arguments. 
    Example: geo(N=567, a=0.2)

        
    '''

    def pmf(self, n):
        '''
        Geometric pmf

        Parameters
        ----------
        n : ndarray or int
            Values of n for which to calculate pmf

        Pass in __init__
        -----------------
        N : ndarray or int
            Total number of individuals in landscape
        a : ndarray or int
            Ratio of cell size to area of whole landscape

        Returns
        -------
        : ndarray
            Returns array with pmf.
            
        '''
        
        N = self.params.get('N', None)
        a = self.params.get('a', None)
        assert N != None, "N parameter not given"
        assert a != None, "a parameter not given"
        return nbd(N=N, a=a, k=1).pmf(n)
    
    def cdf(self, n):
        '''
        Cumulative distribution function

        Parameters
        ----------
        n : ndarray or int
            Values of n for which to calculate pmf 

        Pass from __init__
        ------------------
        N : ndarray or int
            Total number of individuals in landscape
        a : ndarray or int
            Ratio of cell size to area of whole landscape

        Returns
        -------
        : ndarray
            Returns array with pmf.

        '''

        N = self.params.get('N', None)
        a = self.params.get('a', None)
        assert N != None, "N parameter not given"
        assert a != None, "a parameter not given"
        return nbd(N=N, a=a, k=1).cdf(n)

class fgeo(Distribution):
    '''
    Finite geometric pmf (Zillio and He 2010). Use fnbd object as wrapper with 
    k = 1

    Methods
    -------
    pmf(n); N and a parameters are passed in the __init__ method
        Probability mass function
    cdf(n); N and a parameters are passed in the __init__ method
        Cumulative distribution function

    Parameters
    ----------
    N : ndarray or int
        Total number of individuals in landscape
    a : ndarray or int
        Ratio of cell size to area of whole landscape
    
    Notes
    -----
    N and aare passed into the constructor (__init__) as keyword arguments. 
    Example: fgeo(N=567, a=0.2)

    '''
    
    def pmf(self, n):
        '''
        Finite geometric pmf

        Parameters
        ----------
        n : ndarray or int
            Values of n for which to calculate pmf

        Pass to __init__
        -----------------
        N : ndarray or int
            Total number of individuals in landscape
        a : ndarray or int
            Ratio of cell size to area of whole landscape

        Notes
        -----
        This is not a true geometric distribution - calculate a pmf z and run 
        z[1:]/z[0:-1], noting that the ratio is not constant.

        '''
        
        N = self.params.get('N', None)
        a = self.params.get('a', None)
        assert N != None, "N parameter not given"
        assert a != None, "a parameter not given"
        return fnbd(N=N, a=a, k=1).pmf(n)
    
    def cdf(self, n):
        '''
        Cumulative distribution function

        Parameters
        ----------
        n : ndarray or int
            Values of n for which to calculate pmf

        Pass to __init__
        -----------------
        N : ndarray or int
            Total number of individuals in landscape
        a : ndarray or int
            Ratio of cell size to area of whole landscape

        '''
        
        N = self.params.get('N', None)
        a = self.params.get('a', None)
        assert N != None, "N parameter not given"
        assert a != None, "a parameter not given"
        return fnbd(N=N, a=a, k=1).cdf(n)

class tgeo(Distribution):
    '''
    Truncated geometric distribution (Harte 2011)

    Methods
    -------
    pmf(n); N and a parameters are passed in the __init__ method
        Probability mass function
    cdf(n); N and a parameters are passed in the __init__ method
        Cumulative distribution function

    Parameters
    ----------
    N : ndarray or int
        Total number of individuals in landscape
    a : ndarray or int
        Ratio of cell size to area of whole landscape

    Notes
    -----
    N and aare passed into the constructor (__init__) as keyword arguments. 
    Example: tgeo(N=567, a=0.2)


    '''

    
    def pmf(self, n, param_out=False):
        '''
        Truncated geometric pmf (Harte et al 2008).

        Parameters
        ----------
        n : ndarray or int
            Values of n for which to calculate pmf or likelihood
        stop : int or float
            The maximum value for the brentq solver

        Pass to __init__
        -----------------
        N : ndarray or int
            Total number of individuals in landscape
        a : ndarray or int
            Ratio of cell size to area of whole landscape

        Returns
        -------
        : ndarray
            Returns array with pmf.  If param_out is True, returns a tuple with
            the pmf as well as the lagrange multiplier lambda_pi within a
            dictionary.

        Notes
        -----
        N and a must be the same length.

        This is a true geometric distribution in which p = exp(-lambda). 
        Confirmed with z[1:]/z[0:-1].
`
        Plotting eq for various realistic values of p, N, and a confirms that 
        this is a smooth, well-behaved function that should be amenable to 
        using a root finder. When a is large and N is small (i.e. a = .9 and
        N = 32) the lagrange multiplier has no solution.  For large a's (.8 or
        greater), N needs to be sufficiently large for lambda pi to have a 
        solution.

        '''
        N = self.params.get('N', None)
        a = self.params.get('a', None)
        assert N != None, "N parameter not given"
        assert a != None, "a parameter not given"
        try:
            len(n); n = np.array(n)
        except:
            n = np.array([n])
        assert np.max(n) <= N, "n maximum must be less than or equal to N"
        eq = lambda x: ((x / (1 - x)) - (((N + 1) * x ** (N + 1)) / \
                            (1 - x ** (N + 1)))) - (N * a)
        x = scipy.optimize.brentq(eq, 0, min((sys.float_info[0] *\
                                            a)**(1/float(N)), 2), disp=False)
        z = (1 - x ** (N + 1)) / (1 - x)
        pmf = (1 / z) * (x ** n)
        if param_out:
            return pmf, {'lambda_pi' :  -np.log(x)}
        else:
             return pmf 

    def cdf(self, n):
        '''
        Cumulative distribution function

        Parameters
        ----------
        n : ndarray or int
            Values of n for which to calculate pmf 

        Pass from __init__
        ------------------
        N : ndarray or int
            Total number of individuals in landscape
        a : ndarray or int
            Ratio of cell size to area of whole landscape

        Returns
        -------
        : ndarray
            Returns array with pmf.
        '''
        return super(tgeo, self).cdf(n, lower_bound=0)

class Mete_sar(object):
    '''
    This class explores the METE generated SAR

    Methods
    -------
    mete_sar_method1(anchor_area, upscale=0, downscale=0, target_area=None)
    univ_curve(num_iter=10)


    Parameters passed to __init__
    ------------------------------
    S : int
        Total number of species at the given anchor area
    N : int
        Total number of individuals at the given anchor area

    Notes
    -----
    S and N should be passed as keyword arguments to the constructor.
    Example: Mete_sar(S=34, N=235)

    '''

    def __init__(self, **kwargs):
        self.params = kwargs
    
    def mete_sar_method1(self, anchor_area, upscale=0, downscale=0,\
                     target_area=None):
        '''
        Predict the universal SAR curve for the given S and N found at 
        the given anchor scale

        Parameters
        ----------
        anchor_area : float
            The area from which the SAR will be upscaled or downscaled.
        upscale : int
            Number of iterations up from the anchor scale.  Each iteration 
            doubles the previous area.
        downscale : int
            Number of iterations down from the anchor scale. Each iteration 
            halves the previous area.
        target_area : float
            The desired area for the species-area relationship.  If not None, 
            this keyword argument overrides the upscale and downscale 
            arguements.

        Passed to __init__
        -------------------
        S : int
            Total number of species at the given anchor area
        N : int
            Total number of individuals at the given anchor area


        Returns
        -------
        : 1D structured np.array
            The structured array has fields 'species' and 'area'

        Notes
        -----
        This function uses equations 3, 7, 8, and 9 found in Harte et al. 
        (2009). When possible, the approximation 
        sum(x**n / n) ~= log(1 / log( 1/x)) was used to decrease runtime.  
    
    
        '''
        S = self.params.get('S', None)
        N = self.params.get('N', None)
        assert S != None, "S parameter not given"
        assert N != None, "N parameter not given"
        assert S < N, "S must be less than N"
        assert S > 1, "S must be greater than 1"
        assert N > 0, "S must be greater than 0"
        if target_area != None:
            if target_area == anchor_area:
                upscale = 0; downscale = 0
            elif target_area > anchor_area:
                upscale = np.int(np.ceil(np.log2(target_area / anchor_area)))
                downscale = 0
            elif target_area < anchor_area:
                downscale = np.int(np.ceil(np.abs(np.log2(target_area /\
                                                      anchor_area)))) 
                upscale = 0
    
        if upscale == 0 and downscale == 0:
            return np.array((S, anchor_area), dtype=[('species', np.float),\
                                                ('area', np.float)])
        areas = _generate_areas_(anchor_area, upscale, downscale)
        sar = np.empty(len(areas), dtype=[('species', np.float),\
                                      ('area', np.float)])
        sar['area'] = areas
        if upscale != 0: 
            sar['species'][downscale:] = _upscale_sar_(areas[downscale:], N, S)
        if downscale != 0:
            sar['species'][:downscale + 1] =\
                                   _downscale_sar_(areas[:downscale + 1], N, S)
        return sar

    def univ_curve(self, num_iter=10):
        '''
        This function calculates the universal SAR curve.

        Parameters
        ----------
        num_iter : int
            Number of iterations.
            WARNING: Running more than 10 ten begins to take along time

        Passed into __init__
        --------------------
        S : int
            Total number of species at the given anchor area
        N : int
            Total number of individuals at the given anchor area
    
        Returns
        -------
        :1D structured np.array
            The structured array has fields 'z' and 'N/S' for plotting 
            convenience

        Notes
        -----
        This function calculates the universal SAR curve.  Different ratios
        of N/S will only put you at different places along the same curve. 
        Iterations of more than 15 take a long time. The equations used in this
        function were taken from Harte et al. (2009) and Harte (2011). This uses 
        method 1 in Harte (2011)

        '''
        
        S = self.params.get('S', None)
        N = self.params.get('N', None)
        assert S != None, "S parameter not given"
        assert N != None, "N parameter not given"
        assert S < N, "S must be less than N"
        assert S > 1, "S must be greater than 1"
        assert N > 0, "S must be greater than 0"
        num_ind = np.empty(num_iter + 1)
        spp = np.empty(num_iter + 1)
        univ_SAR = np.empty(num_iter + 1, dtype=[('z', np.float),\
                                             ('N/S', np.float)])
        #From Harte et al. (2009)
        z = lambda beta: 1 / (np.log(2) * np.log(1 / beta))

        for i in xrange(num_iter + 1):
            if i == 0:
                num_ind[i] = N
                spp[i] = S
                n = np.linspace(1, N, num=N)
                eq = lambda x: (((x - x ** (N + 1)) / ( 1 - x)) / \
                                                   sum((x ** n) / n)) - (N / S)
                x = scipy.optimize.brentq(eq, 1e-10, min((sys.float_info[0] / S)\
                                      **(1/float(N)), 2), disp=True)
                univ_SAR['z'][i] = z(-np.log(x))
                univ_SAR['N/S'][i] = N / S
            else:
                num_ind[i] = 2 * num_ind[i - 1]
                N2A = num_ind[i]
                n = np.linspace(1, N2A, num=N2A)
                eq1 = lambda x: (sum((x**n)/n) * ((N2A) / ((x - x**(N2A + 1)) / \
                            (1 - x)) * (1 / x))) - ((N2A) * ((1 - x) / \
                            (x - x**(N2A + 1))) * (1 - (x**N2A / (N2A + 1))))\
                            - spp[i - 1]
                x = scipy.optimize.brentq(eq1, 1e-10, 1 - 1e-10, disp=True)
                univ_SAR['z'][i] = z(-np.log(x))
                #Using approximations for summations. See module notes
                eq2 = lambda S2A: np.log(1 / np.log(1 / x)) * ((N2A) / \
                                          ((x - x**(N2A + 1)) / (1 - x))) - S2A
                S2A = scipy.optimize.brentq(eq2, spp[i - 1], 2 * spp[i - 1], \
                                                                     disp=True)
                spp[i] = S2A
                univ_SAR['N/S'][i] = N2A / S2A
        return univ_SAR

class SAR(object):
    '''
    This class contains different functions for examining an SAR.  All \
    METE-related SAR functions can be found in Mete_SAR class.

    Methods
    -------
    power_law(area_list, S, anchor_area, z)

    '''
    def __init__(self, **kwargs):
        self.params = kwargs

    def power_law(self, area_list, S, anchor_area, z):
        '''
        Generate a power law SAR with a given z for some S at an anchor area

        Parameters
        ----------
        area_list : array-like object
            An array-like object containing SAR areas to be computed
        S : int
            Total number of species at the given anchor area
        anchor_area : float
            The area which contains S species
        z : int
            The power of the power law
    
        Returns
        -------
        : structured np.array
            A structured np.array with dtype=[('species', np.float),
            ('area', np.float)]. 
    
        '''
        try:
            len(area_list); area_list = np.array(area_list)
        except:
            area_list = np.array([area_list])
        area_list = np.concatenate((np.array([anchor_area]), area_list))
        output_array = np.empty(len(area_list), dtype=[('species', np.float),\
                                                     ('area', np.float)])
        output_array['area'] = area_list
        c = S / (anchor_area ** z)
        p_law = lambda x: c * (x ** z)
        output_array['species'] = p_law(area_list)
        return output_array

    def predict_sar(self, sad, S, a_list, ssad):
        '''
        A generic sar function the utilizes the relationship between the sad
        and the ssad to generate the sar.  Can take any combination of sad and
        ssad. 

        Parameters
        ----------
        sad : ndarray
            Species abundance distribution, should sum to 1 or nearly so. Support 
            must be >= 1 (ie, no P(0) at start)
        S : int or float
            Number of species in landscape
        a_list : list
            List of area fractions at which to calculate SAD
        ssad : ssad distribution object
            Spatial abundance distribution object distributions module. N and a
            parameters are filled in the function.  Any additional parameters
            to a ssad distribution object need to be filled before the object
            is passed. See examples 

        Notes
        -----
        Example
        
        #Fill k parameter before passing
        sar1 = SAR().predict_sar(sad, S, [.1,.2,.3,.4,.5], nbd(k=.3))

        #No parameters to fill in tgeo. N and a are filled in function
        sar2 = SAR().predict_sar(sad, S, [.1,.2,.6], tgeo())

        '''
        sar = []
        N_range = np.arange(1, len(sad) + 1)
        for i, a in enumerate(a_list):
            assert a < 1, "a must be less than 1"
            p_pres_list = []
            for n in N_range:
                ssad.params['N'] = n; ssad.params['a'] = a
                p_pres_list.append(1 - ssad.pmf(0)[0])
            sar.append(sum(S * sad * np.array(p_pres_list)))
        return np.array(sar)

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

def canonical_lognorm_pmf(r, S, param_ret=False):
    '''
    canonical_lognorm_pmf(r, S, param_ret=False)

    Preston's canonical lognormal

    Parameters
    ----------
    r : int or array-like object
        The number of octaves the octave of interest is away from the the modal
        octave (see May 1975 or Preston 1948)
    S : int
        Total number of species in landscape

    Returns
    -------
    : np.ndarray
        If param_ret == False, returns an np.array with the probability that a
        given species is in the octave r units from the modal octave.

    Notes
    -----
    The canonical lognormal distribution is a one-parameter distribution that
    depends only on the total number of species in a landscape. The modal
    octave of this distribution is defined as the log2 abundance octave in
    which the most species fall.  This was parameterized using May 1975.


    '''
    #Going to try to change this assumption so that it can work for S <= 12
    assert S > 12, "S must be greater than 12"
    try:
        len(r); r = np.array(r)
    except:
        r = np.array([r])

    pisr = m.pi ** 0.5
    #NOTE: Equation from May 1975.  Approximation is made so works only with S
    #greater than 12
    eq = lambda s: ((np.log(2) / (2 * ((s * pisr) / S))) * np.log(s) ** 0.5)\
                                                                            - 1
    s0 = scipy.optimize.brentq(eq, 2, S, disp=True)
    a = (s0 * pisr) / S
    pmf = (s0 / S) * np.exp(-(a ** 2) * (r ** 2)) 
    return pmf, s0, a

def _ln_choose(n, k):
    '''
    Log binomial coefficient with extended gamma factorials. n and k may be int 
    or array - if both array, must be the same length.
    '''
    gammaln = scipy.special.gammaln
    return gammaln(n + 1) - (gammaln(k + 1) + gammaln(n - k + 1))

def _upscale_sar_(up_areas, N, S):
    '''
    This function is used to upscale from the anchor area.

    up_areas -- A 1D area of length upscale + 1.  The areas to be upscaled to.

    N -- Number of individuals at anchor scale (int)

    S -- Number of species at anchor scale (int)

    returns:
        1D array of species at a given upscaled area
    '''

    
    num_ind = np.empty(len(up_areas))
    spp = np.empty(len(up_areas))
    for i in xrange(len(up_areas)):
        if i == 0:
            num_ind[i] = N
            spp[i] = S
        else:
            num_ind[i] = 2 * num_ind[i - 1]
            N2A = num_ind[i]
            n = np.linspace(1, N2A, num=N2A)
            eq1 = lambda x: (sum((x**n)/n) * ((N2A) / ((x - x**(N2A + 1)) / \
                            (1 - x)) * (1 / x))) - ((N2A) * ((1 - x) / \
                            (x - x**(N2A + 1))) * (1 - (x**N2A / (N2A + 1))))\
                            - spp[i - 1]
            x = scipy.optimize.brentq(eq1, 1e-10, 1 - 1e-10, disp=True)
            eq2 = lambda S2A: np.log(1 / np.log(1 / x)) * ((N2A) / \
                                          ((x - x**(N2A + 1)) / (1 - x))) - S2A
            S2A = scipy.optimize.brentq(eq2, spp[i - 1], 2 * spp[i - 1],\
                                                                     disp=True)
            spp[i] = S2A
    return spp

def _downscale_sar_(down_areas, N, S):
    '''
    This function is used to downscale from the anchor area.

    up_areas -- A 1D area of length downscale + 1.  The areas to be downscaled
                to.

    N -- Number of individuals at anchor scale (int)

    S -- Number of species at anchor scale (int)

    returns:
        1D array of species at a given downscaled areas
    '''

    num_ind = np.empty(len(down_areas))
    spp = np.empty(len(down_areas))
    for i in xrange(len(down_areas)):
        if i == 0:
            num_ind[i] = N
            spp[i] = S
        else:
            num_ind[i] = 0.5 * num_ind[i - 1]
            if num_ind[i] <= 1:
                raise DownscaleError('Cannot downscale %i iterations from ' +\
                                     'anchor scale. One or less individuals' +\
                                     ' per cell.' % (len(down_areas) - 1))
            N_A = num_ind[i - 1]
            S_A = spp[i - 1]
            n = np.linspace(1, N_A, num=N_A)
            eq1 = lambda x: sum((x**n) / n) - ((S_A / N_A) * sum(x**n))
            x = scipy.optimize.brentq(eq1, 1e-10, min((sys.float_info[0] / S)\
                                      **(1/float(N)), 2), disp=True)
            ShalfA = (S_A * (1 / x)) - ((N_A) * ((1 - x) / (x - x**(N_A + 1)))\
                     * (1 - ((x**N_A) / (N_A + 1))))
            spp[i] = ShalfA
    return spp[::-1]



def _generate_areas_(anchor_area, upscale, downscale):
    '''
    Utility function that makes the area list
    
    '''

    areas = np.empty(upscale + downscale + 1)
    areas[downscale] = anchor_area
    for i in range(downscale)[::-1]:
        areas[i] = areas[i + 1] / 2
    for i in xrange(downscale + 1, len(areas)):
        areas[i] = areas[i - 1] * 2

    return areas









        






