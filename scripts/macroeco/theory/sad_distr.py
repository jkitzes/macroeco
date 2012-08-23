#!/usr/bin/python

'''
Species abundance distributions.

Distributions
-------------
- `lgsr` -- Fisher's log series (Fisher et al. 1943)
- `neg_binom` -- Negative binomial
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

Misc Functions
--------------
- `make_rank_abund` -- convert any SAD pmf into a rank abundance curve

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

Sugihara, G. 1980. Minimal community structure: an explanation of species
abundance patterns. The American Naturalist, 116:770-787.

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

#TODO: Add truncated negative binomial and truncated log-normal

class RootError(Exception):
    '''Error if no root or multiple roots exist for the equation generated
    for specified values of S and N'''

    def __init__(self, value=None):
        Exception.__init__(self)
        self.value = value
    def __str__(self):
        return '%s' % self.value

class lgsr:
    '''
    Fisher's log series distribution (Fisher et al. 1943, Hubbel 2001)

    Methods
    -------
    
    pmf(n, x)
        Probability mass function
    pmf_SN(n, S, N, param_out=False)
        Probability mass function generated with S and N
    cdf(n, S, N)
        Cumulative distribution function
    rad(S, N)
        Rank abundance distribution

    '''

    def pmf(self, n, S, N, param_out=False):
        '''
        Probability mass function of Fisher log series generated with S and N

        Parameters
        ----------
        n : int, float or array-like object
            Abundances at which to calculate the pmf
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

        assert S < N, "S must be less than N"
        assert S > 1, "S must be greater than 1"
        assert N > 0, "N must be greater than 0"
        try:
            len(n)
            n = np.array(n)
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
    
    def cdf(self, n, S, N):
        '''
        Cumulative distribution function for the Fisher logseries

        Parameters
        ----------
        n : int, float or array-like object
            Abundances at which to calculate the pmf
        S : int
            Total number of species in landscape
        N : int
            Total number of indviduals in landscape

        Returns
        -------
        : ndarray (1D)
            Returns array with cdf values for the given values of n.

        '''

        try:
            len(n)
            n = np.array(n)
        except:
            n = np.array([n])
        max_n = np.max(n)
        cdf = np.cumsum(self.pmf(np.arange(1, max_n + 1), S, N))
        return np.array([cdf[x - 1] for x in n])

    def rad(self, S, N):
        '''
        Rank abundance distribution for the Fisher logseries

        Parameters
        ----------
        S : int
            Total number of species in landscape
        N : int
            Total number of indviduals in landscape
        
        Returns
        -------
        : ndarray (1D)
            Returns and array of length S with the expected abundances given
            the Fisher log series
            
        '''

        full_pmf = self.pmf(np.arange(1, N + 1), S, N)
        return make_rank_abund(full_pmf, S)                                    

class neg_binom:
    '''
    Negative binomial distribution with infinite support.  Because this
    distribution has infinite support the fit the cdf fron 1:N will be less
    than one.

    Methods
    -------
    pmf(n, S, N, k, param_out=False)
        Probability mass function
    cdf(n, S, N, k):
        Cumulative distribution function
    rad(S, N, k)
        Rank abundance distribution
    fit(sad, guess_for_k=1)
        Maximum likelihood fitting for neg_binom

    '''

    def pmf(self, n, S, N, k, param_out=False):
        '''
        Probablity mass function of the Negative Bionomial generated with S and
        N

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
            param_out = True, returns the array as well as the parameter 
            estimates.

        '''
        assert S < N, "S must be less than N"
        assert S > 1, "S must be greater than 1"
        assert N > 0, "N must be greater than 0"
        try:
            len(n)
            n = np.array(n)
        except:
            n = np.array([n])
        assert np.max(n) <= N, "Maximum n must be less than or equal to N"

        mu = float(N) / S
        p = 1 / (mu / k + 1)  # See Bolker book Chapt 4
        pmf = stats.nbinom.pmf(n, k, p)

        if param_out == True:
            return (pmf, {'k' : k, 'p' : p})
        else:
            return pmf

    def cdf(self, n, S, N, k):
        '''
        Cumulative distribution function for the Negative Binomial

        Parameters
        ----------
        n : int, float or array-like object
            Abundances at which to calculate the pmf
        S : int
            Total number of species in landscape
        N : int
            Total number of indviduals in landscape

        Returns
        -------
        : ndarray (1D)
            Returns array with cdf values for the given values of n.

        '''
        try:
            len(n)
            n = np.array(n)
        except:
            n = np.array([n])
        max_n = np.max(n)
        cdf = np.cumsum(self.pmf(np.arange(1, max_n + 1), S, N, k))
        return np.array([cdf[x - 1] for x in n])

    def rad(self, S, N, k):
        '''
        Rank abundance distribution for the Infinite Negative Binomial

        Parameters
        ----------
        S : int
            Total number of species in landscape
        N : int
            Total number of indviduals in landscape
        k : float
            Aggregation paramter
        
        Returns
        -------
        : ndarray (1D)
            Returns and array of length S with the expected abundances given by
            the Negative Binomial

        '''
        
        full_pmf = self.pmf(np.arange(1, N + 1), S, N, k)
        return make_rank_abund(full_pmf, S)

    def fit(self, sad, guess_for_k=1):
        '''
        Fits the Infinite Negative Bionial to the given sad

        Parameters
        ----------
        sad : array-like object
            Observed species abundances 
        guess_for_k : float
            Default parameter for the approximate k given the data

        Returns
        -------
        : dict
            The maximum likelihood estimator (MLE) for k and the negative
            log-likelihood (NLL) of the negative binomial with MLE k. 

        '''
        def nll_nb(k, sad):
            return -sum(np.log(self.pmf(sad, len(sad), np.sum(sad), k)))
        mlek = scipy.optimize.fmin(nll_nb, np.array([guess_for_k]), args=\
                                   (sad,), disp=0)[0]
        return {'k' : mlek}

class geo_series:
    '''
    Geometric series distribution (Motomura 1932 and Magurran 1988).

    Methods
    -------
    rad(S, N, k)
        Rank abundance distribution

    '''

    #TODO:  Need to derive the pmf and cdf for the Geometric series
    def rad(self, S, N, k):
        '''
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
            The fraction of resources that each species acquires. Range is 
            (0, 1].

        Returns
        -------
        : ndarray (1D)
            An array containing the the predicted SAD with element 0 containing
            thespecies with the most individuals and element S - 1 containing 
            the species with the least individuals.
        '''
        assert S < N, "S must be less than N"
        assert S > 1, "S must be greater than 1"
        assert N > 0, "N must be greater than 0"
        assert k > 0 and k <= 1, "k must be between on the interval (0, 1]"

        C = (1 - (1 - k )** S) ** - 1
        pred_rank_abund = N * C * k * (1 - k) ** (np.arange(1, S + 1) - 1)
        return pred_rank_abund

    def fit_k(self, sad):
        '''
        Fit k parameter of geometric series (May 1975)

        Parameters
        ----------
        sad : array-like object
            The observed SAD to be fit to a geometric series

        Returns
        : float
            An estimate for k given the observed SAD
        '''
        sad = np.array(sad)
        S = len(sad)
        N = np.sum(sad)
        Nmin = np.min(sad)
        #Equation from May (1975)
        eq = lambda x: ((x / (1 - x)) * ((1 - x) ** S / (1 - (1 - x) ** S))) -\
                       (Nmin / N)
        k = scipy.optimize.brentq(eq, 1e-10, 1 - 1e-10, disp=True)
        return k

class broken_stick:
    '''
    McArthur's broken stick species abundance distribution (May 1975)

    Methods
    -------
    pmf(n, S, N, param_out=False)
        Probability mass function
    cdf(n, S, N)
        Cumulative distribution function
    rad(S, N)
        Rank abundance distribution

    '''

    def pmf(self, n, S, N, param_out=False):
        '''
        Probability mass function for MacArthur's broken stick

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
            param_out = True, returns a tuple with the pmf as an array as well 
            as a dictionary of the parameter estimates.
        '''
        assert S < N, "S must be less than N"
        assert S > 1, "S must be greater than 1"
        assert N > 0, "N must be greater than 0"
        try:
            len(n)
            n = np.array(n)
        except:
            n = np.array([n])
        assert np.max(n) <= N, "Maximum n must be less than or equal to N"

        eq = lambda x: ((S - 1) / N) * (1 - (x / N)) ** (S - 2)
        pmf = eq(n)

        if param_out == True:
            return (pmf, {'S' : S})
        else:
            return pmf

    def cdf(self, n, S, N):
        '''
        Cumulative distribution function for MacArthur's broken stick model

        Parameters
        ----------
        n : int, float or array-like object
            Abundances at which to calculate the pmf
        S : int
            Total number of species in landscape
        N : int
            Total number of indviduals in landscape

        Returns
        -------
        : ndarray (1D)
            Returns array with cdf values for the given values of n.

        '''
        try:
            len(n)
            n = np.array(n)
        except:
            n = np.array([n])
        max_n = np.max(n)
        cdf = np.cumsum(self.pmf(np.arange(1, max_n + 1), S, N))
        return np.array([cdf[x - 1] for x in n])
    
    def rad(self, S, N):
        '''
        Rank-abundance distribution for McArthur's broken-stick distribution
        (May 1975).

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
        assert S < N, "S must be less than N"
        assert S > 1, "S must be greater than 1"
        assert N > 0, "N must be greater than 0"

        pred_rank_abund = np.empty(S)
        for i in xrange(S):
            n = np.arange(i + 1, S + 1) 
            pred_rank_abund[i] = (N / S) * sum(1 / n)
        return pred_rank_abund

class plognorm:
    '''
    Poisson log-normal distribution (Bulmer 1974)

    Methods
    -------
    pmf(ab, mu, sigma, param_out=False)
        Probability mass function
    cdf(ab, mu, sigma)
        Cumulative distribution function
    rad(S, N, mu, sigma)
        Rank abundance distribution
    fit(sad)
        Maximum likelihood estimator for poisson lognormal

    '''
    
    def pmf(self, ab, mu, sigma, param_out=False):
        '''
        Probability mass function
        
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
            param_out = True, returns the array as well as the parameter 
            estimates.

        Notes
        -----
        This fuction was adopted directly from the VGAM package in R by Mark
        Wilber. The VGAM R package was adopted directly from Bulmer (1974).

        '''
        try:
            len(ab)
            ab = np.array(ab)
        except:
            ab = np.array([ab])
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
            #NOTE: How many parameters should I be returning, 1 or 2?
            return (pmf, {'mu' : mu, 'sigma' : sigma})
        else:
            return pmf

    def cdf(self, ab, mu, sigma):
        '''
        Cumulative distribution function

        Parameters
        ----------
        ab : int, float or array-like object
            Abundances at which to calculate the cdf 
        mu : float
            the mu parameter of the poisson log normal
        sigma : float
            the sigma parameter of the poisson log normal

        Returns
        -------
        : ndarray (1D)
            Returns array with cdf values for the given values of ab.

        Notes
        -----
        Any lognormal distribution is continuous which makes is slightly
        problematic when modeling SADs which are discrete (see Blackburn and 
        Gaston article). We condsider the plognorm to be discrete and sum to 
        find the cdf rather than integrating.

        '''

        try:
            len(ab)
            ab = np.array(ab)
        except:
            ab = np.array([ab])
        max_ab = np.max(ab)
        cdf = np.cumsum(self.pmf(np.arange(1, max_ab + 1), mu, sigma))
        return np.array([cdf[x - 1] for x in ab])

    def rad(self, S, N, mu, sigma):
        '''
        Rank abundance distribution

        Parameters
        ----------
        S : int
            Total number of species in landscape
        N : int
            Total number of indviduals in landscape
        mu : float
            the mu parameter of the poisson log normal
        sigma : float
            the sigma parameter of the poisson log normal

        Returns
        -------
        : ndarray (1D)
            Returns and array of length S with the expected abundances given
            by the poisson log-normal

        '''

        full_pmf = self.pmf(np.arange(1, N + 1), mu, sigma)
        return make_rank_abund(full_pmf, S)

    def fit(self, sad):
        '''
        Maximum likelihood Estimates for Poisson log normal

        Parameter
        ---------
        sad : array-like object
            The observed abundances which will be fit to a poison lognormal

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
            return -sum(np.log(self.pmf(sad, x[0], x[1])))
        mu, sigma = scipy.optimize.fmin(pln_func, x0 = [mu0, sigma0], disp=0)
        return { 'mu' : mu, 'sigma' : sigma}

class trun_plognorm:
    '''
    Truncated Poisson lognormal (Bulmer 1974)

    Methods
    -------
    pmf(ab, mu, sigma, param_out=False)
        Probability mass function
    cdf(ab, mu, sigma)
        Cumulative distribution function
    rad(S, N, mu, sigma)
        Rank abundance distribution
    fit(sad)
        Maximum likelihood estimator for distribution parameters

    '''

    def pmf(self, ab, mu, sigma, param_out=False):
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
            param_out = True, returns the array as well as a dict with
            parameter estimates.

        '''

        if param_out == True:
            untr_pmf, params = plognorm().pmf(ab, mu, sigma,\
                                                          param_out=param_out)
            pmf0 = plognorm().pmf(0, mu, sigma)
            tr_pmf = (untr_pmf / (1 - pmf0))#Truncating based on Bulmer eq. A1
            return (tr_pmf, {'mu' : params['mu'], 'sigma' : params['sigma']})
        else:
            untr_pmf = plognorm().pmf(ab, mu, sigma)
            pmf0 = plognorm().pmf(0, mu, sigma)
            tr_pmf = (untr_pmf / (1 - pmf0))
            return tr_pmf

    def cdf(self, ab, mu, sigma):
        '''
        Cumulative distribution function

        Parameters
        ----------
        ab : int, float or array-like object
            Abundances at which to calculate the cdf 
        mu : float
            the mu parameter of the poisson log normal
        sigma : float
            the sigma parameter of the poisson log normal

        Returns
        -------
        : ndarray (1D)
            Returns array with cdf values for the given values of ab.

        '''
        try:
            len(ab)
            ab = np.array(ab)
        except:
            ab = np.array([ab])
        max_ab = np.max(ab)
        cdf = np.cumsum(self.pmf(np.arange(1, max_ab + 1), mu, sigma))
        return np.array([cdf[x - 1] for x in ab])

    def rad(self, S, N, mu, sigma):
        '''
        Rank abundance distribution

        Parameters
        ----------
        S : int
            Total number of species in landscape
        N : int
            Total number of indviduals in landscape
        mu : float
            the mu parameter of the poisson log normal
        sigma : float
            the sigma parameter of the poisson log normal

        Returns
        -------
        : ndarray (1D)
            Returns and array of length S with the expected abundances given
            by the poisson log-normal

        '''

        full_pmf = self.pmf(np.arange(1, N + 1), mu, sigma)
        return make_rank_abund(full_pmf, S)

    def fit(self, sad):
        '''
        Maximum likelihood Estimates for truncated poisson log normal

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
            return -sum(np.log(self.pmf(sad, x[0], x[1])))
        mu, sigma = scipy.optimize.fmin(pln_func, x0 = [mu0, sigma0], disp=0)
        return { 'mu' : mu, 'sigma' : sigma}

class lognorm:
    '''
    Lognormal distribution

    Methods
    -------
    pmf(ab, mu, sigma, param_out=False)
        Probability mass function
    cdf(ab, mu, sigma, param_out=False)
        Cumulative distribution function
    rad(S, N, mu, sigma)
        Rank abundance distribution
    fit(sad)
        Maximum likelihood estimator for distribution parameters


    '''

    def pmf(self, ab, mu, sigma, param_out=False):
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
            param_out = True, returns the array as well as the parameter 
            estimates.

        Notes
        -----
        scipy.stats.lognorm is coded very poorly and the docstring is not 
        helpful so we coded our own lognormal distribution
        '''
        try:
            len(ab)
            ab = np.array(ab)
        except:
            ab = np.array([ab])
        pmf = stats.norm.pdf(np.log(ab), loc=mu, scale=sigma) / ab
        if param_out == True:
            return (pmf, {'mu' : mu, 'sigma' : sigma})
        else:
            return pmf

    def cdf(self, ab, mu, sigma, param_out=False):
        '''
        Cumulative distribution function

        Parameters
        ----------
        ab : int, float or array-like object
            Abundances at which to calculate the cdf 
        mu : float
            the mu parameter of the log normal
        sigma : float
            the sigma parameter of the log normal

        Returns
        -------
        : ndarray (1D)
            Returns array with cdf values for the given values of ab.

        '''

        try:
            len(ab)
            ab = np.array(ab)
        except:
            ab = np.array([ab])
        max_ab = np.max(ab)
        cdf = np.cumsum(self.pmf(np.arange(1, max_ab + 1), mu, sigma))
        return np.array([cdf[x - 1] for x in ab])


    def rad(self, S, N, mu, sigma):
        '''
        Rank abundance distribution

        Parameters
        ----------
        S : int
            Total number of species in landscape
        N : int
            Total number of indviduals in landscape
        mu : float
            the mu parameter of the log normal
        sigma : float
            the sigma parameter of the log normal

        Returns
        -------
        : ndarray (1D)
            Returns and array of length S with the expected abundances given
            by the log-normal

        '''
        full_pmf = self.pmf(np.arange(1, N + 1), mu, sigma)
        return make_rank_abund(full_pmf, S)
    
    def fit(self, sad):
        '''
        Maximum likelihood Estimates for log normal

        Parameter
        ---------
        sad : array-like object
            The observed abundances

        Returns
        -------
        : dict
            The maximum likelihood estimates of mu and sigma

        '''
        assert len(sad) >= 1, "len(sad) must be greater than or equal to 1"
        sad = np.array(sad)
        mu0 = np.mean(np.log(sad))
        sigma0 = np.std(np.log(sad), ddof=1)
        def ln_func(x):
            return -sum(np.log(self.pmf(sad, x[0], x[1])))
        mu, sigma = scipy.optimize.fmin(ln_func, x0 = [mu0, sigma0], disp=0)
        return {'mu' : mu, 'sigma' : sigma}


class sugihara:
    '''
    Sugihara Rank Abundance Distribution (Sugihara 1980)

    Methods
    -------
    rad(S, N, sample_size=10000)

    '''
    #TODO: Back-derive the pmf?

    def rad(self, S, N, sample_size=10000):
        '''
        Parameters
        ----------
        S : int
            Total number of species in the landscape
        N : int
            Total number of individuals in the lanscape

        Returns
        -------
        : 1D np.ndarray
            The predicted rank abundance distribution of the species in the
            landscape

        Notes
        -----
        This function is a bit questionable because we are using a randomly
        generated breakage sequence. As S gets large, we will start to under 
        sample some of the possible breakage sequences and this model begins to
        fail. Adapted from Sugihara (1980)

        '''
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

class mete_lgsr:
    '''
    METE truncated log series (Harte 2011)

    Methods
    -------
    pmf(n, S, N, param_out=False)
        Probability mass function
    cdf(n, S, N)
        Cumulative mass function 
    rad(S, N)
        Rank abundance distribution

    '''

    def pmf(self, n, S, N, param_out=False):
        '''
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
        are in the range (1/e, 1) (exclusive). Therefore, the start and stop 
        parameters for the brentq procedure are close to these values. However,
        x can occasionally be greater than one so the maximum stop value of the 
        brentq optimizer is 2.

        '''
        assert S < N, "S must be less than N"
        assert S > 1, "S must be greater than 1"
        assert N > 0, "N must be greater than 0"
        try:
            len(n)
            n = np.array(n)
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
    
    def cdf(self, n, S, N):
        '''
        Cumulative distribution function

        Parameters
        ----------
        n : int, float or array-like object
            Abundances at which to calculate the cdf 
        mu : float
            the mu parameter of the log normal
        sigma : float
            the sigma parameter of the log normal

        Returns
        -------
        : ndarray (1D)
            Returns array with cdf values for the given values of n.

        '''
        try:
            len(n)
            n = np.array(n)
        except:
            n = np.array([n])
        max_n = np.max(n)
        cdf = np.cumsum(self.pmf(np.arange(1, max_n + 1), S, N))
        return np.array([cdf[x - 1] for x in n])
        
    def rad(self, S, N):
        '''
        Rank abundance distribution for the METE truncated log series

        Parameters
        ----------
        S : int
            Total number of species in landscape
        N : int
            Total number of indviduals in landscape
        
        Returns
        -------
        : ndarray (1D)
            Returns and array of length S with the expected abundances given
            by the METE truncated log series
        '''

        full_pmf = self.pmf(np.arange(1, N + 1), S, N)
        return make_rank_abund(full_pmf, S)

class mete_lgsr_approx:
    '''
    METE log series distribution with approximation (Harte 2011)

    Methods
    -------
    pmf(n, S, N, param_out=False), root=2)
        Probability mass function
    cdf(n, S, N, root=2)
        Cumulative mass function 
    rad(S, N, root=2)
        Rank abundance distribution

    '''

    def pmf(self, n, S, N, param_out=False, root=2):
        '''
        Truncated log series using approximation (7.30) and (7.32) in Harte 
        2011

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
        assert S < N, "S must be less than N"
        assert S > 1, "S must be greater than 1"
        assert N > 0, "N must be greater than 0"
        try:
            len(n)
            n = np.array(n)
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

    def cdf(self, n, S, N, root=2):
        '''
        Cumulative distribution function

        Parameters
        ----------
        n : int, float or array-like object
            Abundances at which to calculate the cdf 
        mu : float
            the mu parameter of the log normal
        sigma : float
            the sigma parameter of the log normal

        Returns
        -------
        : ndarray (1D)
            Returns array with cdf values for the given values of n.

        '''
        try:
            len(n)
            n = np.array(n)
        except:
            n = np.array([n])
        max_n = np.max(n)
        cdf = np.cumsum(self.pmf(np.arange(1, max_n + 1), S, N, root=root))
        return np.array([cdf[x - 1] for x in n])

    def rad(self, S, N, root=2):
        '''
        Rank abundance distribution for the METE truncated log series with
        approximation

        Parameters
        ----------
        S : int
            Total number of species in landscape
        N : int
            Total number of indviduals in landscape
        
        Returns
        -------
        : ndarray (1D)
            Returns and array of length S with the expected abundances given
            by the METE truncated log series with approximation

        '''

        full_pmf = self.pmf(np.arange(1, N + 1), S, N, root=root)
        return make_rank_abund(full_pmf, S)

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
        len(r)
        r = np.array(r)
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






        






