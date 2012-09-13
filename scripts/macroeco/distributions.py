#!/usr/bin/python

'''
Macroecological distributions and curves.

Distributions
-------------

SAD
- `logser` -- Fisher's log series (Fisher et al. 1943) (ie, logarithmic)
- `logser_ut` -- Upper-truncated log series (Harte 2011)
- `geo_series' -- Geometric series distribution (Motomura 1932)
- `broken_stick` -- McArthur's broken stick distribution (May 1975)
- `plognorm` -- Poisson log-normal (Bulmer 1974)
- `plognorm_lt` -- Truncated poisson log-normal (Bulmer 1974)
- `lognorm` -- Lognormal distribution
- `canonical_lognormal_pmf` -- Preston's canonical lognormal parameterized by
  May (1975)
- `sugihara` -- Sugihara's sequential breakage model (Sugihara 1980)
- `logser_ut_appx` -- METE log series using approximation (Harte 2011)
- `nbd_lt` - Lower truncated negative binomial

SAR
- `mete_sar_iter` - METE sar functions (Harte 2011)
- `powerlaw` - Power law sar
- `gen_sar` - A generic sar that supports any sad and ssad distribution.  Six
  unique sars are derived classes from gen_sar.  More can be made.

SSAD
- `binm` - Binomial distribution (Random Placement Model)
- `pois` - Poisson distribution
- `nbd` - Negative binomial distribution
- `fnbd` - Finite negative binomial (Zillio and He 2010)
- `geo` - Geometric distribution
- `fgeo` - Finite geometric distribution (Zillio and He 2010)
- `tgeo` - Truncated geometric distribution (Harte et al. 2008)

Misc Functions
--------------
- `make_array` 
- `make_rank_abund` 
- `_ln_choose`
- `_downscale_sar_`
- `_upscale_sar_`
- `_generate_areas_`
- `expand_n`

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
from docinherit import DocInherit
#from macroeco.utils.docinherit import DocInherit

doc_inherit = DocInherit

__author__ = "Justin Kitzes and Mark Wilber"
__copyright__ = "Copyright 2012, Regents of the University of California"
__credits__ = ["John Harte"]
__license__ = None
__version__ = "0.1"
__maintainer__ = "Justin Kitzes and Mark Wilber"
__email__ = "jkitzes@berkeley.edu"
__status__ = "Development"

# TODO: Add truncated log-normal?

# TODO: For all subclass inits - what to do if fit method later tries to
# overwrite these? Error, do it with warning?


# ----------------------------------------------------------------------------
# Define base classes Curve, Distribution, and custom errors
# ----------------------------------------------------------------------------


class Curve(object):
    '''
    Class for curves. Curves are quantitative predictions with no requirement 
    to sum to one, such as SAR's.
    
    Attributes
    ----------
    params : dict
        Contains all keyword arguments used to initialize distribution object

    Methods
    -------
    val(n)
        Value of curve at given points
    fit(data)
        Uses data to populate params attribute
    '''

    def __init__(self, **kwargs):
        '''
        Initialize distribution object.

        Stores keyword arguments as params attribute and defines minimum 
        support for distribution.

        Parameters
        ----------
        kwargs : comma separate list of keyword arguments
            Parameters needed for distribution. Stored in dict as params 
            attribute.

        See class docstring for more specific information on this distribution.
        '''
        # This method does nothing, but exists so that derived class pmf
        # methods can inherit this docstring.
        pass

    def vals(self, n):
        '''
        Value of curve method.

        Each derived class requires certain parameters to be available in the 
        params attribute. See class docstring for list.

        Parameters
        ----------
        n : int, float or array-like object
            Values at which to calculate curve. Note that only one list of 
            points may be given, not a list, no matter how long other required 
            parameter lists are.

        Returns
        -------
        vals : list of ndarrays
            List of 1D arrays of value of curve at points n.
        vars : dict containing lists of floats
            Intermediate parameter variables calculated for pmf, as 
            dict.

        See class docstring for more specific information on this distribution.
        '''
        # This method does nothing, but exists so that derived class pmf
        # methods can inherit this docstring.
        pass

    def fit(self, data, full_sad):
        '''
        This fit method fills the required parameters for an SARCurve object.

        Parameters
        ----------
        data :  tuple of array-like objects
            data conatains two array-like object.  The first is a list of area
            fractions (area / anchor area) and the second is a list of
            species/items numbers corresponding to those area fractions.
        full_sad : array-like object
            The full_sad at the anchor area

        '''
        
        #Assertion statements
        assert len(data) == 2, "data must contain two objects"
        assert len(data[0]) == len(data[1]), "Area and species number " + \
                                        "arrays must be of the same length"
        full_sad = make_array(full_sad)
        self.params['S'] = len(full_sad)
        self.params['N'] = sum(full_sad)
        return self


class Distribution(object):
    '''
    Class for statistical distributions.

    Attributes
    ----------
    params : dict
        Contains all keyword arguments used to initialize distribution object
    min_supp : int
        Minimum support for distribution, usually 0 or 1
    par_num : int
        Number of free parameters of distribution, used for AIC calculations

    Methods
    -------
    pmf(n)
        Probability mass function
    cdf(n)
        Cumulative distribution function
    rad()
        Rank abundance distribution, calculated from cdf
    fit(data)
        Uses data to populate params attribute
    '''

    def __init__(self):
        '''
        Initialize distribution object.

        Stores keyword arguments as params attribute and defines minimum 
        support for distribution.

        Parameters
        ----------
        kwargs : comma separate list of keyword arguments
            Parameters needed for distribution. Stored in dict as params 
            attribute.

        See class docstring for more specific information on this distribution.
        '''
        # This method does nothing, but exists so that derived class init
        # methods can inherit this docstring.
        pass

        # TODO: Check that all params are same length

    def pmf(self, n):
        '''
        Probability mass function method.

        Parameters
        ----------
        n : int, float or array-like object
            Values at which to calculate pmf. May be a list of same length as 
            parameters, or single iterable.

        Returns
        -------
        pmf : list of ndarrays
            List of 1D arrays of probability of observing sample n.
        vars : dict containing lists of floats
            Intermediate parameter variables calculated for pmf, as 
            dict.

        See class docstring for more specific information on this distribution.
        '''
        # This method does nothing, but exists so that derived class pmf
        # methods can inherit this docstring.
        pass


    def cdf(self, n):
        '''
        Cumulative distribution method.  

        Parameters
        ----------
        n : int, float or array-like object
            Values at which to calculate cdf. May be a list of same length as 
            parameters, or single iterable.

        Returns
        -------
        cdf : list of ndarrays
            List of 1D arrays of probability of observing sample n.
        vars : dict containing lists of floats
            Intermediate parameter variables calculated for cdf, as dict.

        See class docstring for more specific information on this distribution.
        '''
    
        for kw in self.params.iterkeys():
            if not np.iterable(self.params[kw]):
                self.params[kw] = make_array(self.params[kw])

        # Expand n argument if needed, assumes all params same length
        n = expand_n(n, len(self.params.values()[0]))

        # Calculate pmfs
        max_n = [np.max(tn) for tn in n]
        n_in = [np.arange(self.min_supp, i + 1) for i in max_n]
        pmf_list, var = self.pmf(n_in)

        # Calculate cdfs
        cdf = []
        for tpmf, tn in zip(pmf_list, n):
            full_cdf = np.cumsum(tpmf)
            tcdf = np.array([full_cdf[x - self.min_supp] for x in tn])
            cdf.append(tcdf)

        return cdf, var


    def rad(self):
        '''
        Rank abundance distribution method, calculates rad using pmf.

        Returns
        -------
        rad : list of ndarrays
            List of 1D arrays of predicted abundance for each species

        See class docstring for more specific information on this distribution.
        '''

        # Get parameters
        S, N = self.get_params()

        # TODO: Add error or warning if N too large for memory
        # FIX: Python throws a MemoryError when it runs out of RAM.

        # Calculate pmfs, going up to N for upper limit
        n_arrays = [np.arange(self.min_supp, 1*(i + 1)) for i in N]
        pmf, var = self.pmf(n_arrays)
        
        # Calculate rad
        rad = []
        for tS, tN, tpmf in zip(S, N, pmf):
            trad = make_rank_abund(tpmf, tS, min_supp=self.min_supp)
            rad.append(trad)

        return rad


    def fit(self, data):
        '''
        Fit method.

        Uses input data to get best fit parameters for distribution, and stores 
        these parameters in params attribute.
        
        Parameters
        ----------
        data : list of ndarrays
            Data to use to fit parameters of distribution. Even if only one 
            data array, must be in a list with one element.

        See class docstring for more specific information on this distribution.
        '''

        # By default, loop through ndarrays in data and extract n_samp
        # and tot_obs for each one.

        #Check data argument
        if type(data) != type([]):
            raise TypeError('Data must be a list of iterables')
        if not np.all(np.array([np.iterable(dt) for dt in data])):
            raise TypeError('Objects in data must be iterable')
        
        #Convert data argument to arrays
        data = [np.array(data) for data in data]

        n_samp = []
        tot_obs = []
        
        for tdata in data:
            n_samp.append(len(tdata))
            tot_obs.append(np.sum(tdata))

        self.params['n_samp'] = n_samp
        self.params['tot_obs'] = tot_obs

        return self


    def get_params(self, n=None):
        '''
        Gets and validates basic distribution parameters

        Parameters
        ----------
        n : int, float or array-like object
            Values at which to calculate pmf. May be a list of same length as 
            parameters, or single iterable.

        Returns
        -------
        : tuple
            Validated n_samp/S, tot_obs/N, and n parameters
 
        '''
        n_samp = make_array(self.params.get('n_samp', None))
        if n_samp[0] is None:
            n_samp = make_array(self.params.get('S', None))

        tot_obs = make_array(self.params.get('tot_obs', None))
        if tot_obs[0] is None:
            tot_obs = make_array(self.params.get('N', None))

        if n != None:
            n = expand_n(n, len(n_samp))

        # Validate parameters
        assert len(n_samp) == len(tot_obs), 'Length of n_samp/S and' +\
                                            ' tot_obs/N must be ' + 'the same'
        assert n_samp[0] != None, 'n_samp/S parameter not given'
        assert tot_obs[0] != None, 'tot_obs/N parameter not given'
        assert np.all(n_samp > 1), 'n_samp/S must be greater than 1'
        assert np.all(tot_obs > 0), 'tot_obs/N must be greater than 0'
    
        if n != None:
            return n_samp, tot_obs, n
        else:
            return n_samp, tot_obs


class RootError(Exception):
    '''Error if no or multiple roots exist when only one should exist.'''

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


# ----------------------------------------------------------------------------
# Distributions
# ----------------------------------------------------------------------------


class logser(Distribution):
    __doc__ = Distribution.__doc__ + \
    '''
    Description
    -----------
    Fisher's log series distribution (Fisher et al. 1943). Also known as 
    logarithmic distribution.

    Parameters
    ----------
    S or n_samp : int or iterable
        Total number of species / samples
    N or tot_obs: int or iterable
        Total number of individuals / observations

    Vars
    ----
    p : list of floats
        p parameter of standard logseries distribution

    Notes
    -----
    To use a known mean of the distribution as the parameter, set n_samp = 1 
    and tot_obs = mean.
    '''
    
    @doc_inherit
    def __init__(self, **kwargs):
        self.params = kwargs
        self.min_supp = 1
        self.par_num = 1
    
    @doc_inherit
    def pmf(self, n):
        
        # Get parameters
        S, N, n = self.get_params(n=n)
        assert np.all(S < N), 'n_samp/S must be less than tot_obs/N'
        
        # Calculate pmf
        stop = 1 - 1e-10
        start = -2
        eq = lambda p, S, N: (((N/p) - N) * (-(np.log(1 - p)))) - S

        pmf = []
        var = {}
        var['p'] = []

        for tS, tN, tn in zip(S, N, n):
            tp = scipy.optimize.brentq(eq, start, stop, args=(tS,tN), 
                                       disp=True)
            tpmf = stats.logser.pmf(tn, tp)
            var['p'].append(tp)
            pmf.append(tpmf)
   
        return pmf, var

    # TODO: Add custom method for cdf based on equation
    

class logser_ut(Distribution):
    __doc__ = Distribution.__doc__ + \
    '''
    Description
    -----------
    Upper-truncated log series (Harte et al 2008, Harte 2011). Like Fisher's, 
    but normalized to finite support.

    Parameters
    ----------
    S or n_samp : int or iterable
        Total number of species / samples
    N or tot_obs: int or iterable
        Total number of individuals / observations
   
    Vars
    ----
    x : list of floats
         exp(-beta) parameter

    Notes
    -----
    This distribution is the truncated logseries described in Eq 7.32 of Harte 
    2011. Eq. 7.27 is used to solve for the Lagrange multiplier.

    Realistic values of x where x = e**(-beta) are in the range (1/e, 1). The 
    start and stop parameters for the brentq procedure are set close to these 
    values. However, x can occasionally be greater than one, so the maximum 
    stop value of the brentq optimizer is 2.
    '''
    
    @doc_inherit
    def __init__(self, **kwargs):
        self.params = kwargs
        self.min_supp = 1
        self.par_num = 0 #Number of parameters

    @doc_inherit    
    def pmf(self, n):

        # Get parameters
        S, N, n = self.get_params(n=n)
        assert np.all(S < N), 'n_samp/S must be less than tot_obs/N'

        # Calculate pmf
        start = 0.3
        stop = 2
        flmax = sys.float_info[0]
        eq = lambda x,k,N,S: sum(x ** k / float(N) * S) -  sum((x ** k) / k)

        pmf = []
        var = {}
        var['x'] = []

        for tS, tN, tn in zip(S, N, n):
            k = np.linspace(1, tN, num=tN)
            tx = scipy.optimize.brentq(eq, start,
                                       min((flmax/tS)**(1/float(tN)), stop), 
                                       args = (k, tN, tS), disp=True)
            tnorm = np.sum(tx ** k / k)  # From Ethan White's trun_logser_pmf
            tpmf = (tx ** tn / tn) / tnorm
            var['x'].append(tx)
            pmf.append(tpmf)
   
        return pmf, var

    # TODO: Add exact cdf from JK dissertation
    # FIX: For saving computational time?  Cumulative sums of pmfs are the
    # exact cdfs...need to check that the analytical cdf is faster.


class logser_ut_appx(Distribution):
    __doc__ = Distribution.__doc__ + \
    '''
    Description
    -----------
    Upper-truncated log series (Harte et al 2008, Harte 2011). Like Fisher's, 
    but normalized to finite support.  Uses approximation from Harte (2011).

    Parameters
    ----------
    S or n_samp : int or iterable
        Total number of species / samples
    N or tot_obs: int or iterable
        Total number of individuals / observations
   
    Vars
    ----
    x : list of floats
         exp(-beta) parameter
        
    Notes:
    ------
    This distribution is the truncated logseries described in Eq 7.32 of Harte 
    2011. Eq. 7.30, the approximte equation, is used to solve for the Lagrange 
    multiplier. This class is faster than the logser with no approximation.

    Realistic values of x where x = e**(-beta) are in the range (1/e, 1). The 
    start and stop parameters for the brentq procedure are set close to these 
    values. However, x can occasionally be greater than one, so the maximum 
    stop value of the brentq optimizer is 2.

    The pmf method has an additional optional keyword argument 'roots'. In the 
    approximation equation, there are two roots (solutions) in the solution for 
    the lagrange multiplier.  Root 2 is the root typically used in calculations 
    and is the default.  If root=1, the first root will be used and this is not 
    a true pmf.
    '''
    
    @doc_inherit
    def __init__(self, **kwargs):
        self.params = kwargs
        self.min_supp = 1
        self.par_num = 0
    
    @doc_inherit
    def pmf(self, n, root=2):

        # Get parameters
        S, N, n = self.get_params(n=n)
        assert np.all(S < N), 'n_samp/S must be less than tot_obs/N'

        # Calculate pmf
        start = 0.3
        stop = 1 - 1e-10
        eq = lambda x, S, N: (((-m.log(x))*(m.log(-1/(m.log(x))))) - 
                                                                (float(S)/N))
        pmf = []
        var = {}
        var['x'] = []

        for tS, tN, tn in zip(S, N, n):
            
            # First try using brentq solver
            try:
                x = scipy.optimize.brentq(eq, start, stop, args=(tS, tN),
                                                                    disp=True)
            # If that fails, try a more complex decision tree
            except ValueError:
                eq1 = lambda x: -1 * eq(x, tS, tN)
                xmax = scipy.optimize.fmin(eq1, .5, disp=0)[0]
                ymax = eq(xmax, tS, tN)
                if ymax > 0:
                    if root == 1:
                        x = scipy.optimize.brentq(eq, start, xmax, args=(tS,
                                                                tN), disp=True)
                    if root == 2:
                        x = scipy.optimize.brentq(eq, xmax, stop, args=(tS,
                                                                tN), disp=True)
                if ymax < 0:
                    raise RootError('No solution to constraint equation with '
                                                  + 'given values of S and N') 

            g = -1/np.log(x)
            tpmf = (1/np.log(g)) * ((x**tn)/tn)

            var['x'].append(x)
            pmf.append(tpmf)

        return pmf, var


class plognorm(Distribution):
    __doc__ = Distribution.__doc__ + \
    '''
    Poisson log-normal distribution (Bulmer 1974)

    Parameters
    ----------
    mu : float
        The mu parameter of the poisson log normal
    sigma : float
        The sigma parameter of the poisson log normal
    S or n_samp : int or iterable (optional)
        Total number of species / samples
    N or tot_obs: int or iterable (optional)
        Total number of individuals / observations

    Vars
    ----
    None

    Notes
    -----
    The pmf method was adopted directly from the VGAM package in R by Mark
    Wilber. The VGAM R package was adopted directly from Bulmer (1974). The fit 
    function was adapted from Ethan White's pln_solver function in 
    weecology.
    '''
    
    # @doc_inherit cannot be used here because of derived plognorm_lt
    def __init__(self, **kwargs):
        '''
        Initialize distribution object.

        Stores keyword arguments as params attribute and defines minimum 
        support for distribution.

        Parameters
        ----------
        kwargs : comma separate list of keyword arguments
            Parameters needed for distribution. Stored in dict as params 
            attribute.

        See class docstring for more specific information on this distribution.
        '''
        self.params = kwargs
        self.min_supp = 0
        self.par_num = 2

    
    # @doc_inherit cannot be used here because of derived plognorm_lt
    def pmf(self, n):
        '''
        Probability mass function method.

        Parameters
        ----------
        n : int, float or array-like object
            Values at which to calculate pmf. May be a list of same length as 
            parameters, or single iterable.

        Returns
        -------
        pmf : list of ndarrays
            List of 1D arrays of probability of observing sample n.
        vars : dict containing lists of floats
            Intermediate parameter variables calculated for pmf, as 
            dict.

        See class docstring for more specific information on this distribution.
        '''

        # TODO: Allow mu to be calculated from S and N

        # Get parameters
        mu = make_array(self.params.get('mu', None))
        sigma = make_array(self.params.get('sigma',  None))

        # Validate parameters
        assert mu[0] != None, 'mu paramater not given'
        assert sigma[0] != None, 'sigma parameter not given'
        assert len(mu) == len(sigma), 'Length of mu and sigma must be the same'

        n = expand_n(n, len(mu))
        n_uniq = [np.unique(tn) for tn in n]  # Speed up by calc for uniq vals

        # Calculate pmf, no intermediate vars
        pmf = []

        eq = lambda t,x,mu,sigma: np.exp(t*x - np.exp(t) - 0.5*((t-mu) / 
                                                                sigma)**2)

        for tmu, tsigma, tn_uniq, tn in zip(mu, sigma, n_uniq, n):
            
            # If mu negative, pmf 0
            if tmu <= 0 or tsigma <= 0:
                tpmf_uniq = np.repeat(1e-120, len(tn_uniq))

            # Calculate unique pmf values
            else:
                # TODO: Throwing overflow warning but not affecting result
                tpmf_uniq = np.empty(len(tn_uniq), dtype=np.float)
                for i, g in enumerate(tn_uniq):
                    if g <= 170:
                        integ = integrate.quad(eq, -np.inf, np.inf, 
                                               args=(g,tmu,tsigma))[0]
                        norm = np.exp(-0.5 * m.log(2 * m.pi * tsigma**2) -
                                                m.lgamma(g + 1))
                        tpmf_uniq[i] = norm * integ
                    else:
                        z = (m.log(g) - tmu) / tsigma
                        tpmf_uniq[i] = ((1 + (z**2 + m.log(g) - tmu - 1) /
                                   (2 * g * tsigma**2)) * np.exp(-0.5 * z**2) /
                                    (m.sqrt(2 * m.pi) * tsigma * g))

            # Expand to full pmf
            tpmf = np.empty(len(tn))
            for i, g in enumerate(tn_uniq):
                tpmf[np.where(tn == g)[0]] = tpmf_uniq[i]

            pmf.append(tpmf)

        return pmf, None

    # TODO: Is there a known cdf?
    
    # @doc_inherit cannot be used here because of derived plognorm_lt
    def fit(self, data):
        '''
        Fit method.

        Uses input data to get best fit parameters for distribution, and stores 
        these parameters in params attribute.
        
        Parameters
        ----------
        data : list of ndarrays
            Data to use to fit parameters of distribution. Even if only one 
            data array, must be in a list with one element.

        See class docstring for more specific information on this distribution.
        '''

        # TODO: Check that data is a list of iterables
        # FIX: Check that data is list and objects in data are iterable 
        if type(data) != type([]):
            raise TypeError('Data must be a list of iterables')
        if not np.all(np.array([np.iterable(dt) for dt in data])):
            raise TypeError('Objects in data must be iterable')

        # Check data
        data = [np.array(data) for data in data]

        # Calculate and store parameters
        temp_mu = []
        temp_sigma = []
        self.params['n_samp'] = []
        self.params['tot_obs'] = []

        for tdata in data:
            mu0 = np.mean(np.log(tdata))  # Starting guesses for mu and sigma
            sigma0 = np.std(np.log(tdata), ddof=1)

            def pln_func(x):
                self.params['mu'] = x[0]
                self.params['sigma'] = x[1]
                return -sum(np.log(self.pmf(tdata)[0][0]))

            mu, sigma = scipy.optimize.fmin(pln_func, x0=[mu0, sigma0],
                                            disp=0)
            temp_mu.append(mu)
            temp_sigma.append(sigma)

            self.params['n_samp'].append(len(tdata))
            self.params['tot_obs'].append(np.sum(tdata))

        self.params['mu'] = temp_mu
        self.params['sigma'] = temp_sigma

        return self


class plognorm_lt(plognorm):
    __doc__ = Distribution.__doc__ + \
    '''
    Lower truncated Poisson lognormal (Bulmer 1974)

    Parameters
    ----------
    mu : float
        the mu parameter of the poisson log normal
    sigma : float
        the sigma parameter of the poisson log normal
    S or n_samp : int or iterable (optional)
        Total number of species / samples
    N or tot_obs: int or iterable (optional)
        Total number of individuals / observations

    Vars
    ----
    None

    Notes
    -----
    The pmf method was adopted directly from the VGAM package in R by Mark
    Wilber. The VGAM R package was adopted directly from Bulmer (1974). The fit 
    function was adapted from Ethan White's pln_solver function in weecology.

    Truncation calculation based on Bulmer Eq. A1.
    '''

    # @doc_inherit cannot be used here because class is derived from plognorm
    def __init__(self, **kwargs):
        '''
        Initialize distribution object.

        Stores keyword arguments as params attribute and defines minimum 
        support for distribution.

        Parameters
        ----------
        kwargs : comma separate list of keyword arguments
            Parameters needed for distribution. Stored in dict as params 
            attribute.

        See class docstring for more specific information on this distribution.
        '''
        self.params = kwargs
        self.min_supp = 1
        self.par_num = 2
        
    
    # @doc_inherit cannot be used here because class is derived from plognorm
    def pmf(self, n):
        '''
        Probability mass function method.

        Parameters
        ----------
        n : int, float or array-like object
            Values at which to calculate pmf. May be a list of same length as 
            parameters, or single iterable.

        Returns
        -------
        pmf : list of ndarrays
            List of 1D arrays of probability of observing sample n.
        vars : dict containing lists of floats
            Intermediate parameter variables calculated for pmf, as 
            dict.

        See class docstring for more specific information on this distribution.
        '''

        # Get parameters
        mu = make_array(self.params.get('mu', None))
        sigma = make_array(self.params.get('sigma',  None))

        # Validate parameters
        assert mu[0] != None, 'mu paramater not given'
        assert sigma[0] != None, 'sigma parameter not given'
        assert len(mu) == len(sigma), 'Length of mu and sigma must be the same'

        # Calculate pmf, using plognorm as aid
        reg_plog = plognorm(mu=mu, sigma=sigma)
        reg_pmf, reg_var = reg_plog.pmf(n)
        reg_pmf0, reg_var0 = reg_plog.pmf(0)

        trunc_pmf = [(pr / (1 - p0)) for pr, p0 in zip(reg_pmf, reg_pmf0)]

        return trunc_pmf, None 

    # TODO: Write cdf method based on cdf of plognorm, similar to above


class lognorm(Distribution):
    __doc__ = Distribution.__doc__ + \
    '''
    Description
    ------------
    Lognormal distribution

    Parameters
    ----------
    mu : float
        The mu parameter of the log normal
    sigma : float
        The sigma parameter of the log normal
    S or n_samp : int or iterable (optional)
        Total number of species / samples
    N or tot_obs: int or iterable (optional)
        Total number of individuals / observations

    Vars
    ----
    None
    
    Notes
    -----
    Log-normal distributions are continuous and therefore the cdf should be an
    integral. However, the cdf of an integral from a to a is 0 and the
    probability of there being a single individual within a species given an
    SAD is certainly not 0.  Therefore, we consider the lognormal "discrete"
    and calcuate the cdf by summing.  Note that this is one of the many 
    problems that Blackburn and Gaston have with using the lognormal for SADs.
    '''

    # TODO: MW - what about using exact equation for cdf/quantile function?
    # FIX:  To JK from MW: Refactored lognorm with scipy function for pdf and
    # cdf
    
    @doc_inherit
    def __init__(self, **kwargs):
        self.params = kwargs
        self.min_supp = 0
        self.par_num = 2

    
    @doc_inherit
    def pmf(self, n):

        # Get parameters
        mu = make_array(self.params.get('mu', None))
        sigma = make_array(self.params.get('sigma',  None))
        n = expand_n(n, len(mu))

        # Validate parameters
        assert mu[0] != None, 'mu paramater not given'
        assert sigma[0] != None, 'sigma parameter not given'
        assert len(mu) == len(sigma), 'Length of mu and sigma must be the same'

        # Calculate pmf
        pmf = []
        for tmu, tsigma, tn in zip(mu, sigma, n):
            # TODO: Why divided by tn?
            # FIX: In the way I had it coded previously, you had to divide by
            # tn to get the correct answer (See Wikipedia page of lognormal). I
            # have now coded it using scipy.
            tpmf = stats.lognorm.pdf(tn, tsigma, scale=np.exp(tmu))
            pmf.append(tpmf)

        return pmf, None

    @doc_inherit
    def cdf(self, n):

        # Get parameters
        mu = make_array(self.params.get('mu', None))
        sigma = make_array(self.params.get('sigma',  None))
        n = expand_n(n, len(mu))

        # Validate parameters
        assert mu[0] != None, 'mu paramater not given'
        assert sigma[0] != None, 'sigma parameter not given'
        assert len(mu) == len(sigma), 'Length of mu and sigma must be the same'

        #Calculate cdf
        cdf = []
        for tmu, tsigma, tn in zip(mu, sigma, n):
            tcdf = stats.lognorm.cdf(tn, tsigma, scale=np.exp(tmu))
            cdf.append(tcdf)

        return cdf, None

    
    @doc_inherit
    def fit(self, data):

        # TODO: Should this be done in all fit methods?
        # Fix:  Implemented checking of data argument in base class

        if type(data) != type([]):
            raise TypeError('Data must be a list of iterables')
        if not np.all(np.array([np.iterable(dt) for dt in data])):
            raise TypeError('Objects in data must be iterable')

        # Process data
        data = [np.array(data) for data in data]

        # Calculate and store parameters
        temp_mu = []
        temp_sigma = []
        self.params['n_samp'] = []
        self.params['tot_obs'] = []

        for tdata in data:
            mu0 = np.mean(np.log(tdata))  # Starting guesses for mu and sigma
            sigma0 = np.std(np.log(tdata), ddof=1)

            def ln_func(x):
                self.params['mu'] = x[0]
                self.params['sigma'] = x[1]
                return -sum(np.log(self.pmf(tdata)[0][0]))

            mu, sigma = scipy.optimize.fmin(ln_func, x0=[mu0, sigma0],
                                            disp=0)
            temp_mu.append(mu)
            temp_sigma.append(sigma)

            self.params['n_samp'].append(len(tdata))
            self.params['tot_obs'].append(np.sum(tdata))

        self.params['mu'] = temp_mu
        self.params['sigma'] = temp_sigma

        return self


class geo_ser(Distribution):
    __doc__ = Distribution.__doc__ + \
    '''
    Geometric series distribution (Motomura 1932 and Magurran 1988).

    Parameters
    ----------
    S or n_samp : int or iterable
        Total number of species / samples
    N or tot_obs: int or iterable
        Total number of individuals / observations
    k : float
        The fraction of resources that each species acquires. Range is 
        (0, 1].

    Var
    ----
    None
    
    Notes
    -----
    Equation for pmf and fit from May (1975).
    
    The support of this distribution is from [1, N] where N is the total number
    of individuals/observations. Therefore, the cdf of this function is one at
    N.  The empirical cdf of observed data reaches one at N - S + 1 (if not
    sooner) and therefore comparisons of geo_ser cdfs and empirical cdfs will
    often deviate near as n approaches N - S + 1 (if not sooner). 
    
    
    '''
    # TODO: @MW - Last sentence of docstring hard to understand - clarify.
    # FIX: Clarified!


    @doc_inherit
    def __init__(self, **kwargs):
        self.params = kwargs
        self.min_supp = 1
        self.par_num = 2


    @doc_inherit
    def pmf(self, n):

        # Get parameters
        S, N, n = self.get_params(n=n)
        k = make_array(self.params.get('k', None))
        assert k[0] != None, 'k parameter not given'
        assert np.all(k > 0) and np.all(k <= 1), ('k must be in the '
                                                  'interval (0, 1]')

        # Calculate pmf
        pmf = []
        eq = lambda x, k: (1 / x) * (1 / np.log(1 / (1 - k)))
        for tS, tN, tk, tn in zip(S, N, k, n):
            sumg = sum(eq(np.arange(1, tN + 1), tk))
            tpmf = eq(tn, tk) / sumg
            pmf.append(tpmf)

        return pmf, None


    @doc_inherit
    def rad(self):

        # Get parameters
        S, N = self.get_params()
        k = make_array(self.params.get('k', None))
        assert k[0] != None, "k parameter not given"
        assert np.all(k > 0) and np.all(k <= 1), ('k must be in the '
                                                  'interval (0, 1]')
        # Calculate rad
        rad = []
        for tS, tN, tk in zip(S, N, k):
            C = (1 - (1 - k )** S) ** - 1
            trad = N * C * k * (1 - k) ** (np.arange(1, S + 1) - 1)
            rad.append(trad)

        return rad


    @doc_inherit
    def fit(self, data):

        # Get parameters
        super(geo_ser, self).fit(data)  # Run Distribution.fit method
        S = self.params['n_samp']
        N = self.params['tot_obs']
        assert np.all(S < N), 'n_samp/S must be less than tot_obs/N'

        # Calculate fit
        self.params['k'] = []
        for tdata, tS, tN in zip(data, S, N):
            tNmin = np.min(tdata)
            eq = lambda x: (((x / (1 - x)) *
                             ((1 - x) ** tS / (1 - (1 - x) ** tS)))
                            - (tNmin / tN))
            tk = scipy.optimize.brentq(eq, 1e-10, 1 - 1e-10, disp=True)
            self.params['k'].append(tk)

        return self


class broken_stick(Distribution):
    __doc__ = Distribution.__doc__ + \
    '''
    McArthur's broken stick distribution (May 1975)

    Parameters
    ----------
    S or n_samp : int or iterable
        Total number of species / samples
    N or tot_obs: int or iterable
        Total number of individuals / observations

    Vars
    ----
    None
    '''

    @doc_inherit
    def __init__(self, **kwargs):
        self.params = kwargs
        self.min_supp = 1
        self.par_num = 1  # TODO: Should this be 2?
                          # FIX:  May (1975) p. 114 says the broken stick has
                          # only 1 parameter, S. 


    @doc_inherit
    def pmf(self, n):
        # TODO:  PMF is not quite summing to one

        # Get parameters
        S, N, n = self.get_params(n=n) 
        assert np.all(S < N), 'n_samp/S must be less than tot_obs/N'
        
        # Calculate pmf
        # TODO: If x and n are the same, change x to n in eq for clarity?
        eq = lambda x, S, N: ((S - 1) / N) * (1 - (x / N)) ** (S - 2)
        pmf = []

        for tS, tN, tn in zip(S, N, n):
            tpmf = eq(tn, tS, tN)
            pmf.append(tpmf)

        return pmf, None


    @doc_inherit
    def rad(self):

        # Get parameters
        S, N = self.get_params()
        assert np.all(S < N), 'n_samp/S must be less than tot_obs/N'

        # Calculate rad
        rad = []
        for tS, tN in zip(S, N):
            trad = np.empty(tS)
            for i in xrange(tS):
                n = np.arange(i + 1, tS + 1) 
                trad[i] = (tN / tS) * sum(1 / n)
            rad.append(trad)

        return rad
            

class sugihara(Distribution):
    __doc__ = Distribution.__doc__ + \
    '''
    Description
    -----------
    Sugihara Rank Abundance Distribution (Sugihara 1980)

    Parameters
    ----------
    S or n_samp : int or iterable
        Total number of species / samples
    N or tot_obs: int or iterable
        Total number of individuals / observations

    Notes
    -----
    This is a sampled rank abundance distribution.  It is not derived
    analytically.

    As S gets bigger, the number of possible breakage sequences for a given
    set of species will increase factorially.  To address this issue, one 
    should increase the number of samples via the parameter sample_size (in
    the rad method) to sample all possible breakage sequences. However, for S
    greater than 20 it is computationally difficult to sample all possible
    breakage sequences and one must realize that the resulting rank
    abundance distribution may not be a perfect representation of sequential
    breaking.  

    The rad method has an additional optional argument for sample_size, which 
    is set to 10000 by default. 

    '''
    # TODO: Explain in docstring that this is sampled, not exact
    # TODO: Back-derive pmf?
    # TODO: Make Error/Warning if pmf or cdf called
    # TODO: Can the docstring be clarified? What is "large"?
    # Fix: Addressed most of these issues in the doc string.
    
    @doc_inherit
    def __init__(self, **kwargs):
        self.params = kwargs
        self.min_supp = 1
        self.par_num = 1
    

    @doc_inherit
    def rad(self, sample_size=10000):
        
        # Get parameters
        S, N = self.get_params()
        assert np.all(S < N), 'n_samp/S must be less than tot_obs/N'

        # Calculate rad
        rad = []
        for tS, tN in zip(S, N):
            total = []
            for i in xrange(sample_size):
                U = np.random.triangular(0.5, 0.75, 1, size=tS - 1)
                p = []
                # TODO: Could this be refactored to perform better?
                for i in xrange(tS):
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
            means = np.array([np.mean(total_array[:,i]) for i in xrange(tS)])
            rad.append(tN * means)

        return rad


class binm(Distribution):
    __doc__ = Distribution.__doc__ + \
    ''' 
    Description
    -----------
    Binomial distribution (ie, random placement model)

    Parameters
    ----------
    tot_obs : int or array-like object
        Total number of individuals in landscape
    n_samp : int or array-like object
        Number of bins/cells sampled.

    Var
    ---
    p : list of floats
        p parameter of is equal to 1 / n_samp, the 'p' parameter of the
        binomial distribution
    tot_obs : list of ints
        tot_obs is the N parameter of the binomial distribution

    Notes
    -----
    The inverse of n_samp is 'a', the ratio of cell size to the total area.
    Also 1 / n_samp can be considered the 'p' parameter of the binomial 
    distribution.

    '''

    @doc_inherit
    def __init__(self, **kwargs):
        self.params = kwargs
        self.min_supp = 0
        self.par_num = 2
    
    @doc_inherit
    def pmf(self, n):
        n_samp, tot_obs, n = self.get_params(n=n)
        
        pmf = []
        var = {}
        var['p'] = []
        var['tot_obs'] = []
        for tn_samp, ttot_obs, tn in zip(n_samp, tot_obs, n):
            ta = 1 / tn_samp
            pmf.append(stats.binom.pmf(tn, ttot_obs, ta))
            var['p'].append(ta)
            var['tot_obs'].append(ttot_obs)
        return pmf, var
    
    @doc_inherit
    def cdf(self, n):
        n_samp, tot_obs, n = self.get_params(n=n)

        cdf = []
        var = {}
        var['p'] = []
        var['tot_obs'] = []
        for tn_samp, ttot_obs, tn in zip(n_samp, tot_obs, n):
            ta = 1 / tn_samp
            cdf.append(stats.binom.cdf(tn, ttot_obs, ta))
            var['p'].append(ta)
            var['tot_obs'].append(ttot_obs)
        return cdf, var

class pois(Distribution):
    __doc__ = Distribution.__doc__ + \
    '''
    Description
    -----------
    Poisson distribution

    Parameters
    ----------
    tot_obs : int or array-like object
        Total number of individuals in landscape
    n_samp : int or array-like object
        Number of bins/cells sampled.

    Var
    ---
    mu : list of floats
        the mu parameter of the poisson distribution

    '''

    @doc_inherit
    def __init__(self, **kwargs):
        self.params = kwargs
        self.min_supp = 0
        self.par_num = 1
    
    @doc_inherit
    def pmf(self, n):

        n_samp, tot_obs, n = self.get_params(n=n)
        
        pmf = []
        var = {}
        var['mu'] = []
        for tn_samp, ttot_obs, tn in zip(n_samp, tot_obs, n):
            tmu = ttot_obs * (1 / tn_samp)
            pmf.append(stats.poisson.pmf(tn, tmu))
            var['mu'].append(tmu)
        return pmf, var 
    
    @doc_inherit
    def cdf(self, n): 
        
        n_samp, tot_obs, n = self.get_params(n=n)

        cdf = []
        var = {}
        var['mu'] = []
        for tn_samp, ttot_obs, tn in zip(n_samp, tot_obs, n):
            tmu = ttot_obs * (1 / tn_samp)
            cdf.append(stats.poisson.cdf(tn, tmu))
            var['mu'].append(tmu)
        return cdf, var

class nbd(Distribution):
    __doc__ = Distribution.__doc__ + \
    '''
    Description
    -----------
    Negative binomial distribution

    Parameters
    ----------
    tot_obs : int or array-like object
        Total number of individuals in landscape
    n_samp : int or array-like object
        Number of bins/cells sampled.
    k : int
        Aggregation parameter

    Var
    ---
    Parameterization differs for different forms of the nbd.  We use the
    standard ecological form as described by Ben Bolker. Parameters 'a' (1 /
    n_samp), 'tot_obs', and k are used to derive the nbd parameter p (see code
    for details).  Parameters k and p are used to generate distribution.
        
    k : list of floats
        k parameter of nbd
    p : list of floats 
        p parameters of nbd

    '''
    
    @doc_inherit
    def __init__(self, **kwargs):
        self.params = kwargs
        self.min_supp = 0
        self.par_num = 2
    
    def pmf(self, n):
        '''
        Probability mass function method.

        Parameters
        ----------
        n : int, float or array-like object
            Values at which to calculate pmf. May be a list of same length as 
            parameters, or single iterable.

        Returns
        -------
        pmf : list of ndarrays
            List of 1D arrays of probability of observing sample n.
        vars : dict containing lists of floats
            Intermediate parameter variables calculated for pmf, as 
            dict.

        See class docstring for more specific information on this distribution.
        '''
        
        n_samp, tot_obs, n = self.get_params(n=n)
        k = make_array(self.params.get('k', None))
        assert k[0] != None, "k parameter not given"
        
        pmf = []
        var = {}
        var['p'] = []
        var['k'] = []

        for tn_samp, ttot_obs, tk, tn in zip(n_samp, tot_obs, k, n):
            tmu = ttot_obs * (1 / tn_samp)
            tp = 1 / (tmu / tk + 1) # See Bolker book Chapt 4
            pmf.append(scipy.stats.nbinom.pmf(tn, tk, tp))
            var['p'].append(tp)
            var['k'].append(tk)
        return pmf, var 

    def cdf(self, n):
        '''
        Cumulative distribution method.  

        Parameters
        ----------
        n : int, float or array-like object
            Values at which to calculate cdf. May be a list of same length as 
            parameters, or single iterable.

        Returns
        -------
        cdf : list of ndarrays
            List of 1D arrays of probability of observing sample n.
        vars : dict containing lists of floats
            Intermediate parameter variables calculated for cdf, as dict.

        See class docstring for more specific information on this distribution.
        '''

        n_samp, tot_obs, n = self.get_params(n=n)
        k = make_array(self.params.get('k', None))
        assert k[0] != None, "k parameter not given"
        
        cdf = []
        var = {}
        var['p'] = []
        var['k'] = []

        for tn_samp, ttot_obs, tk, tn in zip(n_samp, tot_obs, k, n):
            tmu = ttot_obs * (1 / tn_samp)
            tp = 1 / (tmu / tk + 1) # See Bolker book Chapt 4
            cdf.append(scipy.stats.nbinom.cdf(tn, tk, tp))
            var['p'].append(tp)
            var['k'].append(tk)
        return cdf, var 
    
    def fit(self, data, guess_for_k=1):
        '''
        Fit method.

        Uses input data to get best fit parameters for distribution, and stores 
        these parameters in params attribute.
        
        Parameters
        ----------
        data : list of ndarrays
            Data to use to fit parameters of distribution. Even if only one 
            data array, must be in a list with one element.
        guess_for_k : float
            Initial guess for parameter k in solver

        See class docstring for more specific information on this distribution.
        '''

        super(nbd, self).fit(data)
        n_samp, tot_obs = self.get_params()

        data = [np.array(tdata) for tdata in data]
        tempk = []

        for tdata, tn_samp, ttot_obs in zip(data, n_samp, tot_obs): 

            def nll_nb(k):
                self.params['tot_obs'] = ttot_obs
                self.params['n_samp'] = tn_samp
                self.params['k'] = k
                return -sum(np.log(self.pmf(tdata)[0][0]))

            mlek = scipy.optimize.fmin(nll_nb, np.array([guess_for_k]), 
                                                                    disp=0)[0]
            tempk.append(mlek)
        self.params['k'] = tempk
        self.params['n_samp'] = n_samp
        self.params['tot_obs'] = tot_obs
        return self

class nbd_lt(nbd):
    '''
    Description
    -----------
    Zero Truncated Negative Binomial

    Parameters
    ----------
    tot_obs : int or array-like object
        Total number of individuals in landscape
    n_samp : int or array-like object
        Number of bins/cells sampled.
    k : int
        Aggregation parameter

    Var
    ---
    Parameterization differs for different forms of the nbd.  We use the
    standard ecological form as described by Ben Bolker. Parameters 'a' (1 /
    n_samp), 'tot_obs', and k are used to derive the nbd parameter p (see code
    for details).  Parameters k and p are used to generate distribution.
        
    k : list of floats
        k parameter of nbd
    p : list of floats 
        p parameters of nbd

    '''

    def __init__(self, **kwargs):
        self.params = kwargs
        self.min_supp = 1
        self.par_num = 2
    
    @doc_inherit
    def pmf(self, n):
        n_samp, tot_obs, n = self.get_params(n=n)
        k = make_array(self.params.get('k', None))
        assert k[0] != None, "k parameter not given"

        reg_nbd = nbd(n_samp=n_samp, tot_obs=tot_obs, k=k)
        reg_pmf, reg_var = reg_nbd.pmf(n)
        reg_pmf0, reg_var0 = reg_nbd.pmf(0)

        trunc_pmf = [(pr / (1 - p0)) for pr, p0 in zip(reg_pmf, reg_pmf0)]

        return trunc_pmf, reg_var         

class fnbd(Distribution):
    __doc__ = Distribution.__doc__ + \
    '''
    Description
    -----------
    Finite negative binomial (Zillio and He 2010).

    Parameters
    ----------
    tot_obs : int or array-like object
        Total number of individuals in landscape
    n_samp : int or array-like object
        Number of bins/cells sampled.
    k : int
        Aggregation parameter

    Var
    ---
    p : list of floats
        p parameter is 1 / n_samp
    k : list of floats
        k parameter is aggregation parameter

    '''
    
    @doc_inherit
    def __init__(self, **kwargs):
        self.params = kwargs
        self.min_supp = 0
        self.par_num = 2
    
    @doc_inherit
    def pmf(self, n):

        # TODO: Fix to work if n and N are one value
        #    if not (n <= N).all():
        #        raise Exception, "All values of n must be <= N."
        #    elif (a <= 0) or (a >= 1):
        #        raise Exception, "a must be between 0 and 1"

        n_samp, tot_obs, n = self.get_params(n=n)
        k = make_array(self.params.get('k', None))
        assert k[0] != None, "k parameter not given"
        
        pmf = []
        var = {}
        var['p'] = []
        var['k'] = []

        for tn_samp, ttot_obs, tk, tn in zip(n_samp, tot_obs, k, n):

            ln_L = lambda n_i,N,a,k: _ln_choose(n_i+k-1,n_i) + \
                _ln_choose(N-n_i+(k/a)-k-1,N-n_i) - _ln_choose(N +(k/a)-1,N)
            ta = 1 / tn_samp
            tpmf = ln_L(tn, ttot_obs, ta, tk) # Already log
            pmf.append(np.exp(tpmf))
            var['p'].append(ta)
            var['k'].append(tk)
        return pmf, var
    
    def fit(self, data, upper_bnd=5):
        '''
        Fit method.

        Uses input data to get best fit parameters for distribution, and stores 
        these parameters in params attribute.
        
        Parameters
        ----------
        data : list of ndarrays
            Data to use to fit parameters of distribution. Even if only one 
            data array, must be in a list with one element.
        upper_bnd : int
            upper_bnd for parameter k in solver

        See class docstring for more specific information on this distribution.
        '''
        super(fnbd, self).fit(data)
        n_samp, tot_obs = self.get_params()

        data = [np.array(tdata) for tdata in data]
        tempk = []

        for tdata, tn_samp, ttot_obs in zip(data, n_samp, tot_obs): 

            def nll_nb(k):
                self.params['tot_obs'] = ttot_obs
                self.params['n_samp'] = tn_samp
                self.params['k'] = k
                return -sum(np.log(self.pmf(tdata)[0][0]))

            mlek = scipy.optimize.brute(nll_nb, ((1e-10, upper_bnd),))
            tempk.append(mlek)

        self.params['k'] = tempk
        self.params['n_samp'] = n_samp
        self.params['tot_obs'] = tot_obs
        return self

class geo(Distribution):
    __doc__ = Distribution.__doc__ + \
    '''
    Description
    -----------
    Geometric distribution.  Uses nbd object as wrapper with k = 1

    Parameters
    ----------
    tot_obs : int or array-like object
        Total number of individuals in landscape
    n_samp : int or array-like object
        Number of bins/cells sampled.

    Var
    ---
    p :  list of floats
        p parameter is equal to 1 / n_samp

    '''
    
    @doc_inherit
    def __init__(self, **kwargs):
        self.params = kwargs
        self.min_supp = 0
        self.par_num = 1
    
    @doc_inherit
    def pmf(self, n):

        n_samp, tot_obs, n = self.get_params(n=n)
        k = np.repeat(1, len(n_samp))
        pmf, nvar = nbd(tot_obs=tot_obs, n_samp=n_samp, k=k).pmf(n)
        var = {}
        var['p'] = 1 / n_samp
        return pmf, var
    
    @doc_inherit
    def cdf(self, n):

        n_samp, tot_obs, n = self.get_params(n=n)
        k = np.repeat(1, len(n_samp))
        cdf, nvar = nbd(tot_obs=tot_obs, n_samp=n_samp, k=k).cdf(n)
        var = {}
        var['p'] = 1 / n_samp
        return cdf, var
        
class fgeo(Distribution):
    __doc__ = Distribution.__doc__ + \
    '''
    Description
    -----------
    Finite geometric pmf (Zillio and He 2010). Use fnbd object as wrapper with 
    k = 1

    Parameters
    ----------
    tot_obs : int or array-like object
        Total number of individuals in landscape
    n_samp : int or array-like object
        Number of bins/cells sampled.

    Var
    ---
    p : list of floats
        a parameter is equal to 1 / n_samp

    Notes
    -----
    This is not a true geometric distribution - calculate a pmf z and run 
    z[1:]/z[0:-1], noting that the ratio is not constant.

    '''
    
    @doc_inherit
    def __init__(self, **kwargs):
        self.params = kwargs
        self.min_supp = 0
        self.par_num = 1
    
    @doc_inherit
    def pmf(self, n):
        
        n_samp, tot_obs, n = self.get_params(n=n)
        k = np.repeat(1, len(n_samp))
        pmf, nvar = fnbd(tot_obs=tot_obs, n_samp=n_samp, k=k).pmf(n)
        var = {}
        var['p'] = nvar['p']
        return pmf, var
    
    @doc_inherit
    def cdf(self, n):
        
        n_samp, tot_obs, n = self.get_params(n=n)
        k = np.repeat(1, len(n_samp))
        cdf, nvar = fnbd(tot_obs=tot_obs, n_samp=n_samp, k=k).cdf(n)
        var = {}
        var['p'] = nvar['p']
        return cdf, var

class tgeo(Distribution):
    __doc__ = Distribution.__doc__ + \
    '''
    Description
    -----------
    Truncated geometric distribution (Harte 2011)

    Parameters
    ----------
    tot_obs : int or array-like object
        Total number of individuals in landscape
    n_samp : int or array-like object
        Number of bins/cells sampled.

    Var
    ---
    x : list of floats
        -np.log(x) is the lagrange multiplier for the tgeo distribution

    Notes
    -----
    This is a true geometric distribution in which p = exp(-lambda). 
    Confirmed with z[1:]/z[0:-1].

    'a' is equal to 1 / n_samp
`
    Plotting eq for various realistic values of p, N, and a confirms that 
    this is a smooth, well-behaved function that should be amenable to 
    using a root finder. When a is large and N is small (i.e. a = .9 and
    N = 32) the lagrange multiplier has no solution.  For large a's (.8 or
    greater), N needs to be sufficiently large for lambda pi to have a 
    solution.

    ''' 
    
    @doc_inherit
    def __init__(self, **kwargs):
        self.params = kwargs
        self.min_supp = 0
        self.par_num = 1
    
    @doc_inherit
    def pmf(self, n):

        n_samp, tot_obs, n = self.get_params(n=n)
        #NOTE: Overflow warning but not affecting results
        eq = lambda x, N, a: ((x / (1 - x)) - (((N + 1) * x ** (N + 1)) / \
                            (1 - x ** (N + 1)))) - (N * a)
        pmf = []
        var = {}
        var['x'] = []
        for tn_samp, ttot_obs, tn in zip(n_samp, tot_obs, n):
            ta = 1 / tn_samp
            x = scipy.optimize.brentq(eq, 0, min((sys.float_info[0] *
                ta)**(1/float(ttot_obs)), 2), args=(ttot_obs, ta), disp=False)
            z = (1 - x ** (ttot_obs + 1)) / (1 - x)
            tpmf = (1 / z) * (x ** tn)
            pmf.append(tpmf)
            var['x'].append(x)

        return pmf, var

class mete_sar_iter(Curve):
    __doc__ = Distribution.__doc__ + \
    '''
    Description
    -----------
    This class explores the METE generated SAR

    Parameters
    ----------
    S : int
        Total number of species at the anchor area
    N : int
        Total number of individuals at the anchor area

    Notes
    -----
    This class uses method 1 in Harte (2011) to calculate the SAR
    
    '''
    
    @doc_inherit
    def __init__(self, **kwargs):
        self.params = kwargs

    def vals(self, a_list, upscale=0, downscale=0):
        '''
        Predict the universal SAR curve for the given S and N found at 
        the given anchor scale

        Parameters
        ----------
        a_list : array-like or None
            Target areas for which to calculate SAR
        upscale : int
            Number of iterations up from the anchor scale.  Each iteration 
            doubles the previous area. Only active if a_list is None.
        downscale : int
            Number of iterations down from the anchor scale. Each iteration 
            halves the previous area. Only active if a_list is None.

        Returns
        -------
        : 1D structured np.array
            The structured array has fields 'items' and 'area'

        Notes
        -----
        With this method of the METE SAR, one cannot calculate the SAR at exact
        areas.  Rather this method iterates up and down by powers of 2.
        Therefore, the output of this function will contain all the SAR
        calcutions in between ~min(a_list) ~max(a_list).


        '''
        #Get and check params
        S = self.params.get('S', None)
        N = self.params.get('N', None)
        assert S != None, "S parameter not given"
        assert N != None, "N parameter not given"
        assert S < N, "S must be less than N"
        assert S > 1, "S must be greater than 1"
        assert N > 0, "S must be greater than 0"
        if not np.iterable(a_list) and a_list is not None:
            raise TypeError('a_list is not an array-like object')
        
        anch = 1
        if a_list != None:
            mint = np.min(a_list)
            maxt = np.max(a_list)
            if mint == anch and maxt == anch: 
                upscale = 0; downscale = 0
            elif (mint > anch or mint == anch) and maxt > anch:
                upscale = np.int(np.ceil(np.log2(maxt / anch)))
                downscale = 0
            elif mint < anch and (maxt < anch or maxt == anch):
                downscale = np.int(np.ceil(np.abs(np.log2(mint / anch)))) 
                upscale = 0
            elif mint < anch and maxt > anch:
                upscale = np.int(np.ceil(np.log2(maxt / anch)))
                downscale = np.int(np.ceil(np.abs(np.log2(mint / anch)))) 
    
        if upscale == 0 and downscale == 0:
            return np.array((S, anch), dtype=[('items', np.float),
                                                ('area', np.float)])
        areas = _generate_areas_(anch, upscale, downscale)
        sar = np.empty(len(areas), dtype=[('items', np.float),
                                      ('area', np.float)])
        sar['area'] = areas
        if upscale != 0: 
            sar['items'][downscale:] = _upscale_sar_(areas[downscale:], N, S)
        if downscale != 0:
            sar['items'][:downscale + 1] =\
                                   _downscale_sar_(areas[:downscale + 1], N, S)
        return sar
   

    def univ_curve(self, num_iter=10):
        '''
        This function calculates the universal SAR curve.

        Parameters
        ----------
        num_iter : int
            Number of iterations.
            WARNING: Running more than 10 begins to take along time 
    
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
        function were taken from Harte et al. (2009) and Harte (2011). This 
        uses method 1 in Harte (2011)

        '''
        
        S = self.params.get('S', None)
        N = self.params.get('N', None)
        assert S != None, "S parameter not given"
        assert N != None, "N parameter not given"
        assert S < N, "S must be less than N"
        assert S > 1, "S must be greater than 1"
        assert N > 0, "S must be greater than 0"
        num_ind = np.empty(num_iter + 1); spp = np.empty(num_iter + 1)
        univ_SAR = np.empty(num_iter + 1, dtype=[('z', np.float),
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

class powerlaw(Curve):
    __doc__ = Distribution.__doc__ + \
    '''
    Description
    -----------
    A power law curve to describe an SAR

    Parameters
    ---------- 
    c : int
        intercept of the loglog power law.  If passing in area fractions, c can
        be considered S at the anchor scale.
    z : int
        Slope of the loglog power law
    
    '''

    @doc_inherit
    def __init__(self, **kwargs):
        self.params = kwargs
    
    @doc_inherit
    def vals(self, a_list):
        '''
        Generate a power law SAR with a given z and c.

        Parameters
        ----------
       a_list : array-like object
            List of area fractions.
    
        Returns
        -------
        : structured np.array
            A structured np.array with dtype=[('items', np.float),
            ('area', np.float)]. 
    
        '''
        z = self.params.get('z', None)
        c = self.params.get('c', None)
        assert z != None, "z parameter not given"
        assert c != None, "c parameter not given"
        a_list = make_array(a_list)
        output_array = np.empty(len(a_list), dtype=[('items', np.float),
                                                     ('area', np.float)])
        output_array['area'] = a_list
        p_law = lambda x: c * (x ** z)
        output_array['items'] = p_law(a_list)
        return output_array
    
    @doc_inherit
    def fit(self, data, full_sad):

        super(powerlaw, self).fit(data, full_sad)
        reg = stats.linregress(np.log(data[0]), np.log(data[1]))
        self.params['z'] = reg[0]
        self.params['c'] = np.exp(reg[1])
        return self

class gen_sar(Curve):
    __doc__ = Distribution.__doc__ + \
    '''
    A generic sar function the utilizes the relationship between the sad
    and the ssad to generate the sar.  Can take any combination of sad and
    ssad. 

    Parameters
    ----------
    S : float
        Total number of species at the anchor area
    N : float
        Total number of individuals at the anchor area

    Notes
    -----
    In order to use gen_sar method vals(), self.params['sad_pmf'] must exist.
    Running gen_sar method fit() fills self.params['sad_pmf'] or it can be
    filled manually.  Generic sar can only downscale.  

    '''
    #NOTE: Might be making this too limiting be forcing an sad object in.
    #Just going to do it for now
    def __init__(self, sad, ssad, **kwargs):
        '''

        Parameters
        ----------
        sad : a sad distribution object
            A distribution object with minimum support equal to 1. pmf of sad
            from 1 to N should sum to approximately 1.
        ssad : a ssad distribution object
            A distribution object with minimum support equal to 0.

        Notes
        -----
        Generic sar must take in a sad and ssad distribution object upon
        instantiation

        '''
        self.sad = sad
        self.ssad = ssad
        self.params = kwargs

    def vals(self, a_list):
        '''

        Parameters
        ----------
        a_list : array-like object
            List of area fractions at which to calculate the SAR

        Returns
        -------
        : structured np.array
            A structured np.array with dtype=[('items', np.float),
            ('area', np.float)]. Items can be species.

        '''
        sad = self.params.get('sad_pmf', None)
        assert sad != None, "params['sad_pmf'] does not exist.  Try fitting" \
                    + " gen_sar object or initialize self.params['sad_pmf']"
        ssad = self.ssad
        S = self.params.get('S', None)
        assert S != None, "S parameter not given"
        sar = []
        N_range = np.arange(1, len(sad) + 1)
        ssad.params['tot_obs'] = N_range
        for i, a in enumerate(a_list):
            if not a <= 1:
                raise ValueError("a must be less than or equal to 1")
            if a == 1:
                p_pres_list = np.repeat(1, len(N_range))
            else:
                ssad.params['n_samp'] = np.repeat(1 / a, len(N_range))
                #import pdb; pdb.set_trace()
                p_pres_list = [1 - absnt[0] for absnt in ssad.pmf(0)[0]]
            sar.append(sum(S * sad * np.array(p_pres_list)))
        return np.array(zip(sar, a_list), dtype=[('items', np.float), 
                                                  ('area', np.float)])
    
    def fit(self, data, full_sad):
        '''
        This fit method fills the required parameters for an SARCurve object.
        For the gen_sar object, the fit method will fill self.params['sad_pmf']
        which is required by the vals function. 

        Parameters
        ----------
        data :  tuple of array-like objects
            data contains two array-like object.  The first is a list of area
            fractions and the second is a list of species/items numbers 
            corresponding to those areas.
        full_sad : array-like object
            The full_sad at the anchor area (species/items counts)

        '''
        super(gen_sar, self).fit(data, full_sad)
        self.sad.fit([full_sad])
        self.params['sad_pmf'] = self.sad.pmf(np.arange(1, self.params['N'] + 
                                                                    1))[0][0]
        return self

###########################
#---Generic SAR classes---#
###########################

class logser_ut_tgeo(gen_sar):
    '''
    Generic sar: 
    SAD : Truncated logseries (logser_ut)
    SSAD: Truncated geometric (tgeo)

    See class gen_sar for more information
    '''
    def __init__(self, **kwargs):
        self.params = kwargs
        self.sad = logser_ut()
        self.ssad = tgeo()

class logser_ut_fgeo(gen_sar):
    '''
    Generic sar: 
    SAD : Truncated logseries (logser)
    SSAD: Finite geometric (fgeo)

    See class gen_sar for more information
    '''
    def __init__(self, **kwargs):
        self.params = kwargs
        self.sad = logser_ut()
        self.ssad = fgeo()

class logser_ut_binm(gen_sar):
    '''
    Generic sar: 
    SAD : Truncated logseries (logser)
    SSAD: Binomial (binm)

    See class gen_sar for more information
    '''
    def __init__(self, **kwargs):
        self.params = kwargs
        self.sad = logser_ut()
        self.ssad = binm()

class plognorm_lt_binm(gen_sar):
    '''
    Generic sar: 
    SAD : Truncated Poisson lognormal (plognorm_lt)
    SSAD: Binomial (binm)

    See class gen_sar for more information
    '''
    def __init__(self, **kwargs):
        self.params = kwargs
        self.sad = plognorm_lt()
        self.ssad = binm()

class plognorm_lt_tgeo(gen_sar):
    '''
    Generic sar: 
    SAD : Truncated Poisson lognormal (plognorm_lt)
    SSAD: Truncated geometric (tgeo)

    See class gen_sar for more information
    '''
    def __init__(self, **kwargs):
        self.params = kwargs
        self.sad = plognorm_lt()
        self.ssad = tgeo()

class plognorm_lt_fgeo(gen_sar):
    '''
    Generic sar: 
    SAD : Truncated Poisson lognormal (plognorm_lt)
    SSAD: Finite geometric (fgeo)

    See class gen_sar for more information
    '''
    def __init__(self, **kwargs):
        self.params = kwargs
        self.sad = plognorm_lt()
        self.ssad = fgeo()

#########################
## Energy Distributions##
#########################

#TODO: Another base class for continuous distributions?
#Should not inherit Distribution
class psi(Distribution):
    '''
    The community energy distribution (psi) described in Harte (2011)

    Parameters
    ----------
    S or n_samp : int or iterable
        Total number of species / samples
    N or tot_obs: int or iterable
        Total number of individuals / observations
    E : int or iterable
        Total energy output of community

    Var
    ---
    beta : list of floats
        The beta lagrange multiplier
    l2 : list of floats
        The lambda2 lagrange multiplier

    Notes
    -----
    All other lagrange multipliers can be calculated from beta and l2.

    '''


    def __init__(self, **kwargs):
        '''
        Initialize distribution object.

        Stores keyword arguments as params attribute and defines minimum 
        support for distribution.

        Parameters
        ----------
        kwargs : comma separate list of keyword arguments
            Parameters needed for distribution. Stored in dict as params 
            attribute.

        See class docstring for more specific information on this distribution.
        '''
    
        self.params = kwargs
        self.min_supp = 1


    def pdf(self, e):
        '''
        The probability density function for the community energy distribution.
        This distribution is described by Harte (2011).

        Parameters
        ----------
        e : int or array-like object
            The values at which to calculate the pdf

        Returns
        -------
        : tuple
            A tuple of two objects.  The fist object is a dict that has
            keywords 'beta' and 'l2'.  These keywords look up arrays that have
            the same length as self.params['N'].  The second object is a list 
            of arrays with the same length as the parameters in self.params.
            The members of the arrays are the pdf calculates as the given
            values of e. 

        '''
        #Get and check parameters
        S, N, e = self.get_params(n=e)
        E = make_array(self.params.get('E', None))
        assert E[0] != None, 'E parameter not given'
        assert np.all(E > N), 'E must be greater than N'

        start = 0.3
        stop = 2
        flmax = sys.float_info[0]
        eq = lambda x,k,N,S: sum(x ** k / float(N) * S) -  sum((x ** k) / k)

        pdf = []
        var = {}
        var['beta'] = []
        var['l2'] = []

        for tS, tN, tE, te in zip(S, N, E, e):
            k = np.linspace(1, N, num=N)
            tx = scipy.optimize.brentq(eq, start,
                                       min((flmax/tS)**(1/float(tN)), stop), 
                                       args = (k, tN, tS), disp=True)
            tbeta = -np.log(tx)
            tl2 = float(tS) / (tE - tN) # Harte (2011) 7.26
            tl1 = tbeta - tl2
            tsigma = tl1 + (tE * tl2)

            norm = (float(tS) / (tl2 * tN)) * (((np.exp(-tbeta) - \
                    np.exp(-tbeta*(tN + 1))) / (1 - np.exp(-tbeta))) - \
                    ((np.exp(-tsigma) - np.exp(-tsigma*(tN + 1))) / \
                    (1 - np.exp(-tsigma)))) #Harte (2011) 7.22

            #Notation from E.W.
            exp_neg_gamma = np.exp(-(tbeta + (te - 1) * tl2))

            tpdf = (float(tS) / (tN * norm)) * ((exp_neg_gamma  / \
                    (1 - exp_neg_gamma)**2) - ((exp_neg_gamma**tN / \
                    (1 - exp_neg_gamma)) * (tN + (exp_neg_gamma / \
                    (1 - exp_neg_gamma)))))
                    # Harte (2011) 7.24
            pdf.append(tpdf)
            var['beta'].append(tbeta)
            var['l2'].append(tl2)
        
        return pdf, var
    
    def red(self):
        '''
        Rank energy distribution derived from the psi distribution.  Predicts
        energy values for entire community (N individuals). 

        Parameters
        ----------
        None

        Returns
        -------
        : list
            A list of arrays containing the rank energy distribution for a
            given set of S, N, and E.  The length of the list is the same
            length as self.params['N']. 

        '''

        S, N = self.get_params()
        E = make_array(self.params.get('E', None))
        assert E[0] != None, 'E parameter not given'
        assert np.all(E > N), 'E must be greater than N'

        start = 0.3
        stop = 2
        flmax = sys.float_info[0]
        eq = lambda x,k,N,S: sum(x ** k / float(N) * S) -  sum((x ** k) / k)

        n_arrays = [np.arange(1, i + 1) for i in N]
        prad = lambda beta, r, N, l1, l2: (1 / l2) * np.log(((beta * N) + \
                                  r - 0.5) / (r - 0.5)) - (l1 / l2)

        rad = []
        for tS, tN, tE, tn, in zip(S, N, E, n_arrays):

            k = np.linspace(1, N, num=N)
            tx = scipy.optimize.brentq(eq, start,
                                       min((flmax/tS)**(1/float(tN)), stop), 
                                       args = (k, tN, tS), disp=True)
            tbeta = -np.log(tx)
            tl2 = float(tS) / (tE - tN) # Harte (2011) 7.26
            tl1 = tbeta - tl2

            trad = prad(tbeta, tn, tN, tl1, tl2)
            rad.append(trad)

        return rad

#TODO: Theta should inherit a different base class
class theta(Distribution):
    '''
    The species specific energy distribution (theta) as described by Harte
    (2011).

    Parameters
    ----------
    S or n_samp : int or iterable
        Total number of species / samples
    N or tot_obs: int or iterable
        Total number of individuals / observations
    E : int or iterable
        Total energy output of community
    n : int or iterable
        Number of individuals in a given species

    '''

    def __init__(self, **kwargs):
        '''
        Initialize distribution object.

        Stores keyword arguments as params attribute and defines minimum 
        support for distribution.

        Parameters
        ----------
        kwargs : comma separate list of keyword arguments
            Parameters needed for distribution. Stored in dict as params 
            attribute.

        See class docstring for more specific information on this distribution.
        '''

        self.params = kwargs
        self.min_supp = 1

    def pdf(self, e):
        '''
        Probaility density function for the species specific energy
        distribution.

        Parameters
        ----------
        e : int or array-like object
            The values at which to calculate the pdf

        Returns
        -------
        : tuple
            A tuple of two objects.  The fist object is a dict that has
            keyword 'l2'.  This keyword looks up an array that has
            the same length as self.params['N'].  The second object is a list 
            of arrays with the same length as the parameters in self.params.
            The members of the arrays are the pdf calculates as the given
            values of e. 
            
        '''
        S, N, e = self.get_params(n=e)
        E = make_array(self.params.get('E', None))
        n = make_array(self.params.get('n', None))
        assert E[0] != None, 'E parameter not given'
        assert np.all(E > N), 'E must be greater than N'
        assert np.all(n <= N), 'n must be less than or equal to N'

        pdf = []
        var = {}
        var['l2'] = []

        for tS, tN, tE, tn, te in zip(S, N, E, n, e):

            tl2 = float(tS) / (tE - tN)
            tpdf = (tn * tl2 * np.exp(-tl2 * tn * te)) / (np.exp(-tl2 * tn)\
                                - np.exp(-tl2 * tn * tE)) #Harte (2011) 7.25

            pdf.append(tpdf)
            var['l2'].append(tl2)
        
        return pdf, var

    def red(self):
        '''
        Rank energy distribution for species specific energy distribution.
        Predicts energy values for every individual in a given species with
        abundance n. 

        Parameters
        ----------
        None

        Returns
        -------
        : list
            A list of arrays containing the rank energy distribution for a
            given set of S, N, E, n.  The length of the list is the same
            length as self.params['N']. 

        '''

        S, N = self.get_params()
        E = make_array(self.params.get('E', None))
        n = make_array(self.params.get('n', None))
        assert E[0] != None, 'E parameter not given'
        assert np.all(E > N), 'E must be greater than N'
        assert np.all(n <= N), 'n must be less than or equal to N'
        
        n_arrays = [np.arange(1, i + 1) for i in n]

        pred = lambda r, n, l2 : 1 + (1 / (l2 * n)) * np.log( n / (r - 0.5))
        red = []
        for tS, tN, tE, tn, tn_arr in zip(S, N, E, n, n_arrays):
           
           tl2 = float(tS) / (tE - tN)
           tred = pred(tn_arr, tn, tl2)
           red.append(tred)

        return red






        


    



            






        










        



        
        
        



def make_array(n):
    '''Cast n as iterable array.'''
    if np.iterable(n):
        return np.array(n)
    else:
        return np.array([n])


def expand_n(n, size):
    '''Check dimensions of n and expand to match size if necessary.'''
    if np.iterable(n) and np.iterable(n[0]):  # If n is iterable of iterables
        if len(n) != size:
            raise TypeError('If n is a list of iterables, list must be ' +
                                'the same length as parameter lists')
        else:
            new_n = [np.array(tn) for tn in n]
    
    elif np.iterable(n):  # If n is iterable of non-iterables
        new_n = [np.array(np.copy(n)) for i in xrange(size)]

    else:  # If n is not iterable
        new_n = [np.array([n]) for i in xrange(size)]

    return new_n


def make_rank_abund(pmf, n_samp, min_supp=1):
    '''
    Convert any pmf into a rank abundance curve for S species using cumulative 
    distribution function.
 
    Parameters
    ----------
    pmf : ndarray
        Probability of observing a species from 1 to length pmf individs.
    n_samp : int
        Total number of samples 

    Returns
    -------
    S_abunds : ndarray
        1D array of predicted abundance for each species

    Notes
    -----
    Function actually implements (philosophically) a step quantile function.

    '''

    points = np.arange(1/(2*n_samp), 1, 1/n_samp)
    counts = np.zeros(n_samp)
    
    if min_supp == 1:
        pmf = np.array([0] + list(pmf)) # Add 0 to start of pmf
    cum_pmf = np.cumsum(pmf)
    
    for cutoff in cum_pmf:
        greater_thans = (points >= cutoff)
        counts[greater_thans] += 1

        if not greater_thans.any():  # If no greater thans, done with samples
            break
    
    return counts


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
