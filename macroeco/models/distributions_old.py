#!/usr/bin/python

'''
Macroecological distributions and curves.

Distributions
-------------

SAD
- `most_even` -- Most even distribution of balls into bins
- `most_uneven` -- Most uneven distribution of balls into bins
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
- `gen_sar` - A generic sar that supports any sad and ssad distribution. 

SSAD
- `binm` - Binomial distribution (Random Placement Model)
- `pois` - Poisson distribution
- `nbd` - Negative binomial distribution
- `fnbd` - Finite negative binomial (Zillio and He 2010)
- `geo` - Geometric distribution
- `fgeo` - Finite geometric distribution (Zillio and He 2010)
- `tgeo` - Truncated geometric distribution (Harte et al. 2008)

Energy Distributions
- `psi` -- The individual energy distribution described in Harte (2011)
- `theta` -- The species energy distribution (SED) as described by Harte (2011).
- `nu` -- The average species energy distribution (ASED) as described by Harte
  (2011)

Misc Functions
--------------
- `make_array` 
- `make_rank_abund` 
- `_ln_choose`
- `_downscale_sar_`
- `_upscale_sar_`
- `_generate_areas_`
- `expand_n`
- `check_list_of_iterables`
- `set_up_and_down`
- `unpack`

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
import scipy.special as spec
from copy import deepcopy
import math as m
import scipy.integrate as integrate
import sys
#from docinherit import DocInherit
from utils.docinherit import DocInherit

doc_inherit = DocInherit


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

    def get_name(self):
        '''
        Returns the class name
        '''
        return self.__class__.__name__

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

        See class docstring for more specific information on this distribution.
        '''
        # This method does nothing, but exists so that derived class pmf
        # methods can inherit this docstring.
        pass

    def univ_curve(self, num_iter=5, direction='down', param='tot_obs',
                                      iterative=False, base=2, use_rad=False):
        '''
        Generating a univsersal curve for different curves.  A universal curve
        is defined as the slope value (z) at a function y = f(x) at a given x
        plotted against (x / y). This function halfs and doubles x and
        calculates the slope using the average. x = 1 is considered the
        anchor scale.

        Parameters
        ----------
        num_iter : int
            Number of halvings at which to calculate z.
        direction : string
            Either 'down' or 'up'.  The direction to iterate the curve.
        param : string
            If not None, will look for the given string in the self.params
            dict and set multiplier to this value. If None, multiplier is one.
        iterative : bool
            If False, uses the one-shot method to calculate z via Curve.vals.
            If iterative is True, uses the iterative method to calculate z via
            Curve.iter_vals.
        base : int
            Specifies whether you would like each iteration to double/half
            (base=2), triple/third (base=3), etc
        use_rad : bool
            If False, uses the sad pmf to calculate the SAR.  If True, uses the
            sad rank abundance distribution to calculate the SAR. An SAR class
            inherits curve.

        Returns
        -------
        : tuple
            the tuple contains two structured arrays.  The first structured
            array is an array with dtype = [('z', np.float), ('x_over_y',
            np.float)] that contains the slope and the quantity x * multiplier
            / y. This array is the universal curve.  The second array has dtype
            = [('items', np.float), ('area', np.float)] where 'items' is often
            equivalent to species and 'area' is the area fraction.  This second
            array is the species/item area relationship that gives the
            universal curve specified in the first array in the tuple.

        '''
        
        # Allow multiplier to be set to examine different relationships
        if param != None:
            multiplier = self.params.get(param, None)
            assert multiplier != None, "%s not found in self.params" % param
        else:
            multiplier = 1

        
        if iterative:
            def z(a):
                na = [((1. / base) * a[0])] + list(a) + [(base * a[-1])]
                complete_a = self.iter_vals(na, non_iter=True, base=base,
                                                      use_rad=use_rad)['items']
                sz = len(complete_a)
                a1 = complete_a[np.arange(1,sz - 1)]
                a2 = complete_a[np.arange(sz - 2)]
                a3 = complete_a[np.arange(2, sz)]
                return \
                 (0.5 * (np.log(a3/ a2))) / np.log(base), a1, zip(complete_a,
                 na) #Added last return for UNI sar analysis

        else:
            def z(a):
                a1 = self.vals(a, use_rad=use_rad)['items']

                a2 = np.array(list(self.vals(((1./base) * a)[0],
                        use_rad=use_rad)['items']) + list(a1[:-1]))
                a3 = np.array(list(a1[1:]) +  list(self.vals((base *
                                            a)[-1], use_rad=use_rad)['items']))
                return (0.5 * (np.log(a3 / a2))) / np.log(base), a1, zip(a1, a)
        
        # Get the area list
        if direction == 'down':
            a_list = [1 / (base**(i)) for i in np.arange(num_iter + 1)]

        elif direction == 'up':
            a_list = [base**(i) for i in np.arange(num_iter + 1)]
        else:
            raise ValueError('%s not a recognized direction' % direction)

        # Compute universal curve parameters
        a_list.sort()
        a_list = np.array(a_list)
        zs, base_a, spp_area = z(a_list)
        x_over_y = (a_list * multiplier) / base_a 

        uni =  np.array(zip(zs, x_over_y), dtype=
                                      [('z', np.float),('x_over_y', np.float)])
        uni.sort(order=['x_over_y'])

        return uni, np.array(spp_area, dtype=[('items', np.float),
                                                  ('area', np.float)])

    def get_params(self, parameter_list):
        '''
        Gets and validates basic distribution parameters

        Parameters
        ----------
        parameter_list : list
            A list of strings where each string is a keyword for the parameter
            in the self.params dictionary.  

        Returns
        -------
        : tuple
            Validated parameters with the sample length as parameter_list
 
        '''

        retrieved_params = []

        # Get params or None from self.params
        for i, param in enumerate(parameter_list):
            retrieved_params.append(self.params.get(param, None))

            # If parameter not found, raise error
            if retrieved_params[i] is None:
                raise TypeError('%s not found in self.params' % param)
        
        return tuple(retrieved_params)

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
    var : dict
        A dictionary of useful variables that are computed internally to
        generate pmf, pdf, cdf, or rad. If parameters are calculated in the fit
        method, they are included in var dict

    Methods
    -------
    pdf(n)
        Probability density function
    pmf(n)
        Probability mass function
    cdf(n)
        Cumulative distribution function
    rad()
        Rank abundance distribution, calculated from cdf
    fit(data)
        Uses data to populate params attribute

    Examples
    --------

    # logser is a Distribution object
    dist_object = logser(tot_obs=500, n_samp=30)
    pmf_vals = dist_object.pmf(3)
    dist_object = dist_object.pmf([1,1,1,2,3,4])
    fitted = logser().fit([1,1,1,1,2,3,4,5])
    rad = fitted.rad()

    # nbd is a Distribution object
    dist_object = nbd(tot_obs=120, n_samp=20, k=4)
    cdf_vals = dist_object.cdf([1,2,3,4])
    dist_object = nbd(tot_obs=[345,678], n_samp=[23, 45], k=1)
    pmf_vals = dist_object.pmf([[2,34,2], [15, 16, 17]])

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
        # Initializes self.var and self.params
        self.params = kwargs
        self.var = {}



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

        See class docstring for more specific information on this distribution.
        '''
        # This method raises an error if a derived class does not have its own
        # pmf method.
        raise NotImplementedError('PMF is not implemented for this' + 
                                  ' Distribution class')

    def pdf(self, n):
        '''
        Probability density function method.

        Parameters
        ----------
        n : int, float or array-like object
            Values at which to calculate pdf. May be a list of same length as 
            parameters, or single iterable.

        Returns
        -------
        pdf : list of ndarrays
            List of 1D arrays of probability of observing sample n.

        See class docstring for more specific information on this distribution.
        '''
        # This method raises an error if a derived class does not have its own
        # pdf method

        raise NotImplementedError('PDF is not implemented for this' + 
                                  ' Distribution class')


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
            List of 1D arrays of cumulative probability of observing sample n.

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

        # Extend for pdf or pmf
        try:
            pmf_list = self.pdf(n_in)

        except(NotImplementedError):
            pmf_list = self.pmf(n_in)

        # Calculate cdfs
        cdf = []
        for tpmf, tn in zip(pmf_list, n):
            full_cdf = np.cumsum(tpmf)
            tcdf = np.array([full_cdf[x - self.min_supp] for x in tn])
            cdf.append(tcdf)

        return cdf 


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
        n_samp, tot_obs = self.get_params(['n_samp', 'tot_obs'])

        # TODO: Add error or warning if tot_obs is large enough that python is slow
        # FIX: Throw warning that tot_obs is large and that calculation will take
        # awhile.  tot_obs > 100,000. Base CDF should thrown the same warning. 

        # Calculate pmfs, going up to tot_obs for upper limit
        n_arrays = [np.arange(self.min_supp, 1*(i + 1)) for i in tot_obs]
        pmf = self.pmf(n_arrays)
        
        # Calculate rad
        rad = []
        for tn_samp, ttot_obs, tpmf in zip(n_samp, tot_obs, pmf):
            trad = make_rank_abund(tpmf, tn_samp, min_supp=self.min_supp)
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

        data = check_list_of_iterables(data) 
        
        # Check if distribution can support the fitted data
        num_zeros = np.array([sum(dt == 0) for dt in data])
        if np.any(num_zeros != 0) and self.min_supp == 1:
            raise ValueError('%s does not support data with zeros' %
                                                    self.__class__.__name__)

        n_samp = []
        tot_obs = []
        
        for tdata in data:
            n_samp.append(len(tdata))
            tot_obs.append(np.sum(tdata))

        self.params['n_samp'] = n_samp
        self.params['tot_obs'] = tot_obs

        return self


    def get_params(self, parameter_list):
        '''
        Gets and validates basic distribution parameters

        Parameters
        ----------
        parameter_list : list
            A list of strings where each string is a keyword for the parameter
            in the self.params dictionary.  

        Returns
        -------
        : tuple
            Validated parameters with the sample length as parameter_list
 
        '''

        retrieved_params = []

        # Get params or None from self.params
        for i, param in enumerate(parameter_list):
            retrieved_params.append(make_array(self.params.get(param, None)))

            # If parameter not found, raise error
            if retrieved_params[i][0] is None:
                raise TypeError('%s not found in self.params' % param)

        # Check that all params are the same length. If they don't, extend
        # parameters with only one item
        len_ind = [len(p) for p in retrieved_params]
        unq = np.unique(len_ind)
        if len(unq) != 1:
            if len(unq) == 2 and (unq[0] == 1 or unq[1] == 1):
                max_len = np.max(len_ind)
                ones = np.where(np.array(len_ind) == 1)[0]
                for i in ones:
                    retrieved_params[i] = np.repeat(retrieved_params[i],
                                                                       max_len)
            else:
                raise ValueError('Parameters do not have the same length')

        return tuple(retrieved_params)

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

class most_even(Distribution):
    '''
    Description
    -----------
    The most even species abundance distribution where every species get N/S
    individuals

    Parameters
    ----------
    n_samp : int or iterable
        Total number of species / samples
    tot_obs: int or iterable
        Total number of individuals / observations

    self.var keywords
    -----------------
    None

    '''

    @doc_inherit
    def __init__(self, **kwargs):
        
        self.params = kwargs
        self.min_supp = 1
        self.par_num = 2
        self.var = {}

    @doc_inherit
    def pmf(self, n):
        
        # Get parameters
        n_samp, tot_obs = self.get_params(['n_samp', 'tot_obs'])
        n = expand_n(n, len(n_samp))
        
        # TODO: Additional checks?
        assert np.all(n_samp <= tot_obs), 'n_samp must be <= tot_obs'
        
        pmf = []

        for tn_samp, ttot_obs, tn in zip(n_samp, tot_obs, n):
            
            # The mean needs to be a discrete value
            tmean = np.round(ttot_obs / tn_samp, decimals=0)
            tpmf = np.array([0 if k != tmean else 1 for k in tn])
            pmf.append(tpmf)

        return pmf
            
    @doc_inherit
    def rad(self):
        # Get parameters
        n_samp, tot_obs = self.get_params(['n_samp', 'tot_obs' ])

        # Calculate rad
        rad = []
        for tn_samp, ttot_obs in zip(n_samp, tot_obs):
            trad = np.repeat(ttot_obs / tn_samp, tn_samp)
            rad.append(trad)

        return rad

class most_uneven(Distribution):
    '''
    Description
    -----------
    The most uneven distribution possible. All species have one individual
    except one which has N - S + 1 individuals

    Parameters
    ----------
    n_samp : int or iterable
        Total number of species / samples
    tot_obs: int or iterable
        Total number of individuals / observations

    self.var keywords
    -----------------
    None

    '''

    @doc_inherit
    def __init__(self, **kwargs):
        
        self.params = kwargs
        self.min_supp = 1
        self.par_num = 2
        self.var = {}

    @doc_inherit
    def pmf(self, n):
        
        # Get parameters
        n_samp, tot_obs = self.get_params(['n_samp', 'tot_obs'])
        n = expand_n(n, len(n_samp))
        
        # TODO: Additional checks?
        assert np.all(n_samp <= tot_obs), 'n_samp must be <= tot_obs'
        
        pmf = []

        for tn_samp, ttot_obs, tn in zip(n_samp, tot_obs, n):

            tpmf = []
            tmax = ttot_obs - tn_samp + 1
            for k in tn:
                if k == tmax:
                    tpmf.append(1. / tn_samp)
                elif k == 1:
                    tpmf.append((tn_samp - 1) / tn_samp)
                else:
                    tpmf.append(0)

            pmf.append(np.array(tpmf))

        return pmf
            
    @doc_inherit
    def rad(self):
        # Get parameters
        n_samp, tot_obs = self.get_params(['n_samp', 'tot_obs' ])

        # Calculate rad
        rad = []
        for tn_samp, ttot_obs in zip(n_samp, tot_obs):
            trad = [ttot_obs - tn_samp + 1]
            trad = np.array(list(np.repeat(1, tn_samp - 1)) + trad)
            rad.append(trad)

        return rad


class logser(Distribution):
    __doc__ = Distribution.__doc__ + \
    '''
    Description
    -----------
    Fisher's log series distribution (Fisher et al. 1943). Also known as the 
    logarithmic distribution.

    Parameters
    ----------
    n_samp : int or iterable
        Total number of species / samples
    tot_obs: int or iterable
        Total number of individuals / observations

    self.var keywords
    ------------------
    p : list of floats
        p parameter of standard logseries distribution

    Notes
    -----
    To use a known mean of the distribution as the parameter, set n_samp = 1 
    and tot_obs = mean.

    The total species (S) is equivalent to n_samp and the total
    individuals (N) is equivalent to tot_obs.
    '''
    
    @doc_inherit
    def __init__(self, **kwargs):
        self.params = kwargs
        self.min_supp = 1
        self.par_num = 1
        self.var = {}
    
    @doc_inherit
    def pmf(self, n):
        
        # Get parameters
        n_samp, tot_obs = self.get_params(['n_samp', 'tot_obs'])
        n = expand_n(n, len(n_samp))
        
        # TODO: Additional checks?
        assert np.all(n_samp <= tot_obs), 'n_samp must be <= tot_obs'
        
        # Calculate pmf
        stop = 1 - 1e-10
        start = -2
        eq = lambda x, n_samp, tot_obs: (((tot_obs/x) - tot_obs) * 
                                                (-(np.log(1 - x)))) - n_samp

        pmf = []
        self.var['p'] = []

        for tn_samp, ttot_obs, tn in zip(n_samp, tot_obs, n):
            # Catching cryptic brentq error
            try:
                tp = scipy.optimize.brentq(eq, start, stop, 
                                            args=(tn_samp,ttot_obs), disp=True)
            except(ValueError):
                raise ValueError("No solution to %s.pmf when tot_obs = %.2f"\
                                  % (self.__class__.__name__, ttot_obs) + 
                                  " and n_samp = %.2f" % (tn_samp)) 
            tpmf = stats.logser.pmf(tn, tp)
            self.var['p'].append(tp)
            pmf.append(tpmf)
        
        self.var['p'] = np.array(self.var['p'])
        return pmf

    @doc_inherit
    def cdf(self, n):
        
        # Get parameters
        n_samp, tot_obs = self.get_params(['n_samp', 'tot_obs'])
        n = expand_n(n, len(n_samp))
        
        # TODO: Additional checks?
        assert np.all(n_samp <= tot_obs), 'n_samp must be <= tot_obs'
        
        # Calculate pmf
        stop = 1 - 1e-10
        start = -2
        eq = lambda x, n_samp, tot_obs: (((tot_obs/x) - tot_obs) * 
                                                (-(np.log(1 - x)))) - n_samp

        cdf = []
        self.var['p'] = []

        for tn_samp, ttot_obs, tn in zip(n_samp, tot_obs, n):
            # Catching cryptic brentq error
            try:
                tp = scipy.optimize.brentq(eq, start, stop, 
                                            args=(tn_samp,ttot_obs), disp=True)
            except(ValueError):
                raise ValueError("No solution to %s.cdf when tot_obs = %.2f"\
                                  % (self.__class__.__name__, ttot_obs) + 
                                  " and n_samp = %.2f" % (tn_samp)) 
            tcdf = stats.logser.cdf(tn, tp)
            self.var['p'].append(tp)
            cdf.append(tcdf)
   
        self.var['p'] = np.array(self.var['p'])

        return cdf

class logser_ut(Distribution):
    __doc__ = Distribution.__doc__ + \
    '''
    Description
    -----------
    Upper-truncated log series (Harte et al 2008, Harte 2011). Like Fisher's, 
    but normalized to finite support.

    Parameters
    ----------
    n_samp : int or iterable
        Total number of species / samples
    tot_obs: int or iterable
        Total number of individuals / observations
   
    self.vars keywords
    ------------------
    x : list of floats
         exp(-beta) parameter. Beta is lagrange multiplier

    Notes
    -----
    This distribution is the truncated logseries described in Eq 7.32 of Harte 
    2011. Eq. 7.27 is used to solve for the Lagrange multiplier.

    Realistic values of x where x = e**(-beta) are in the range (1/e, 1). The 
    start and stop parameters for the brentq procedure are set close to these 
    values. However, x can occasionally be greater than one, so the maximum 
    stop value of the brentq optimizer is 2.

    The total species (S) is equivalent to n_samp and the total
    individuals (N) is equivalent to tot_obs.
    '''
    
    @doc_inherit
    def __init__(self, **kwargs):
        self.params = kwargs
        self.min_supp = 1
        self.par_num = 2 # This is highly contested
        self.var = {}

    @doc_inherit    
    def pmf(self, n):

        # Get parameters
        n_samp, tot_obs = self.get_params(['n_samp', 'tot_obs'])
        n = expand_n(n, len(n_samp))
        
        # TODO: Additional checks?
        assert np.all(n_samp <= tot_obs), 'n_samp must be <= tot_obs'
        

        # Calculate pmf
        start = 0.3
        stop = 2
        flmax = sys.float_info[0]

        pmf = []
        self.var['x'] = []

        for tn_samp, ttot_obs, tn in zip(n_samp, tot_obs, n):

            # If n_samp = tot_obs, return 1 for n = 1 and 0 otherwise 
            # (e**-beta = 0)
            if tn_samp == ttot_obs:
                tpmf = np.zeros(len(tn))
                tpmf[tn == 1] = 1
                tx = 0

            else:
                k = np.linspace(1, ttot_obs, num=ttot_obs)
                try:
                    tx = scipy.optimize.brentq(beta_solver, start,
                               min((flmax/tn_samp)**(1/float(ttot_obs)), stop), 
                               args = (k, ttot_obs, tn_samp), disp=True)
                except(ValueError):
                    raise ValueError("No solution to %s.pmf when tot_obs = "
                                  % (self.__class__.__name__) + 
                                  "%.2f and n_samp = %.2f" % (ttot_obs, tn_samp))
                tnorm = np.sum(tx ** k / k)
                tpmf = (tx ** tn / tn) / tnorm

            self.var['x'].append(tx)
            pmf.append(tpmf)
   
        self.var['x'] = np.array(self.var['x'])
        return pmf

    # TODO: Add exact cdf from JK dissertation


class logser_ut_appx(Distribution):
    __doc__ = Distribution.__doc__ + \
    '''
    Description
    -----------
    Upper-truncated log series (Harte et al 2008, Harte 2011). Like Fisher's, 
    but normalized to finite support.  Uses approximation from Harte (2011).

    Parameters
    ----------
    n_samp : int or iterable
        Total number of species / samples
    tot_obs: int or iterable
        Total number of individuals / observations
   
    self.var keywords
    ------------------
    x : list of floats
         exp(-beta) parameter. Beta is lagrange multiplier
        
    Notes:
    ------
    This distribution is the truncated logseries described in Eq 7.32 of Harte 
    2011. Eq. 7.30, the approximate equation, is used to solve for the Lagrange 
    multiplier. This class is faster than the logser_ut with no approximation.

    Realistic values of x where x = e**(-beta) are in the range (1/e, 1). The 
    start and stop parameters for the brentq procedure are set close to these 
    values. However, x can occasionally be greater than one, so the maximum 
    stop value of the brentq optimizer is 2.

    The pmf method has an internal value named root. In the 
    approximation equation, there are two roots (solutions) in the solution for 
    the lagrange multiplier.  Root 2 is the root typically used in calculations 
    and is the default.  If root=1, the first root will be used and this is not 
    a true pmf. Root can only be changes within the code.

    The total species (S) is equivalent to n_samp and the total
    individuals (N) is equivalent to tot_obs.
    '''
    
    @doc_inherit
    def __init__(self, **kwargs):
        self.params = kwargs
        self.min_supp = 1
        self.par_num = 2 # This is highly contested
        self.var = {}
    
    @doc_inherit
    def pmf(self, n):
        
        # Multiple roots. root = 2 makes it a logseries
        root = 2

        # Get parameters
        n_samp, tot_obs = self.get_params(['n_samp', 'tot_obs'])
        n = expand_n(n, len(n_samp))

        # TODO: Additional Checks
        assert np.all(n_samp <= tot_obs), 'n_samp must be <= tot_obs'

        # Calculate pmf
        start = 0.3
        stop = 1 - 1e-10
        eq = lambda x, n_samp, tot_obs: (((-m.log(x))*(m.log(-1/(m.log(x))))) - 
                                                       (float(n_samp)/tot_obs))
        pmf = []
        self.var['x'] = []

        for tn_samp, ttot_obs, tn in zip(n_samp, tot_obs, n):
            
            # TODO: What if tot_obs = n_samp? 
            if tn_samp == ttot_obs:
                tpmf = np.zeros(len(tn))
                tpmf[tn == 1] = 1
                tx = 0
            else:
                # Try normal root finder. Will fail if two roots
                try:
                    tx = scipy.optimize.brentq(eq, start, stop, 
                                            args=(tn_samp, ttot_obs),disp=True)

                # If that fails, try a more complex decision tree
                except ValueError:
                    eq1 = lambda x: -1 * eq(x, tn_samp, ttot_obs)
                    xmax = scipy.optimize.fmin(eq1, .5, disp=0)[0]
                    ymax = eq(xmax, tn_samp, ttot_obs)
                    if ymax > 0:
                        if root == 1:
                            tx = scipy.optimize.brentq(eq, start, xmax,
                                           args=(tn_samp, ttot_obs), disp=True)
                        if root == 2:
                            tx = scipy.optimize.brentq(eq, xmax, stop, 
                                           args=(tn_samp, ttot_obs), disp=True)
                    if ymax < 0:
                        raise ValueError('No solution to ' +
                            ' %s.pmf' % (self.__class__.__name__) +
                            ' when tot_obs = %.2f and n_samp = %.2f ' % 
                            (ttot_obs, tn_samp)) 

                g = -1/np.log(tx)
                tpmf = (1/np.log(g)) * ((tx**tn)/tn)

            self.var['x'].append(tx)
            pmf.append(tpmf)

        self.var['x'] = np.array(self.var['x'])

        return pmf


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
    n_samp : int or iterable (optional)
        Total number of species / samples
    tot_obs: int or iterable (optional)
        Total number of individuals / observations

    self.var keywords
    -----------------
    mu : float or iterable
        The mu parameter of the poisson log normal
    sigma : float or iterable
        The sigma parameter of the poisson log normal

    These parameters are stored in var as well as params if they are calculated
    with the fit method.

    Notes
    -----
    The pmf method was adopted directly from the VGAM package in R by Mark
    Wilber. The VGAM R package was adopted directly from Bulmer (1974). The fit 
    function was adapted from Ethan White's pln_solver function in 
    weecology.

    The total species (S) is equivalent to n_samp and the total
    individuals (N) is equivalent to tot_obs.
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
        self.var = {}

    
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

        See class docstring for more specific information on this distribution.
        '''

        # Get parameters
        mu, sigma = self.get_params(['mu', 'sigma'])
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

        return pmf

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

        super(plognorm, self).fit(data)
        data = check_list_of_iterables(data)

        # Calculate and store parameters
        temp_mu = []
        temp_sigma = []

        for tdata in data:
            mu0 = np.mean(np.log(tdata))  # Starting guesses for mu and sigma
            sigma0 = np.std(np.log(tdata), ddof=1)
            
            # TODO: Can we do this without setting the self.params? Make
            # another plognorm inside?
            def pln_func(x):
                self.params['mu'] = x[0]
                self.params['sigma'] = x[1]
                return -sum(np.log(self.pmf(tdata)[0]))

            mu, sigma = scipy.optimize.fmin(pln_func, x0=[mu0, sigma0],
                                            disp=0)
            temp_mu.append(mu)
            temp_sigma.append(sigma)

        self.params['mu'] = np.array(temp_mu)
        self.params['sigma'] = np.array(temp_sigma)
        self.var['mu'] = np.array(temp_mu)
        self.var['sigma'] = np.array(temp_sigma)

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
    n_samp : int or iterable (optional)
        Total number of species / samples
    tot_obs: int or iterable (optional)
        Total number of individuals / observations

    self.var keywords
    ------------------
    mu : float or iterable
        The mu parameter of the poisson log normal
    sigma : float or iterable
        The sigma parameter of the poisson log normal

    These parameters are stored in var as well as params if they are calculated
    with the fit method.

    Notes
    -----
    The pmf method was adopted directly from the VGAM package in R by Mark
    Wilber. The VGAM R package was adopted directly from Bulmer (1974). The fit 
    function was adapted from Ethan White's pln_solver function in weecology.

    Truncation calculation based on Bulmer Eq. A1.

    The total species (S) is equivalent to n_samp and the total
    individuals (N) is equivalent to tot_obs.

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
        self.var = {}
        
    
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

        See class docstring for more specific information on this distribution.
        '''

        # Get parameters
        mu, sigma = self.get_params(['mu', 'sigma'])

        # TODO: Additional parameter checks

        # Calculate pmf, using plognorm as aid
        reg_plog = plognorm(mu=mu, sigma=sigma)
        reg_pmf = reg_plog.pmf(n)
        reg_pmf0 = reg_plog.pmf(0)
        self.var = reg_plog.var

        trunc_pmf = [(pr / (1 - p0)) for pr, p0 in zip(reg_pmf, reg_pmf0)]

        return trunc_pmf 

    # TODO: Write cdf method based on cdf of plognorm, similar to above

class canonical_lognorm(Distribution):
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
    n_samp : int or iterable (optional)
        Total number of species / samples
    tot_obs: int or iterable (optional)
        Total number of individuals / observations

    self.var keywords
    -----------------
    mu : list of floats
        The mu parameter of the lognormal calculated with
        np.log(tot_obs / n_samp) - (sigma**2 / 2).
    sigma : list of float
        The sigma parameter of the log normal

    Notes
    -----
    Currently, lognormal is implemented so that mu is calculated using tot_obs,
    n_samp, and sigma.  While, mu can be passed in as a keyword argument, this
    mu will be ignored. 
        
    '''
    
    @doc_inherit
    def __init__(self, **kwargs):
        self.params = kwargs
        self.min_supp = 1
        self.par_num = 2
        self.var = {}

    @doc_inherit  
    def pmf(self, n):

        # Get parameters
        tot_obs, n_samp = self.get_params(['tot_obs','n_samp'])
        n = expand_n(n, len(tot_obs))

        # Calculate sigma
        sigma = np.sqrt((2 * np.log(n_samp)) / np.log(2)**2)

        # Calculate mu
        mu = np.log(tot_obs / n_samp) - (sigma**2 / 2)
        self.var['mu'] = mu
        self.var['sigma'] = sigma

        # Calculate pmf
        pmf = []
        for tmu, tsigma, tn in zip(mu, sigma, n):
            tpmf = stats.lognorm.pdf(tn, tsigma, scale=np.exp(tmu))
            pmf.append(tpmf)

        return pmf

    @doc_inherit  
    def cdf(self, n):

        # Get parameters
        tot_obs, n_samp = self.get_params(['tot_obs','n_samp'])
        n = expand_n(n, len(tot_obs))

        # Calculate sigma
        sigma = np.sqrt((2 * np.log(n_samp)) / np.log(2)**2)

        # Calculate mu
        mu = np.log(tot_obs / n_samp) - (sigma**2 / 2)
        self.var['mu'] = mu
        self.var['sigma'] = sigma

        #Calculate cdf
        cdf = []
        for tmu, tsigma, tn in zip(mu, sigma, n):
            tcdf = stats.lognorm.cdf(tn, tsigma, scale=np.exp(tmu))
            cdf.append(tcdf)

        return cdf


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
    n_samp : int or iterable (optional)
        Total number of species / samples
    tot_obs: int or iterable (optional)
        Total number of individuals / observations

    self.var keywords
    -----------------
    mu : list of floats
        The mu parameter of the lognormal calculated with
        np.log(tot_obs / n_samp) - (sigma**2 / 2).
    sigma : list of float
        The sigma parameter of the log normal

    Notes
    -----
    Currently, lognormal is implemented so that mu is calculated using tot_obs,
    n_samp, and sigma.  While, mu can be passed in as a keyword argument, this
    mu will be ignored. 
        
    '''
    
    @doc_inherit
    def __init__(self, **kwargs):
        self.params = kwargs
        self.min_supp = 1
        self.par_num = 2
        self.var = {}

    @doc_inherit  
    def pmf(self, n):

        # Get parameters
        tot_obs, n_samp, sigma = self.get_params(['tot_obs','n_samp','sigma'])
        
        #mu, sigma = self.get_params(['mu','sigma'])
        # TODO: Additional parameter checks
        n = expand_n(n, len(sigma))
        
        # Calculate mu
        mu = np.log(tot_obs / n_samp) - (sigma**2 / 2)
        self.var['mu'] = mu
        self.var['sigma'] = sigma

        # Calculate pmf
        pmf = []
        for tmu, tsigma, tn in zip(mu, sigma, n):
            tpmf = stats.lognorm.pdf(tn, tsigma, scale=np.exp(tmu))
            pmf.append(tpmf)

        return pmf

    @doc_inherit  
    def cdf(self, n):

        # Get parameters
        tot_obs, n_samp, sigma = self.get_params(['tot_obs','n_samp','sigma'])
        
        # TODO: Additional parameter checks
        n = expand_n(n, len(tot_obs))

        # Calculate mu
        mu = np.log(tot_obs / n_samp) - (sigma**2 / 2)
        self.var['mu'] = mu

        #Calculate cdf
        cdf = []
        for tmu, tsigma, tn in zip(mu, sigma, n):
            tcdf = stats.lognorm.cdf(tn, tsigma, scale=np.exp(tmu))
            cdf.append(tcdf)

        return cdf

    @doc_inherit 
    def fit(self, data):

        super(lognorm, self).fit(data)
        n_samp, tot_obs = self.get_params(['n_samp', 'tot_obs'])

        data = check_list_of_iterables(data) 
        tempsig = []

        for tdata, tn_samp, ttot_obs in zip(data, n_samp, tot_obs): 

            def ln_func(sigma):
                self.params['tot_obs'] = ttot_obs
                self.params['n_samp'] = tn_samp
                self.params['sigma'] = sigma 
                return -sum(np.log(self.pmf(tdata)[0]))

            mle_sigma = scipy.optimize.fmin(ln_func,
                        np.array([np.std(np.log(tdata), ddof=1)]), disp=0)[0]
            tempsig.append(mle_sigma)

        self.params['sigma'] = np.array(tempsig)
        self.var['sigma'] = np.array(tempsig)
        self.params['n_samp'] = n_samp
        self.params['tot_obs'] = tot_obs

        return self

class dgamma(Distribution):
    __doc__ = Distribution.__doc__ + \
    '''
    Discrete Gamma distribution as specified by Ugland et al. 2007 and Frank et
    al 2011.

    Parameters
    ----------
    n_samp : int or iterable
        Total number of species / samples
    tot_obs: int or iterable
        Total number of individuals / observations
    alpha : float or iterable 
        The alpha parameter of the discrete gamma distribution
    theta : float or iterable
        The theta parameter of the discrete gamma distribution

    self.var keywords
    -----------------
    alpha : float or iterable 
        The alpha parameter of the discrete gamma distribution
    theta : float or iterable
        The theta parameter of the discrete gamma distribution

    Alpha and theta are included in vars if they are calculated from the fit
    function.


    '''

    @doc_inherit
    def __init__(self, **kwargs):
        self.params = kwargs
        self.min_supp = 1
        self.par_num = 2 
        self.var = {}

    @doc_inherit
    def pmf(self, n):

        # Get parameters
        n_samp, tot_obs, alpha, theta =\
                self.get_params(['n_samp', 'tot_obs', 'alpha', 'theta'])
        n = expand_n(n, len(n_samp))

        assert np.all(n_samp <= tot_obs), 'n_samp must be <= tot_obs'

        # Calculate pmf
        eq = lambda x, alpha, theta: x**(alpha - 1) * theta**x       

        pmf = []
        for tn_samp, ttot_obs, talpha, ttheta, tn in zip(n_samp, tot_obs,
                                                alpha, theta, n):
            # Normalization constant
            sumg = sum(eq(np.arange(1, np.floor(ttot_obs) + 1), talpha,
                                                                    ttheta))
            tpmf = eq(tn, talpha, ttheta) / sumg # Normalizing
            pmf.append(tpmf)

        return pmf

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

        super(dgamma, self).fit(data)
        data = check_list_of_iterables(data)

        # Calculate and store parameters
        temp_alpha = []
        temp_theta = []

        for tdata in data:
            alpha0 = 1 # starting guesses for alpha and theta 
            theta0 = .9 
            
            def dgm_func(x):
                self.params['alpha'] = x[0]
                self.params['theta'] = x[1]
                return -sum(np.log(self.pmf(tdata)[0]))

            alpha, theta = scipy.optimize.fmin(dgm_func, x0=[alpha0, theta0],
                                            disp=0)
            temp_alpha.append(alpha)
            temp_theta.append(theta)

        self.params['alpha'] = np.array(temp_alpha)
        self.params['theta'] = np.array(temp_theta)
        self.var['alpha'] = np.array(temp_alpha)
        self.var['theta'] = np.array(temp_theta)

        return self


class geo_ser(Distribution):
    __doc__ = Distribution.__doc__ + \
    '''
    Geometric series distribution (Motomura 1932 and Magurran 1988).

    Parameters
    ----------
    n_samp : int or iterable
        Total number of species / samples
    tot_obs: int or iterable
        Total number of individuals / observations
    k : float
        The fraction of resources that each species acquires. Range is 
        (0, 1].

    self.var keywords
    -----------------
    k : array of floats
        The fraction of resources that each species acquires. Range is 
        (0, 1].
        
        k is only in vars if it is calculated with the fit function
    
    
    Notes
    -----
    Equation for pmf and fit from May (1975).  Note that the pmf derived by May
    is continuous. In this case, the pmf given in the .pmf() method DOES NOT
    sum to one. If we normalize the distribution, the super(geo_ser,
    self).rad() does not return a rank abundance distribution with the proper
    tot_obs and the the value of fitted value of'k' is inherently changed.
    Because our fitting method for 'k' has been tested and verified, we chose
    not to change the value of 'k' computed in geo_ser.fit().  Therefore, we
    recommend that the rank abundance distribution (geo_ser.rad()) is used for
    serious analyses with geo_ser.  

    The total species (S) is equivalent to n_samp and the total
    individuals (N) is equivalent to tot_obs.

    '''

    @doc_inherit
    def __init__(self, **kwargs):
        self.params = kwargs
        self.min_supp = 1
        self.par_num = 3 # May says 2 parameters, test this
        self.var = {}

    # The pmf, as adapted from May, describes a continuous SAD. The rad
    # generate by the below pmf (super.rad) does not give the rad predicted by
    # geo_ser.rad.
    
    @doc_inherit
    def pmf(self, n):

        # Get parameters

        n_samp, tot_obs, k = self.get_params(['n_samp', 'tot_obs', 'k'])
        n = expand_n(n, len(n_samp))
        assert np.all(n_samp <= tot_obs), 'n_samp must be <= tot_obs'
        assert np.all(k > 0) and np.all(k <= 1), ('k must be in the '
                                                  'interval (0, 1]')

        # Calculate pmf
        pmf = []

        # Equation from May 1975.  
        eq = lambda x, n_samp, k: (1 / x) * (1 / n_samp) * (1 / np.log(1 / 
                                                                (1 - k)))
        ### Other possible equations for calculating the geo_ser pmf ###
        # Normalization for continuous pdf
        #kfxn = lambda N, S: (np.exp((np.log(N) / S)) - 1) / np.exp(np.log(N) /
        #                                                                     S)
        #Normalization for discrete pmf
        #keq = lambda k, N, S: 1  - (((1/S) * (1 / (np.log(1 / (1 - k))))) * sum((1 /
        #                                                 np.arange(1, N + 1))))

        #self.var['k'] = []
        #eq = lambda x, n_samp, k: (-1 / n_samp) * (1 / x) * (1 / np.log(1 - k))

        ### END ###

        for tn_samp, ttot_obs, tk, tn in zip(n_samp, tot_obs, k, n):

            ### Other possible ways in which to calculate the geo_ser pmf ###

            #c = 1 / np.log(ttot_obs)
            #sumg = sum( c / np.arange(1, np.floor(ttot_obs) + 1))
            #tpmf = (c / tn) / sumg # normalize
            #tk = kfxn(ttot_obs, tn_samp)
            #tk = scipy.optimize.brentq(keq, 1e-10, 1 - 1e-10, args=(ttot_obs, tn_samp),
            #                                                        disp=True)
            #self.var['k'] = tk
            #sumg = sum(eq(np.arange(1, ttot_obs + 1), tn_samp, tk))

            ### END ###

            tpmf = eq(tn, tn_samp, tk) #/ sumg
            pmf.append(tpmf)

        return pmf

    @doc_inherit
    def rad(self):

        # Get parameters
        n_samp, tot_obs, k = self.get_params(['n_samp', 'tot_obs', 'k'])
        assert np.all(k > 0) and np.all(k <= 1), ('k must be in the ' + 
                                                  'interval (0, 1]')

        # Calculate rad
        rad = []
        for tn_samp, ttot_obs, tk in zip(n_samp, tot_obs, k):
            tn_samp = np.round(tn_samp, decimals=0)
            tC = (1 - (1 - tk ) ** tn_samp) ** - 1
            trad = ttot_obs * tC * tk * (1 - tk) ** (np.arange(1, tn_samp + 1) 
                                                                           - 1)
            rad.append(np.sort(trad))

        return rad


    @doc_inherit
    def fit(self, data):

        # Get parameters
        super(geo_ser, self).fit(data)  # Run Distribution.fit method
        n_samp = self.params['n_samp']
        tot_obs = self.params['tot_obs']

        # TODO: Additional checks?
        assert np.all(n_samp <= tot_obs), 'n_samp must be <= tot_obs'

        # Calculate fit
        self.params['k'] = []
        for tdata, tn_samp, ttot_obs in zip(data, n_samp, tot_obs):
            ttot_obs_min = np.min(tdata)
            eq = lambda x: (((x / (1 - x)) *
                             ((1 - x) ** tn_samp / (1 - (1 - x) ** tn_samp)))
                            - (ttot_obs_min / ttot_obs))
            try:
                tk = scipy.optimize.brentq(eq, 1e-10, 1 - 1e-10, disp=True)
            except(ValueError):
                raise ValueError("No solution for k in %s.fit with tot_obs = "\
                                 % (self.__class__.__name__) + 
                                 "%.2f and n_samp = %.2f" % (ttot_obs, tn_samp))
            self.params['k'].append(tk)
        self.params['k'] = np.array(self.params['k'])
        self.var['k'] = deepcopy(self.params['k'])
        return self


class broken_stick(Distribution):
    __doc__ = Distribution.__doc__ + \
    '''
    McArthur's broken stick distribution (May 1975)

    Parameters
    ----------
    n_samp : int or iterable
        Total number of species / samples
    tot_obs: int or iterable
        Total number of individuals / observations

    self.var keywords
    -----------------
    None

    The total species (S) is equivalent to n_samp and the total
    individuals (N) is equivalent to tot_obs.

    '''

    @doc_inherit
    def __init__(self, **kwargs):
        self.params = kwargs
        self.min_supp = 1
        self.par_num = 2  # May says 1
        self.var = {}


    @doc_inherit
    def pmf(self, n):
        # TODO:  PMF is not quite summing to one. But it is checking against
        # known results. See test_distributions. 
        
        # Get parameters
        n_samp, tot_obs = self.get_params(['n_samp', 'tot_obs'])
        n = expand_n(n, len(n_samp))

        assert np.all(n_samp <= tot_obs), 'n_samp must be <= tot_obs'
        
        # Calculate pmf
        eq = lambda x, n_samp, tot_obs: ((n_samp - 1) / tot_obs) * \
                                          ((1 - (x / tot_obs)) ** (n_samp - 2))
        pmf = []
        for tn_samp, ttot_obs, tn in zip(n_samp, tot_obs, n):
            ttot_obs = np.round(ttot_obs, decimals=0)
            #sumg = sum(eq(np.arange(1, np.floor(ttot_obs) + 1), tn_samp, ttot_obs))
            tpmf = eq(tn, tn_samp, ttot_obs) #/ sumg # Normalizing
            pmf.append(tpmf)

        return pmf


    @doc_inherit
    def rad(self):
        
        # Get parameters
        n_samp, tot_obs = self.get_params(['n_samp', 'tot_obs'])

        # TODO: Additional checks?
        assert np.all(n_samp <= tot_obs), 'n_samp must be <= tot_obs'

        # Calculate rad
        rad = []
        for tn_samp, ttot_obs in zip(n_samp, tot_obs):
            tn_samp = np.round(tn_samp, decimals=0)
            trad = np.empty(tn_samp)
            for i in np.arange(tn_samp):
                n = np.arange(i + 1, tn_samp + 1) 
                trad[i] = (ttot_obs / tn_samp) * sum(1 / n)
            rad.append(np.sort(trad))

        return rad
            

class sugihara(Distribution):
    __doc__ = Distribution.__doc__ + \
    '''
    Description
    -----------
    Sugihara Rank Abundance Distribution (Sugihara 1980)

    Parameters
    ----------
    n_samp : int or iterable
        Total number of species / samples
    tot_obs: int or iterable
        Total number of individuals / observations

    self.var keyword
    ----------------
    None

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
    
    The total species (S) is equivalent to n_samp and the total
    individuals (N) is equivalent to tot_obs.
    
    '''
    # TODO: Back-derive pmf?
    
    @doc_inherit
    def __init__(self, **kwargs):
        self.params = kwargs
        self.min_supp = 1
        self.par_num = 1
        self.var = {}
    

    @doc_inherit
    def rad(self, sample_size=10000):
        
        # Get parameters
        n_samp, tot_obs = self.get_params(['n_samp', 'tot_obs'])
        assert np.all(n_samp <= tot_obs), 'n_samp must be <= tot_obs'

        # Calculate rad
        rad = []
        for tn_samp, ttot_obs in zip(n_samp, tot_obs):
            total = []
            for i in xrange(sample_size):
                U = np.random.triangular(0.5, 0.75, 1, size=tn_samp - 1)
                p = []
                # TODO: Could this be refactored to perform better?
                for i in xrange(np.int(tn_samp)):
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
            means = np.array([np.mean(total_array[:,i]) for i in 
                                                         xrange(int(tn_samp))])
            rad.append(np.sort(ttot_obs * means))

        return rad

    def cdf(self, n):
        '''
        No cdf exists for this distribution

        '''

        raise NotImplementedError('No CDF exists for object %s' %
                                                    self.__class__.__name__)


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

    self.var keywords
    ------------------
    p : list of floats
        p parameter of is equal to 1 / n_samp, the 'p' parameter of the
        binomial distribution

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
        self.var = {}
    
    @doc_inherit
    def pmf(self, n):
        # Get parameters
        n_samp, tot_obs = self.get_params(['n_samp', 'tot_obs'])
        n = expand_n(n, len(n_samp))
        
        # TODO: Additional checks?
        
        pmf = []
        self.var['p'] = []
        for tn_samp, ttot_obs, tn in zip(n_samp, tot_obs, n):
            ta = 1 / tn_samp
            pmf.append(stats.binom.pmf(tn, ttot_obs, ta))
            self.var['p'].append(ta)
        return pmf
    
    @doc_inherit
    def cdf(self, n):
        n_samp, tot_obs = self.get_params(['n_samp', 'tot_obs'])
        n = expand_n(n, len(n_samp))
        
        # TODO: Additional checks?

        cdf = []
        self.var['p'] = []
        for tn_samp, ttot_obs, tn in zip(n_samp, tot_obs, n):
            ta = 1 / tn_samp
            cdf.append(stats.binom.cdf(tn, ttot_obs, ta))
            self.var['p'].append(ta)
        return cdf

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

    self.var keywords
    -----------------
    mu : list of floats
        the mu parameter of the poisson distribution

    '''

    @doc_inherit
    def __init__(self, **kwargs):
        self.params = kwargs
        self.min_supp = 0
        self.par_num = 1
        self.var = {}
    
    @doc_inherit
    def pmf(self, n):

        n_samp, tot_obs = self.get_params(['n_samp', 'tot_obs'])
        n = expand_n(n, len(n_samp))
        
        # TODO: Additional checks?
        
        pmf = []
        self.var['mu'] = []
        for tn_samp, ttot_obs, tn in zip(n_samp, tot_obs, n):
            tmu = ttot_obs * (1 / tn_samp)
            pmf.append(stats.poisson.pmf(tn, tmu))
            self.var['mu'].append(tmu)
        return pmf
    
    @doc_inherit
    def cdf(self, n): 
        
        n_samp, tot_obs = self.get_params(['n_samp', 'tot_obs'])
        n = expand_n(n, len(n_samp))
        
        # TODO: Additional checks?

        cdf = []
        self.var['mu'] = []
        for tn_samp, ttot_obs, tn in zip(n_samp, tot_obs, n):
            tmu = ttot_obs * (1 / tn_samp)
            cdf.append(stats.poisson.cdf(tn, tmu))
            self.var['mu'].append(tmu)
        return cdf

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

    self.var keywords
    -----------------
    Parameterization differs for different forms of the nbd.  We use the
    standard ecological form as described by Ben Bolker. Parameters 'a' (1 /
    n_samp), 'tot_obs', and k are used to derive the nbd parameter p (see code
    for details).  Parameters k and p are used to generate distribution. k is
    included in self.var if it is calculated in fit.
        
    p : array of floats 
        p parameters of nbd
    k : array of floats
        Aggregation parameter
    '''
    
    @doc_inherit
    def __init__(self, **kwargs):
        self.params = kwargs
        self.min_supp = 0
        self.par_num = 2
        self.var = {}
    
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

        See class docstring for more specific information on this distribution.
        '''
        
        n_samp, tot_obs, k = self.get_params(['n_samp', 'tot_obs', 'k'])
        n = expand_n(n, len(n_samp))
        
        # TODO: Additional checks?
        
        pmf = []
        self.var['p'] = []

        for tn_samp, ttot_obs, tk, tn in zip(n_samp, tot_obs, k, n):
            tmu = ttot_obs * (1 / tn_samp)
            tp = 1 / (tmu / tk + 1) # See Bolker book Chapt 4
            pmf.append(scipy.stats.nbinom.pmf(tn, tk, tp))
            self.var['p'].append(tp)

        self.var['p'] = np.array(self.var['p'])

        return pmf 

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

        See class docstring for more specific information on this distribution.
        '''

        n_samp, tot_obs, k = self.get_params(['n_samp', 'tot_obs', 'k'])
        n = expand_n(n, len(n_samp))
        
        # TODO: Additional checks?
        
        cdf = []
        self.var['p'] = []

        for tn_samp, ttot_obs, tk, tn in zip(n_samp, tot_obs, k, n):
            tmu = ttot_obs * (1 / tn_samp)
            tp = 1 / (tmu / tk + 1) # See Bolker book Chapt 4
            cdf.append(scipy.stats.nbinom.cdf(tn, tk, tp))
            self.var['p'].append(tp)
        
        self.var['p'] = np.array(self.var['p'])
        return cdf
    
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
        n_samp, tot_obs = self.get_params(['n_samp', 'tot_obs'])

        data = check_list_of_iterables(data) 
        tempk = []

        for tdata, tn_samp, ttot_obs in zip(data, n_samp, tot_obs): 

            def nll_nb(k):
                self.params['tot_obs'] = ttot_obs
                self.params['n_samp'] = tn_samp
                self.params['k'] = k
                return -sum(np.log(self.pmf(tdata)[0]))

            mlek = scipy.optimize.fmin(nll_nb, np.array([guess_for_k]), 
                                                                    disp=0)[0]
            tempk.append(mlek)
        self.params['k'] = np.array(tempk)
        self.params['n_samp'] = n_samp
        self.params['tot_obs'] = tot_obs
        self.var['k'] = np.array(tempk)

        return self


class nbd_lt(Distribution):
    __doc__ = Distribution.__doc__ + \
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
    k : float
        Aggregation parameter

    self.var keywords
    -----------------
    p : array of floats 
        p parameters of nbd_lt
    k : array of floats
        Aggregation parameter
        k is included in self.var if it is calculated in fit.

    Notes
    -----
    The total species (S) is equivalent to n_samp and the total
    individuals (N) is equivalent to tot_obs.

    Parameterization based on Sampford 1955 and He and Legendre 2002


    '''

    @doc_inherit
    def __init__(self, **kwargs):
        self.params = kwargs
        self.min_supp = 1
        self.par_num = 2 
        self.var = {}

    def pmf(self, n):
        """
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

        See class docstring for more specific information on this distribution.


        """

        # Get parameters
        n_samp, tot_obs, k =\
                self.get_params(['n_samp', 'tot_obs', 'k'])
        n = expand_n(n, len(n_samp))

        assert np.all(n_samp <= tot_obs), 'n_samp must be <= tot_obs'

        # Calculate pmf
        def pmf_eq(n, p, k):

            norm = np.exp(spec.gammaln(k + n) - ((spec.gammaln(k) + 
                                            spec.gammaln(n + 1))))
            
            kernel = (p / (1 + p))**n * (1 / ((1 + p)**k - 1))
            return norm * kernel

        self.var['p'] = []

        pmf = []
        p_eq = lambda p, k, N, S : (k * p) / (1 - (1 + p)**-k) -\
                                    (float(N) / S) 

        for tn_samp, ttot_obs, tk, tn in zip(n_samp, tot_obs, k, n):
            # Find p
            
            do_it = True
            count = 0
            while do_it and count < 20:

                stop = 10**(count + 1)
                count += 1

                try:
                    tp = scipy.optimize.brentq(p_eq, 1e-10, stop, args=(tk,
                                ttot_obs, tn_samp))
                    do_it = False
                except(ValueError):
                    if count >= 20:
                        tp = np.nan
                    
            self.var['p'].append(tp)
            tpmf = pmf_eq(tn, tp, tk)

            pmf.append(tpmf)
     
        self.var['p'] = np.array(self.var['p'])
        return pmf

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
            List of 1D arrays of cumulative probability of observing sample n.

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

        pmf_list = self.pmf(n_in)

        # Calculate cdfs
        cdf = []
        for tpmf, tn in zip(pmf_list, n):
            full_cdf = np.cumsum(tpmf)
            tcdf = np.array([full_cdf[x - self.min_supp] for x in tn])
            cdf.append(tcdf)

        return cdf 

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

        super(nbd_lt, self).fit(data)
        n_samp, tot_obs = self.get_params(['n_samp', 'tot_obs'])

        data = check_list_of_iterables(data) 
        tempk = []

        for tdata, tn_samp, ttot_obs in zip(data, n_samp, tot_obs): 

            def nll_nb(k):
                self.params['tot_obs'] = ttot_obs
                self.params['n_samp'] = tn_samp
                self.params['k'] = k
                return -sum(np.log(self.pmf(tdata)[0]))

            mlek = scipy.optimize.fmin(nll_nb, np.array([guess_for_k]), 
                                                                    disp=0)[0]
            tempk.append(mlek)
        self.params['k'] = np.array(tempk)
        self.params['n_samp'] = n_samp
        self.params['tot_obs'] = tot_obs
        self.var['k'] = np.array(tempk)

        return self

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
    k : float
        Aggregation parameter

    self.var keyword
    -----------------
    p : array of floats
        p parameter is 1 / n_samp
    k : array of floats
        Aggregation parameter. k is included if it is calculated in the fit
        method.

    '''
    
    @doc_inherit
    def __init__(self, **kwargs):
        self.params = kwargs
        self.min_supp = 0
        self.par_num = 2
        self.var = {}
    
    @doc_inherit
    def pmf(self, n):

        # TODO: Fix to work if n and N are one value
        #    if not (n <= N).all():
        #        raise Exception, "All values of n must be <= N."
        #    elif (a <= 0) or (a >= 1):
        #        raise Exception, "a must be between 0 and 1"

        # Get parameters
        n_samp, tot_obs, k = self.get_params(['n_samp', 'tot_obs', 'k'])
        n = expand_n(n, len(n_samp))
        
        # TODO: Additional checks?
        
        pmf = []
        self.var['p'] = []

        for tn_samp, ttot_obs, tk, tn in zip(n_samp, tot_obs, k, n):

            ln_L = lambda n_i,N,a,k: _ln_choose(n_i+k-1,n_i) + \
                _ln_choose(N-n_i+(k/a)-k-1,N-n_i) - _ln_choose(N +(k/a)-1,N)
            ta = 1 / tn_samp
            tpmf = ln_L(tn, ttot_obs, ta, tk) # Already log
            pmf.append(np.exp(tpmf))
            self.var['p'].append(ta)

        self.var['p'] = np.array(self.var['p'])

        return pmf
    
    def fit(self, data, upper_bnd=10):
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
        n_samp, tot_obs = self.get_params(['n_samp', 'tot_obs'])
        
        data = check_list_of_iterables(data)

        tempk = []

        for tdata, tn_samp, ttot_obs in zip(data, n_samp, tot_obs): 

            def nll_nb(k):
                self.params['tot_obs'] = ttot_obs
                self.params['n_samp'] = tn_samp
                self.params['k'] = k
                return -sum(np.log(self.pmf(tdata)[0]))
            
            mlek = scipy.optimize.brute(nll_nb, ((1e-10, upper_bnd),))
            tempk.append(mlek[0])

        self.params['k'] = np.array(tempk)
        self.var['k'] = np.array(tempk)

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

    self.var keywords
    -----------------
    p :  list of floats
        p parameter is equal to 1 / n_samp

    '''
    
    @doc_inherit
    def __init__(self, **kwargs):
        self.params = kwargs
        self.min_supp = 0
        self.par_num = 1
        self.var = {}
    
    @doc_inherit
    def pmf(self, n):

        # Get parameters
        n_samp, tot_obs = self.get_params(['n_samp', 'tot_obs'])
        n = expand_n(n, len(n_samp))
        
        # TODO: Additional checks?

        k = np.repeat(1, len(n_samp))
        pmf = nbd(tot_obs=tot_obs, n_samp=n_samp, k=k).pmf(n)
        self.var['p'] = 1 / n_samp
        return pmf
    
    @doc_inherit
    def cdf(self, n):

        # Get parameters
        n_samp, tot_obs = self.get_params(['n_samp', 'tot_obs'])
        n = expand_n(n, len(n_samp))
        
        # TODO: Additional checks?

        k = np.repeat(1, len(n_samp))
        cdf = nbd(tot_obs=tot_obs, n_samp=n_samp, k=k).cdf(n)
        self.var['p'] = 1 / n_samp
        return cdf
        
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

    self.var keywords
    ------------------
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
        self.var = {}
    
    @doc_inherit
    def pmf(self, n):
        
        # Get parameters
        n_samp, tot_obs = self.get_params(['n_samp', 'tot_obs'])
        n = expand_n(n, len(n_samp))
        
        # TODO: Additional checks?

        k = np.repeat(1, len(n_samp))
        tfnbd = fnbd(tot_obs=tot_obs, n_samp=n_samp, k=k)
        pmf = tfnbd.pmf(n)
        self.var=  tfnbd.var
        return pmf 
    
    @doc_inherit
    def cdf(self, n):
        
        # Get parameters
        n_samp, tot_obs = self.get_params(['n_samp', 'tot_obs'])
        n = expand_n(n, len(n_samp))
        
        # TODO: Additional checks?

        k = np.repeat(1, len(n_samp))
        tfnbd = fnbd(tot_obs=tot_obs, n_samp=n_samp, k=k).cdf(n)
        cdf = tfnbd.cdf(n)
        self.var = tfnbd.var
        return cdf

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

    self.var keywords
    -----------------
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
        self.var = {}
    
    @doc_inherit
    def pmf(self, n):

        # Get parameters
        n_samp, tot_obs = self.get_params(['n_samp', 'tot_obs'])
        n = expand_n(n, len(n_samp))
       
        # Define normalizing constant and pmf functions
        z_func = lambda x, ttot_obs: (1 - x ** (ttot_obs + 1)) / (1 - x)
        pmf_func = lambda z, x, tn: (1 / z) * (x ** tn)

        pmf = []
        self.var['x'] = []
        for tn_samp, ttot_obs, tn in zip(n_samp, tot_obs, n):
            ta = 1 / tn_samp

            #Compute probability directly to save time
            if ta == 0.5: 
                x = 1
                tpmf = np.repeat(1 / (1 + ttot_obs), len(tn))

            # All values zero except ttot_obs
            elif ta == 1:
                tpmf = np.zeros(len(tn))
                tpmf[np.where(tn == ttot_obs)[0]] = 1
                x = 0 

            elif ta < 0.5:
                try:
                    stop = 1 - 1e-10
                    # This is very brittle for some reason.  Changing the stop
                    # value can make this fail for strange reasons
                    x = scipy.optimize.brentq(l_solver, 0, .999999,
                            args=(ttot_obs, ta), disp=False)
                except:
                    try:
                        x = scipy.optimize.brentq(l_solver, 0, .95,
                                args=(ttot_obs, ta), disp=False)
                    except:
                        raise ValueError("No solution to " + 
                                "%s.pmf when tot_obs = " % 
                                (self.__class__.__name__) +
                                "%.2f, n_samp = %.10f and a = %.10f" % 
                                 (ttot_obs, tn_samp, ta))
                z = z_func(x, ttot_obs)
                tpmf = pmf_func(z, x, tn)
            else:
                try:
                    
                    x = scipy.optimize.brentq(l_solver, 0,
                            min((sys.float_info[0] * ta)**(1/float(ttot_obs)),
                                8), args=(ttot_obs, ta), disp=False,
                            xtol=1e-60, max_iter=200)

                except: 

                    try: # Allows it to pass, but optimizer starts rounding.
                         # Not Sure why it is doing this.
                        x = scipy.optimize.brentq(l_solver, 8.0, 50.0, \
                                   args=(ttot_obs, ta), disp=False, xtol=1e-60,
                                   max_iter=200)
                    except:

                        raise ValueError("No solution to " + 
                                "%s.pmf when tot_obs = " % 
                                (self.__class__.__name__) +
                                "%.2f, n_samp = %.10f and a = %.10f" % 
                                 (ttot_obs, tn_samp, ta))
                z = z_func(x, ttot_obs)
                tpmf = pmf_func(z, x, tn)

            pmf.append(tpmf)
            self.var['x'].append(x)

        return pmf

class mete_sar_iter(Curve):
    __doc__ = Curve.__doc__ + \
    '''
    Description
    -----------
    Computes values for the approximated METE generated SAR iteratively.

    Parameters
    ----------
    S : int
        Total number of species at the anchor area
    N : int
        Total number of individuals at the anchor area

    Notes
    -----
    This class uses method 1 in Harte (2011) to calculate the SAR.  It is much
    faster than the equivalent object gen_sar(logser_ut, tgeo) because of the
    approximations used. 
    
    '''
    
    @doc_inherit
    def __init__(self, **kwargs):
        self.params = kwargs

    def vals(*args, **kwargs):
        '''
        Vals is not implemented for mete_sar_iter
        '''
        raise NotImplementedError("'vals' method not implemented for class" + 
                                  " 'mete_sar_iter'")

    def iter_vals(self, a_list=None, upscale=0, downscale=0, non_iter=False,
                                                                     **kwargs):
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
        non_iter : bool
            If False, returns all iterations.  If True, only returns iterations
            that match a_list.

        Returns
        -------
        : 1D structured np.array
            The structured array has fields 'items' and 'area'. 'items' is
            equivalent to species.

        Notes
        -----
        With this method of the METE SAR, one cannot calculate the SAR at exact
        areas.  Rather this method iterates up and down by powers of 2.
        Therefore, the output of this function will contain all the SAR
        calculations in between ~min(a_list) ~max(a_list).


        '''
        #Get and check params
        S, N = self.get_params(['n_samp', 'tot_obs'])
        assert S < N, "S must be less than N"
        assert S > 1, "S must be greater than 1"
        assert N > 0, "S must be greater than 0"
        
        if a_list is not None:
            a_list = make_array(a_list)
        if not np.iterable(a_list) and a_list is not None:
            raise TypeError('a_list is not an array-like object')

        anch = 1
        if a_list != None:
            upscale, downscale = set_up_and_down(anch, a_list)
        
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

        if non_iter == False:
            return sar
        else:
            ind = np.zeros(len(sar), dtype=bool)
            for a in a_list:
                ind = np.bitwise_or(ind, sar['area'] == a)
            return sar[ind]

    @doc_inherit
    def univ_curve(self, num_iter=5, direction='down', **kwargs):

        return super(mete_sar_iter, self).univ_curve(num_iter=num_iter,
                  direction=direction, param='tot_obs', iterative=True, base=2)
        
    def fit(self, *args):
        '''
        This fit method fills the required parameters for a mete_sar_iter 
        object.

        Parameters
        ----------
        full_sad : array-like object
            The full_sad at the anchor area

        Notes
        -----
        full_sad is the first argument in the tuple args.  If there is a second
        argument, mete_sar_iter.fit() will ignore it.
        

        '''
        full_sad = args[0]
        self.params['n_samp'] = len(full_sad)
        self.params['tot_obs'] = sum(full_sad)

        return self

class powerlaw(Curve):
    __doc__ = Curve.__doc__ + \
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
    def vals(self, a_list, **kwargs):
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
            ('area', np.float)]. 'items' can be thought of as 'species'. 
    
        '''
        z, c = self.get_params(['z', 'c'])
        a_list = make_array(a_list)
        output_array = np.empty(len(a_list), dtype=[('items', np.float),
                                                     ('area', np.float)])
        output_array['area'] = a_list
        p_law = lambda x: c * (x ** z)
        output_array['items'] = p_law(a_list)
        return output_array
    
    def fit(self, *args):
        '''
        This fit method fills the required parameters for a powerlaw 
        object.

        Parameters
        ----------
        full_sad : array-like object
            The full_sad at the anchor area

        data :  tuple of array-like objects
            data contains two array-like objects.  The first is a list of area
            fractions (area / anchor area) and the second is a list of
            species/items numbers corresponding to those area fractions.

        Notes
        -----
        full_sad parameter should be first object in args and data should be
        second object in args.

        '''
        
        # Check and unpack args
        assert len(args) == 2, 'Expected two arguments'
        full_sad, data = args
        
        assert len(data) == 2, "data must contain two objects"
        assert len(data[0]) == len(data[1]), "Area and species number " + \
                                        "arrays must be of the same length"
        full_sad = make_array(full_sad)
        self.params['n_samp'] = len(full_sad)
        self.params['tot_obs'] = sum(full_sad)

        # Fit power law to regression
        reg = stats.linregress(np.log(data[0]), np.log(data[1]))
        self.params['z'] = reg[0]
        self.params['c'] = np.exp(reg[1])
        return self

class gen_sar(Curve):
    __doc__ = Curve.__doc__ + \
    '''
    A generic species area relationship (SAR/sar) function the utilizes the
    relationship between the sad and the ssad to generate the sar.  Can take
    any combination of sad and ssad. This function can also generate a endemic
    species area relationship (EAR)

    Parameters
    ----------
    n_samp : float
        Total number of species at the anchor area
    tot_obs : float
        Total number of individuals at the anchor area

    Notes
    -----
    plognorm and plognorm_lt are not supported by gen_sar. If one would like
    them to be supported, the full pmf for the sad must be calculated in the
    fit method.

    Examples
    --------
    import distributions as dist

    # Make an SAR with a Logseries SAD and Truncated Geometric SSAD.  The
    # community has 500 individuals and 14 species

    sar1 = dist.gen_sar(dist.logser(), dist.tgeo(), tot_obs=500, n_samp=14)

    # Number of species in half the base area and double the base area
    sar1.vals([.5, 2])

    # Make an SAR with Logseries and Truncated NBD
    sar2 = dist.gen_sar(dist.logser(), dist.tnbd(k=.2), tot_obs=500, n_samp=14)

    # Iterated the SAR 2 doublings from the base scale 
    sar2.iter_vals(upscale=2)


    '''

    def __init__(self, sad, ssad, **kwargs):
        '''

        Parameters
        ----------
        sad : a sad distribution object
            A Distribution object with minimum support equal to 1. pmf of sad
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

    def get_name(self):
        '''
        Returns the name of the general SAR which is a concatenation of the sad
        and ssad object

        Returns
        -------
        : str
            The name of the general SAR

        '''
        
        return self.sad.__class__.__name__ + '-' + self.ssad.__class__.__name__

    def iter_vals(self, a_list=None, upscale=0, downscale=0, non_iter=False, 
                                            base=2, use_rad=False, form='sar'):
        '''
        Calculates values in a_list by iteration.

        Parameters
        ----------
        a_list : array-like object
            List of area fractions at which to calculate the SAR
        upscale : int
            Number of iterations up from the anchor scale.  Each iteration 
            multiplies the previous area by base parameter. Only active if 
            a_list is None.
        downscale : int
            Number of iterations down from the anchor scale. Each iteration 
            multiplies the previous area by 1 / base. Only active if a_list 
            is None.
        non_iter : bool
            If False, returns all iterations.  If True, only returns iterations
            that match a_list.
        base : int
            Specifies the base of the logarithm.  In other words, whether you
            would like to iterate via double and half (base = 2), triple and
            third (base = 3), etc.
        use_rad : bool
            If False, uses the pmf of a distribution to calculate the species
            area relationship.  It True, uses a rank abundance distribution to
            calculate the species area relationship. If you chose to use the
            rank abundance distribution, please note that you cannot upscale
            with the rad past one iteration.
        form : string
            Default value is 'sar' which calculates the SAR given the
            parameters. You can also use 'ear' which calculates the EAR with
            the given parameters. 

        Returns
        -------
        : structured np.array
            A structured np.array with dtype=[('items', np.float),
            ('area', np.float)]. Items can be species.

        Notes
        -----
        At the moment, you cannot upscale the EAR
        '''

        S, N = self.get_params(['n_samp', 'tot_obs'])
        assert S < N, "S must be less than N"
        assert S > 1, "S must be greater than 1"
        assert N > 0, "S must be greater than 0"
        
        if a_list is not None:
            a_list = make_array(a_list)

        if not np.iterable(a_list) and a_list is not None:
            raise TypeError('a_list is not an array-like object')

        anch = 1
        if a_list != None:
            upscale, downscale = set_up_and_down(anch, a_list, base=base)

        if upscale == 0 and downscale == 0:
            return np.array((S, anch), dtype=[('items', np.float),
                                                ('area', np.float)])
        areas = _generate_areas_(anch, upscale, downscale, base=base)
        sar = np.empty(len(areas), dtype=[('items', np.float),
                                      ('area', np.float)])
        sar['area'] = areas
        def up_down_scale(areas, up_down):
            N_list = []; S_list = [] 
            
            if up_down == 'up':
                a = base #iterate up given base
            else:
                a = 1. / base # iterate down given base

            for i, da in enumerate(areas):
                if i == 0: # Base area calculation. Not always exactly S.
                           # For example, the broken_stick pmf does sum to
                           # exactly one so S does not quite come out right

                    self.params['tot_obs'] = N 
                    self.params['n_samp'] = S
                    S_list.append(self.vals([da], use_rad=use_rad, 
                                                        form=form)['items'][0])
                    N_list.append(N)

                else:
                    if np.isnan(S_list[i - 1]):
                        self.params['tot_obs'] = N
                        self.params['n_samp'] = S
                        if use_rad:
                            raise ValueError('Cannot calculate species number'
                             + ' at area = %s. Species value at previous' %
                             str(da) + ' iteration is NaN. Set parameter' + 
                             " 'use_rad' to False and try again.")
                        else:
                            raise ValueError('Cannot calculate species number'
                             + ' at area = %s. Species value at previous' %
                             str(da) + ' iteration is NaN. Use less iterations.')

                    self.params['tot_obs'] = N_list[i - 1]
                    self.params['n_samp'] = S_list[i - 1]
                    S_list.append(self.vals([a], use_rad=use_rad, 
                                                        form=form)['items'][0])

                    # Can't have less then one individual
                    if N * da < 1:
                        raise DownscaleError("Can't downscale %i iterations below"
                            % (downscale) + " %.2f individuals at the anchor scale"
                            % (N))
                    N_list.append(N * da)

            # Reset anchor values
            self.params['tot_obs'] = N
            self.params['n_samp'] = S
            self.sad.params['tot_obs'] = N
            self.sad.params['n_samp'] = S

            if up_down == 'down':
                return np.array(S_list)[::-1]
            else:
                return np.array(S_list)
        
        if upscale != 0:
            sar['items'][downscale:] = up_down_scale(areas[downscale:], 'up')

        if downscale != 0:
            sar['items'][:downscale + 1] = up_down_scale(areas[:downscale +
                                                              1][::-1], 'down')
        if non_iter == False or a_list == None:
            return sar
        else:
            ind = np.zeros(len(sar), dtype=bool)
            for a in a_list:
                ind = np.bitwise_or(ind, np.round(sar['area'], decimals=8) ==
                                                       np.round(a, decimals=8))
            return sar[ind]

    def vals(self, a_list, use_rad=False, form='sar'):
        '''

        Calculates sar value at each value in a_list

        Parameters
        ----------
        a_list : array-like object
            List of area fractions at which to calculate the SAR
        use_rad : bool
            If False, uses the sad pmf to calculate the SAR.  If True, uses the
            sad rank abundance distribution to calculate the SAR.
        form : string
            Default value is 'sar' which calculates the SAR given the
            parameters. You can also use 'ear' which calculates the EAR with
            the given parameters. 

        Returns
        -------
        : structured np.array
            A structured np.array with dtype=[('items', np.float),
            ('area', np.float)]. Items can be species.

        Notes
        -----
        Setting use_rad=True makes the calculations much quicker.  However, the
        rad relies on whole species numbers and therefore caution should be
        used when upscaling. Using the rad while upscaling will not give
        exactly the same value at upscaling with the pmf.  However, the SAR
        curves should have the same general pattern.

        At the moment you cannot upscale the EAR.
        '''
        
        # If sad is plognorm or plognorm_lt Throw and error for now
        nm = self.sad.__class__.__name__
        if nm == 'plognorm' or nm == 'plognorm_lt':
            raise ValueError("SAD %s not supported" % (nm))

        # Check that form is correct
        if not(form == 'sar' or form == 'ear'):
            raise ValueError("Parameter 'form' with value '%s' is not supported" %
                            (form))
    
        # Calculating sad in this method, not in fit.  More flexible this way.
        # However, this is a bit slower
        S, N = self.get_params(['n_samp', 'tot_obs'])
        self.sad.params['n_samp'] = S; self.sad.params['tot_obs'] = N

        # Calculate either rad or full pmf
        if use_rad:
            # If n_samp is fractional, need to round
            self.sad.params['n_samp'] = np.round(self.sad.params['n_samp'],
                                                                    decimals=0)

            rad = self.sad.rad()[0]
        else:
            sad = self.sad.pmf(np.arange(1, np.floor(N) + 1))[0]
        ssad = self.ssad
        sar = []

        a_list = make_array(a_list)
        for i, a in enumerate(a_list):
            
            #Setting ssad parameters differently if using rad 
            if use_rad:
                ssad.params['tot_obs'] = rad
            else:
                N_range = np.arange(1, len(sad) + 1)
                ssad.params['tot_obs'] = N_range

            # Upscale
            if a > 1:

                sad_params = deepcopy(self.sad.params)
                def eq(Sbig, abig, S):
                    # Setting distribution parameters for guess at upscale
                    # NOTE: You can't refit plognorm when you upscale. 
                    Nbig = np.round(abig * self.params['tot_obs'], decimals=0)
                    self.sad.params['tot_obs'] = Nbig

                    if use_rad:
                        self.sad.params['n_samp'] = np.round(Sbig, decimals=0)
                        ssad.params['tot_obs'] = self.sad.rad()[0]
                    else:
                        self.sad.params['n_samp'] = Sbig
                        ssad.params['tot_obs'] = np.arange(1, Nbig + 1)
                    ssad.params['n_samp'] = np.repeat(abig,
                                                   len(ssad.params['tot_obs']))
                    
                    # Probability of presence list
                    if form == 'sar':
                        p_pres_list = [1 - absnt[0] for absnt in ssad.pmf(0)]
                    elif form == 'ear':
                        p_pres_list = [fval[0] for fval in 
                                        ssad.pmf(zip(ssad.params['tot_obs']))]
                        
                    if use_rad:
                        return sum(np.array(p_pres_list)) - S
                    else:
                        sadbig = self.sad.pmf(np.arange(1, Nbig + 1))[0]
                        val = sum(Sbig * sadbig * np.array(p_pres_list)) - S
                        return val
                
                #Optimizing to find Sbig. If error set to nan
                try:
                    Sbig = scipy.optimize.brentq(eq, S, a * S, args=(a, S), disp=0)
                    sar.append(Sbig)
                except(ValueError):
                    print 'Could not calculate species number with values' +\
                           ' a = %s and S = %s' % (str(a), str(S))
                    sar.append(np.nan)

                self.sad.params = sad_params # Reset sad params

            elif a == 1:
                if use_rad:
                    sar.append(S)
                else:
                    p_pres_list = np.repeat(1, len(N_range))
                    sar.append(sum(S * sad * np.array(p_pres_list)))

            # Downscale
            else:
                ssad.params['n_samp'] = np.repeat(1 / a,
                                                   len(ssad.params['tot_obs']))
                if form == 'sar':
                    p_pres_list = [1 - absnt[0] for absnt in ssad.pmf(0)]
                elif form == 'ear':
                    p_pres_list = [fval[0] for fval in 
                                        ssad.pmf(zip(ssad.params['tot_obs']))]
                if use_rad:
                    sar.append(sum(np.array(p_pres_list)))
                else:
                    sar.append(sum(S * sad * np.array(p_pres_list)))

        return np.array(zip(sar, a_list), dtype=[('items', np.float), 
                                                  ('area', np.float)])
    
    def fit(self, *args):
        '''
        This fit method fills the required parameters for an gen_sar object.
        For the gen_sar object, the pmf is remade each time the vals method is
        called. A bit slower but more flexible.

        Parameters
        ----------
        full_sad : array-like object
            The full_sad at the anchor area (species/items counts)

        Notes
        -----
        full_sad is the first object in args

        '''
        full_sad = args[0]
        full_sad = make_array(full_sad)

        # Get N/tot_obs and S/n_samp from full sad
        self.params['n_samp'] = len(full_sad)
        self.params['tot_obs'] = sum(full_sad)
        
        self.sad.fit([full_sad])
        return self


#########################
## Energy Distributions##
#########################


class psi(Distribution):
    __doc__ = Distribution.__doc__ + \
    '''
    The individual energy distribution (IED) described in Harte (2011)

    Parameters
    ----------
    n_samp : int or iterable
        Total number of species / samples
    tot_obs: int or iterable
        Total number of individuals / observations
    E : int or iterable
        Total energy output of community

    self.var keywords 
    ------------------
    beta : list of floats
        The beta lagrange multiplier
    lambda_2 : list of floats
        The lambda2 lagrange multiplier

    Notes
    -----
    All other lagrange multipliers can be calculated from beta and lambda_2.

    '''

    @doc_inherit
    def __init__(self, **kwargs):

        self.params = kwargs
        self.par_num = 2        
        self.min_supp = 1
        self.var = {}

    @doc_inherit
    def pdf(self, e):
        #Get and check parameters
        n_samp, tot_obs, E = self.get_params(['n_samp', 'tot_obs', 'E'])
        e = expand_n(e, len(n_samp))

        start = 0.3
        stop = 2
        flmax = sys.float_info[0]

        pdf = []
        self.var['beta'] = []
        self.var['lambda_2'] = []

        for tn_samp, ttot_obs, tE, te in zip(n_samp, tot_obs, E, e):
            k = np.linspace(1, ttot_obs, num=ttot_obs)

            try:
                tx = scipy.optimize.brentq(beta_solver, start,
                               min((flmax/tn_samp)**(1/float(ttot_obs)), stop), 
                               args = (k, ttot_obs, tn_samp), disp=True)
            except(ValueError):
                raise ValueError("No solution to %s.pmf for tot_obs = %.2f"
                                 % (self.__class__.__name__, ttot_obs) + 
                                 " and n_samp = %.2f" % (tn_samp))

            # Set lagrange multipliers
            tbeta = -np.log(tx)
            tl2 = float(tn_samp) / (tE - ttot_obs) # Harte (2011) 7.26
            tl1 = tbeta - tl2
            tsigma = tl1 + (tE * tl2)
            
            norm = (float(tn_samp) / (tl2 * ttot_obs)) * (((np.exp(-tbeta) - \
                    np.exp(-tbeta*(ttot_obs + 1))) / (1 - np.exp(-tbeta))) - \
                    ((np.exp(-tsigma) - np.exp(-tsigma*(ttot_obs + 1))) / \
                    (1 - np.exp(-tsigma)))) #Harte (2011) 7.22

            #Notation from E.W.
            exp_neg_gamma = np.exp(-(tbeta + (te - 1) * tl2))

            #tpdf = tl2 * tbeta * (exp_neg_gamma / (1 - exp_neg_gamma)**2)
            
            tpdf = (float(tn_samp) / (ttot_obs * norm)) * exp_neg_gamma * (1
                        - ((ttot_obs + 1) * exp_neg_gamma ** ttot_obs) + (ttot_obs *
                        exp_neg_gamma ** (ttot_obs + 1))) / ((1 -
                        exp_neg_gamma) ** 2)

            #tpdf = (float(tn_samp) / (ttot_obs * norm)) * ((exp_neg_gamma  / \
            #        (1 - exp_neg_gamma)**2) - ((exp_neg_gamma**ttot_obs / \
            #        (1 - exp_neg_gamma)) * (ttot_obs + (exp_neg_gamma / \
            #        (1 - exp_neg_gamma)))))
                    # Harte (2011) 7.24
            pdf.append(tpdf)
            self.var['beta'].append(tbeta)
            self.var['lambda_2'].append(tl2)
        
        return pdf
    
    @doc_inherit
    def cdf(self, e):

        n_samp, tot_obs, E = self.get_params(['n_samp', 'tot_obs', 'E'])
        e = expand_n(e, len(n_samp))

        start = 0.3
        stop = 2
        flmax = sys.float_info[0]

        cdf = []

        self.var['beta'] = []
        self.var['lambda_2'] = []

        for tn_samp, ttot_obs, tE, te in zip(n_samp, tot_obs, E, e):
            k = np.linspace(1, ttot_obs, num=ttot_obs)

            try:
                tx = scipy.optimize.brentq(beta_solver, start,
                               min((flmax/tn_samp)**(1/float(ttot_obs)), stop), 
                               args = (k, ttot_obs, tn_samp), disp=True)
            except(ValueError):
                raise ValueError("No solution to %s.cdf for tot_obs = %.2f"
                                 % (self.__class__.__name__, ttot_obs) + 
                                 " and n_samp = %.2f" % (tn_samp))

            # Set lagrange multipliers
            tbeta = -np.log(tx)
            tl2 = float(tn_samp) / (tE - ttot_obs) # Harte (2011) 7.26
            tl1 = tbeta - tl2

            # Exact cdf equation. 
            eq2 = lambda x: tbeta * ((1 / (1 - np.exp(tl1 + (tl2 * x)))) - 
                                    (1 / (1 - np.exp(tl1 + tl2))))

            cdf.append(eq2(te))
            self.var['beta'].append(tbeta)
            self.var['lambda_2'].append(tl2)

        return cdf

    @doc_inherit
    def rad(self):
        
        n_samp, tot_obs, E = self.get_params(['n_samp', 'tot_obs', 'E'])

        start = 0.3
        stop = 2
        flmax = sys.float_info[0]

        n_arrays = [np.arange(1, i + 1) for i in tot_obs]
        
        # Define the predicted rad
        prad = lambda beta, r, tot_obs, l1, l2: (1 / l2) * np.log(((beta *\
                                   tot_obs) + r - 0.5) / (r - 0.5)) - (l1 / l2)
        rad = []
        for tn_samp, ttot_obs, tE, tn, in zip(n_samp, tot_obs, E, n_arrays):

            k = np.linspace(1, ttot_obs, num=ttot_obs)
            try:
                tx = scipy.optimize.brentq(beta_solver, start,
                                min((flmax/tn_samp)**(1/float(ttot_obs)), stop), 
                                args = (k, ttot_obs, tn_samp), disp=True)
            except(ValueError):
                raise ValueError("No solution to %s.rad for tot_obs = %.2f"
                                 % (self.__class__.__name__, ttot_obs) + 
                                 " and n_samp = %.2f" % (tn_samp))
            tbeta = -np.log(tx)
            tl2 = float(tn_samp) / (tE - ttot_obs) # Harte (2011) 7.26
            tl1 = tbeta - tl2

            trad = prad(tbeta, tn, ttot_obs, tl1, tl2)
            rad.append(trad)

        return rad
        

    def fit(self, data):
        '''
        Fit the community individual energy distribution data
        
        Parameters
        ----------
        data : list of tuples
            
            A list containing tuples of length two.  The first object in a
            tuple an iterable containing the community individual energy
            distribution.  The second object in a tuple is an iterable
            containing the empirical species abundance distribution. 

        '''

        # Unpack the list of tuples
        ied, sad = unpack(data)

        # Use base class fit
        super(psi, self).fit(sad)

        # Format and check energy data
        data_eng = check_list_of_iterables(ied)

        # Store energy data in self.params
        E = [np.sum(np.array(edata)) for edata in data_eng]
        self.params['E'] = E

        return self

class theta(Distribution):
    __doc__ = Distribution.__doc__ + \
    '''
    The species energy distribution (SED) as described by Harte (2011).

    Parameters
    ----------
    n_samp : int or iterable
        Total number of species / samples
    tot_obs: int or iterable
        Total number of individuals / observations
    E : int or iterable
        Total energy output of community
    n : int or iterable
        Number of individuals in a given species

    self.var keywords
    ------------------
    lambda_2 : list of floats
        The lambda2 lagrange multiplier

    '''

    @doc_inherit
    def __init__(self, **kwargs): 

        self.params = kwargs
        self.par_num = 2
        self.min_supp = 1
        self.var = {}
    
    @doc_inherit
    def pdf(self, e):

        n_samp, tot_obs, E, n = self.get_params(['n_samp', 'tot_obs', 'E','n'])
        e = expand_n(e, len(n_samp))
        
        # TODO: More checks?
        assert np.all(n <= tot_obs), 'n must be less than or equal to tot_obs'


        pdf = []
        self.var['lambda_2'] = []

        for tn_samp, ttot_obs, tE, tn, te in zip(n_samp, tot_obs, E, n, e):

            tl2 = float(tn_samp) / (tE - ttot_obs)
            tpdf = (tn * tl2 * np.exp(-tl2 * tn * te)) / (np.exp(-tl2 * tn)\
                                - np.exp(-tl2 * tn * tE)) #Harte (2011) 7.25

            pdf.append(tpdf)
            self.var['lambda_2'].append(tl2)
        
        return pdf

    @doc_inherit
    def cdf(self, e):

        n_samp, tot_obs, E, n = self.get_params(['n_samp', 'tot_obs', 'E','n'])
        e = expand_n(e, len(n_samp))
        
        # TODO: More checks?
        assert np.all(n <= tot_obs), 'n must be less than or equal to tot_obs'


        cdf = []
        self.var['lambda_2'] = []

        for tn_samp, ttot_obs, tE, tn, te in zip(n_samp, tot_obs, E, n, e):

            tl2 = float(tn_samp) / (tE - ttot_obs)

            # Exact cdf
            tcdf = -np.exp(tl2 * tn) * (np.exp(-tl2 * tn * te) - 
                                                        np.exp(-tl2 * tn))   

            cdf.append(tcdf)
            self.var['lambda_2'].append(tl2)
        
        return cdf

    @doc_inherit
    def rad(self, tol=.1):

        n_samp, tot_obs, E, n = self.get_params(['n_samp', 'tot_obs', 'E','n'])
        
        # TODO: More checks?
        assert np.all(n <= tot_obs), 'n must be less than or equal to tot_obs'

        n_arrays = [np.arange(1, i + 1) for i in n]

        prad = lambda r, n, l2, e_max : 1 + (1 / (l2 * n)) * \
                    np.log(1 / (((r - 0.5) / n) + np.exp(-l2 * n * (e_max))))

        rad = []
        terms= []
        for tn_samp, ttot_obs, tE, tn, tn_arr in zip(n_samp, tot_obs, E, n, 
                                                        n_arrays):
            
            tl2 = float(tn_samp) / (tE - ttot_obs)

            # Exact cdf
            cdf_eq = lambda es: -np.exp(tl2 * tn) * (np.exp(-tl2 * tn * es) - 
                                                        np.exp(-tl2 * tn))

            # Calculate e_max for correction term
            max_obs_cdf = 1 - (1 / (2 * tn))
            min_func = lambda e_max: max_obs_cdf - cdf_eq(e_max)
            e_max = scipy.optimize.brentq(min_func, 1, tE, disp=True)
            terms.append((tn, e_max, np.exp(-tl2 * tn * (e_max))))

            trad = prad(tn_arr, tn, tl2, e_max)
            rad.append(trad)

        return rad

    def fit(self, data):
        '''
        Fits empirical species energy distribution data

        Parameters
        ----------
        data : list of tuples

            A list of tuple where each tuple has length 3.  The first object in
            a tuple is an iterable containing the empirical species energy
            distribution.  The second object in a tuple is an iterable
            containing the community individual energy distribution.  The third
            object in a tuple is an iterable containing the empirical species
            abundance distribution.

        '''
        #TODO: Check format of data? 
        # Unpack the tuples
        sed, ied, sad = unpack(data)

        super(theta, self).fit(sad)

        # Check and set energy data
        data_eng = check_list_of_iterables(ied)
        E = [np.sum(np.array(edata)) for edata in data_eng]
        self.params['E'] = E
        
        # Check and set species abundance data
        n_data = check_list_of_iterables(sed)
        n = [len(np.array(ndata)) for ndata in n_data]
        self.params['n'] = n

        
        return self

# This distribution is a pain
class nu(Distribution):
    '''
    An energy distribution describing the distribution of average energy across
    all species in a community.  The average species energy distribution (ASED)
    
    Parameters
    ----------
    n_samp : int or iterable
        Total number of species / samples
    tot_obs: int or iterable
        Total number of individuals / observations
    E : int or iterable
        Total energy output of community

    self.var keywords
    -----------------
    lambda_2 : list of floats
        The lambda2 lagrange multiplier

    Notes
    -----
    This is a discrete distribution.
    '''

    @doc_inherit
    def __init__(self, **kwargs):
        self.params = kwargs
        self.par_num = 2
        self.min_supp = 1
        self.var = {}

    @doc_inherit
    def pmf(self, e):
        '''
        Notes
        -----
        The nu distribution is only defined at e values given by 
        e = 1 + (1 / (lambda2 * n)). While this function will return a pmf 
        value for all e greater than or equal to one, note that the pmf will
        only sum to one when provided with the proper support. lambda2 can be
        calculated by the equation: n_samp / (E - tot_obs) or S / (E - N)


        '''

        n_samp, tot_obs, E = self.get_params(['n_samp', 'tot_obs', 'E'])
        e = expand_n(e, len(n_samp))

        
        pmf = []
        self.var['lambda_2'] = []

        convert_e = lambda ep, l2: 1 / (l2 * (ep - 1))
        
        for tn_samp, ttot_obs, tE, te in zip(n_samp, tot_obs, E, e):

            # Set lagrange multipliers
            tl2 = float(tn_samp) / (tE - ttot_obs) # Harte (2011) 7.26
            e_max = 1 + (1 / tl2)
            e_min = 1 + (1 / (ttot_obs * tl2))
            
            tpmf = np.empty(len(te), dtype=float)
            tns = np.ceil(convert_e(te, tl2))
            
            # Parse values that aren't in range and set to zero
            ind_tot = np.arange(len(tpmf))
            ind_less = np.where(te >= e_min)[0]
            ind_more = np.where(te <= e_max)[0]
            ind_include = np.intersect1d(ind_more, ind_less)
            ind_exclude = np.array(list(set(ind_tot) - set(ind_include)))
            if len(ind_exclude) != 0:
                tpmf[ind_exclude] = 0
            
            if len(ind_include) != 0:
                tpmf[ind_include] = logser_ut(tot_obs=ttot_obs,
                        n_samp=tn_samp).pmf(tns[ind_include])[0]

            pmf.append(tpmf)
            self.var['lambda_2'].append(tl2)

        return pmf
    
    @doc_inherit
    def cdf(self, e):
        
        n_samp, tot_obs, E = self.get_params(['n_samp', 'tot_obs', 'E'])
        e = expand_n(e, len(n_samp))

        
        cdf = []
        self.var['beta'] = []
        self.var['lambda_2'] = []

        convert_n = lambda n, l2: 1 + (1 / (n *  l2))

        for tn_samp, ttot_obs, tE, te in zip(n_samp, tot_obs, E, e):

            tl2 = float(tn_samp) / (tE - ttot_obs) # Harte (2011) 7.26

            # Set all e so you can sum
            all_e = convert_n(np.arange(1, ttot_obs + 1), tl2)[::-1]

            pmf_for_all_e = nu(tot_obs=ttot_obs, n_samp=tn_samp,
                    E=tE).pmf(all_e)[0]
            cum_sum = np.cumsum(pmf_for_all_e)

            tcdf = np.array([cum_sum[np.sum(e_val >= all_e) - 1] if sum(e_val
                >= all_e) - 1 != -1 else 0 for e_val in te])

            cdf.append(tcdf)
            self.var['lambda_2'].append(tl2)

        return cdf
        
    def rad(self):
        '''
        This rad uses the observed cdf for a given nu distribution and the
        predicted cdf to calculate the rank energy distribution.  

        Returns
        -------
        : list
            A list of rank energy distributions 

        '''
    
        n_samp, tot_obs, E = self.get_params(['n_samp', 'tot_obs', 'E'])
        rad = []
        
        convert_n = lambda n, l2: 1 + (1 / (n *  l2))

        for tn_samp, ttot_obs, tE in zip(n_samp, tot_obs, E):

            # Set temp params
            self.params['n_samp'] = tn_samp
            self.params['ttot_obs'] = ttot_obs
            self.params['E'] = tE

            tl2 = float(tn_samp) / (tE - ttot_obs) # Harte (2011) 7.26
            all_e = convert_n(np.arange(1, ttot_obs + 1), tl2)[::-1]
            tpmf = self.pmf(all_e)[0]
            tcdf = np.cumsum(tpmf)

            # Observed cdf. Not quite true if some energies overlap
            obs_cdf = np.arange(1 / (2 * (tn_samp)), 1, 1/tn_samp)

            trad = [all_e[sum(oc >= tcdf) - 1] if sum(oc >= tcdf) - 1 != -1
                    else all_e[0] for oc in obs_cdf]

            rad.append(trad)

        self.params['n_samp'] = n_samp
        self.params['ttot_obs'] = tot_obs
        self.params['E'] = E

        return rad


    def fit(self, data):
        '''
        Fit the average species energy distribution to data
        
        Parameters
        ----------
        data : list of tuples
            
            A list containing tuples of length two or a list containing tuples
            of length three.  If the tuples are of length two, the first object
            in a tuple is an iterable containing the community individual energy
            distribution.  The second object in a tuple is an iterable
            containing the empirical species abundance distribution. If the
            tuples are of length three, the first object in the tuple is an
            iterable containing the average energy distribution. The second object
            in a tuple an iterable containing the community individual energy
            distribution.  The third object in a tuple is an iterable
            containing the empirical species abundance distribution.  

        '''

        # Unpack the list of tuples
        # Can either take 
        if len(data[0]) == 2:
            ied, sad = unpack(data)
        elif len(data[0]) == 3:
            ased, ied, sad = unpack(data)

        # Use base class fit
        super(nu, self).fit(sad)

        # Format and check energy data
        data_eng = check_list_of_iterables(ied)

        # Store energy data in self.params
        E = [np.sum(np.array(edata)) for edata in data_eng]
        self.params['E'] = E

        return self

class omega(Distribution):
    """
    This distribution is the distribution of total energy within a species
    across all species. The means of this distribution is E / S.

    Parameters
    ----------
    n_samp : int or iterable
        Total number of species / samples
    tot_obs: int or iterable
        Total number of individuals / observations
    E : int or iterable
        Total energy output of community

    self.var keywords
    -----------------
    lambda_2 : list of floats
        The lambda2 lagrange multiplier
    emaxmin : list fo tuples
        Each tuple contains the max total energy and min total energy for the
        given state variables.  

    Notes
    -----
    This is a discrete distribution.

    
    """
    
    def pmf(self, e):
        '''
        Notes
        -----
        The omega distribution is only defined at e values given by 
        e = n + (1 / lambda2). While this function will return a pmf 
        value for all e greater than or equal to one, note that the pmf will
        only sum to one when provided with the proper support. lambda2 can be
        calculated by the equation: n_samp / (E - tot_obs) or S / (E - N)


        '''

        n_samp, tot_obs, E = self.get_params(['n_samp', 'tot_obs', 'E'])
        e = expand_n(e, len(n_samp))
        
        pmf = []
        self.var['lambda_2'] = []
        self.var['emaxmin'] = []

        convert_e = lambda ep, l2: ep - (1 / l2)
        
        for tn_samp, ttot_obs, tE, te in zip(n_samp, tot_obs, E, e):

            # Set lagrange multipliers
            tl2 = float(tn_samp) / (tE - ttot_obs) # Harte (2011) 7.26
            e_max = ttot_obs + (1 / tl2)
            e_min = 1 + (1 / tl2)
            
            tpmf = np.empty(len(te), dtype=float)
            tns = convert_e(te, tl2)
            
            # Parse values that aren't in range and set to zero
            ind_tot = np.arange(len(tpmf))
            ind_less = np.where(te >= e_min)[0]
            ind_more = np.where(te <= e_max)[0]
            ind_include = np.intersect1d(ind_more, ind_less)
            ind_exclude = np.array(list(set(ind_tot) - set(ind_include)))
            if len(ind_exclude) != 0:
                tpmf[ind_exclude] = 0
            
            if len(ind_include) != 0:
                tpmf[ind_include] = logser_ut(tot_obs=ttot_obs,
                        n_samp=tn_samp).pmf(tns[ind_include])[0]

            pmf.append(tpmf)
            self.var['lambda_2'].append(tl2)
            self.var['emaxmin'].append((e_max, e_min))

        return pmf

    @doc_inherit
    def cdf(self, e):
        
        n_samp, tot_obs, E = self.get_params(['n_samp', 'tot_obs', 'E'])
        e = expand_n(e, len(n_samp))

        
        cdf = []
        self.var['lambda_2'] = []

        convert_n = lambda n, l2: n + (1 / l2)

        for tn_samp, ttot_obs, tE, te in zip(n_samp, tot_obs, E, e):

            tl2 = float(tn_samp) / (tE - ttot_obs) # Harte (2011) 7.26

            # Set all e so you can sum
            all_e = convert_n(np.arange(1, ttot_obs + 1), tl2)

            pmf_for_all_e = omega(tot_obs=ttot_obs, n_samp=tn_samp,
                    E=tE).pmf(all_e)[0]
            cum_sum = np.cumsum(pmf_for_all_e)

            tcdf = np.array([cum_sum[np.sum(e_val >= all_e) - 1] if sum(e_val
                >= all_e) - 1 != -1 else 0 for e_val in te])

            cdf.append(tcdf)
            self.var['lambda_2'].append(tl2)

        return cdf
        
    def rad(self):
        '''
        This rad uses the observed cdf for a given omega distribution and the
        predicted cdf to calculate the rank energy distribution.  

        Parameter
        ----------
        tol : float
            Precision interval.  The integral of omega is approximated at the
            interval tol.  Smaller intervals can be more precise, but a tol
            between 0.1 and 0.5 is more effecient and the results are changed
            only marginally.

        Returns
        -------
        : list
            A list of rank energy distributions 

        '''
    
        n_samp, tot_obs, E = self.get_params(['n_samp', 'tot_obs', 'E'])
        rad = []
        
        convert_n = lambda n, l2: n + (1 / l2)

        for tn_samp, ttot_obs, tE in zip(n_samp, tot_obs, E):

            # Set temp params
            self.params['n_samp'] = tn_samp
            self.params['ttot_obs'] = ttot_obs
            self.params['E'] = tE

            tl2 = float(tn_samp) / (tE - ttot_obs) # Harte (2011) 7.26
            all_e = convert_n(np.arange(1, ttot_obs + 1), tl2)
            tpmf = self.pmf(all_e)[0]
            tcdf = np.cumsum(tpmf)

            # Observed cdf. Not quite true if some energies overlap
            obs_cdf = np.arange(1 / (2 * (tn_samp)), 1, 1/tn_samp)

            trad = [all_e[sum(oc >= tcdf) - 1] if sum(oc >= tcdf) - 1 != -1
                    else all_e[0] for oc in obs_cdf]

            rad.append(trad)

        self.params['n_samp'] = n_samp
        self.params['ttot_obs'] = tot_obs
        self.params['E'] = E

        return rad


    def fit(self, data):
        '''
        Fit the average species energy distribution to data
        
        Parameters
        ----------
        data : list of tuples
            
            A list containing tuples of length two or a list containing tuples
            of length three.  If the tuples are of length two, the first object
            in a tuple is an iterable containing the community individual
            energy distribution.  The second object in a tuple is an iterable
            containing the empirical species abundance distribution. If the
            tuples are of length three, the first object in the tuple is an
            iterable containing the total species energy distribution. The
            second object in a tuple an iterable containing the community
            individual energy distribution.  The third object in a tuple is an
            iterable containing the empirical species abundance distribution.  

        '''

        # Unpack the list of tuples
        # Can either take 
        if len(data[0]) == 2:
            ied, sad = unpack(data)
        elif len(data[0]) == 3:
            tsed, ied, sad = unpack(data)

        # Use base class fit
        super(nu, self).fit(sad)

        # Format and check energy data
        data_eng = check_list_of_iterables(ied)

        # Store energy data in self.params
        E = [np.sum(np.array(edata)) for edata in data_eng]
        self.params['E'] = E

        return self

def nu_pmf_eq(es, beta, l2, s):
    '''Nu pmf
    
    Parameters
    ----------
    es : float or np.array
        Energy values at which to calculate pmf
    beta, l2 : float
        Lagrange multipliers
    s : float
        Species number
    
    Return
    ------
    : array or float
    '''
        
    # Nu pmf equation
    return (1 / np.log(s / beta)) * (np.exp(-beta / (l2 * (es - 1)))) / \
                                                                    (es - 1)

def l_solver(x, N, a):
    """
    Used with a solver to get the langrange multiplier for a pi distribution

    Parameters
    ----------
    x : float
        Lagrange multiplier x = e**-lambda
    N : float
        total balls (individuals) in urn (species)
    a : float
        area fraction. 1 / n_samp or 1 / urn_number

    Returns
    -------
    : float
    

    """
    return ((x / (1 - x)) - (((N + 1) * x ** (N + 1)) / \
                            (1 - x ** (N + 1)))) - (N * a)

def beta_solver(x, k, tot_obs, n_samp):
    """ Used with a solver to get the beta lagrange multiplier in the METE
    distributions.  With a solver, this function
    returns x and beta = -np.log(x)

    Parameters
    ----------
    x : float
        Lagrange multiplier x = e**-beta
    k : np.array
        np.arange(1, tot_obs + 1)
    tot_obs : float
        The total number of individuals observed (N in METE, see Harte 2011)
    n_samp : float
        The total number of species observed (S in METE, see Harte 2011)

    Returns
    -------
    : float
    """

    # Beta Solver
    return sum(x ** k / float(tot_obs) * n_samp) -  sum((x ** k) / k)


def make_array(n, dtype=None):
    '''Cast n as iterable array. If dtype not none this will be the dtype of
    the array.  Otherwise it lets python choose.  Must be a valid dtype or an
    error will be thrown'''
    if np.iterable(n):
        if dtype==None:
            return np.array(n)
        else:
            return np.array(n, dtype=dtype)
    else:
        if dtype==None:
            return np.array([n])
        else:
            return np.array([n], dtype=dtype)


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

def check_list_of_iterables(data):
    '''
    Checks if the given object is a list of iterables.  If so, returns a
    list of arrays.  Else, a TypeError is raised.
    
    Parameters
    ----------
    data : object
        The object to be tested.

    Returns
    -------
    : list
        If error is not thrown, returns a list of arrays

    '''
    
    # Check that data is a list of iterables
    if type(data) != type([]):
        raise TypeError('Data must be a list of iterables')
    if not np.all(np.array([np.iterable(dt) for dt in data])):
        raise TypeError('Objects in data must be iterable')

    # Make a list of arrays
    return [np.array(data) for data in data]

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
    min_supp : int
        The minimum support of the distribution. Often either 1 or 0.

    Returns
    -------
    S_abunds : ndarray
        1D array of predicted abundance for each species

    Notes
    -----
    Function actually implements (philosophically) a step quantile function.

    '''
    
    pmf = pmf / np.sum(pmf)  # Ensure distribution is normalized

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
    
    return counts # / (sum(counts) / len(pmf))


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
    gammaln = spec.gammaln
    return gammaln(n + 1) - (gammaln(k + 1) + gammaln(n - k + 1))

def set_up_and_down(anch, a_list, base=2):
    '''
    Sets the number of upscales and downscales given an a_list.
    By setting the base parameter, you can specify whether you would like to
    double and half (base = 2), triple an third (base=3), etc.
    
    Parameters
    ----------
    anch : float
        The anchor area
    a_list : array-like object
        List of areas
    base : int
        Base of the logarithm
    
    Returns
    -------
    : tuple
        Number of upscale and downscales required to cover the whole range
        specified in a_list.

    '''

    mint = np.min(a_list)
    maxt = np.max(a_list)
    if mint == anch and maxt == anch: 
        upscale = 0; downscale = 0
    elif (mint > anch or mint == anch) and maxt > anch:
        upscale = np.int(np.ceil(np.log(maxt / anch) / np.log(base)))
        downscale = 0
    elif mint < anch and (maxt < anch or maxt == anch):
        downscale = np.int(np.ceil(np.abs(np.log(mint / anch) / np.log(base)))) 
        upscale = 0
    elif mint < anch and maxt > anch:
        upscale = np.int(np.ceil(np.log(maxt / anch) / np.log(base)))
        downscale = np.int(np.ceil(np.abs(np.log(mint / anch) / np.log(base))))
    
    return upscale, downscale

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
                raise DownscaleError('Cannot downscale %i iterations from ' 
                                     % (len(down_areas) - 1) +\
                                     'anchor scale. One or less individuals' +\
                                     ' per cell.')
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

def _generate_areas_(anchor_area, upscale, downscale, base=2):
    '''
    Utility function that makes the area list
    
    '''

    areas = np.empty(upscale + downscale + 1)
    areas[downscale] = anchor_area
    for i in range(downscale)[::-1]:
        areas[i] = areas[i + 1] / base
    for i in xrange(downscale + 1, len(areas)):
        areas[i] = areas[i - 1] * base

    return areas

def unpack(zipped_data):
    '''
    Unpacks zipped data

    '''

    unzipped_data = zip(*zipped_data)
    unzipped_data = [list(tup) for tup in unzipped_data]
    return tuple(unzipped_data)
