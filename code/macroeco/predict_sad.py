#!/usr/bin/python


#A script to predict sad calculations
#TODO: Make appropriate python header
'''
Calculate pmf and likelihood of spatial-abundance distributions.

All distributions have an argument summary, which if False returns the entire 
pmf for the inputted values of n, and if true returns the summed negative 
log-likelihood of the inputted values (useful for likelihood ratio tests or 
AIC).

Distributions
-------------
- 'logser_pmf' -- Fisher's log series (Fisher et al. 1943)
- 'plognorm_pmf' -- Poisson lognormal (Bulmer 1974)
- 'mete_logser_pmf' -- METE log series (Harte 2011)
- 'mete_logser_approx_pmf' -- METE log series using approximation (Harte 2011)


References
----------

Bulmer, M. G. 1974. On fitting the poisson lognormal distribution to species
abundance data. Biometrics, 30:101-110.

Fisher, R. A., Corbet, A. S., and C. B. Williams. 1943. The relation between
The number of species and the number of individuals in a random sample
of an animal populatin. Journal of Animal Ecology, 12:42-58.

Harte, J. 2011. Maximum Entropy and Ecology: A Theory of Abundance,
Distribution, and Energetics. Oxford University Press.

Hubbell, S. P. 2001. The unified theory of biodiversity and biogeography. 
Monographs in Population Biology, 32,1:375.

'''

import numpy as np
import scipy as sp
import scipy.stats as stats
import scipy.optimize 
import scipy.special
import math as m
import exceptions
import scipy.integrate as integrate



class RootError(Exception):
    '''Error if no root or multiple roots exist for the equation generated
    for specified values of S and N'''

    def __init__(self, value=None):
        Exception.__init__(self)
        self.value = value
    def __str__(self):
        return '%s' % self.value

def lgser_pmf(S, N, summary=False):
    '''
    lgser_pmf(S, N, summary=False)
    
    Fisher's log series pmf (Fisher et al. 1943, Hubbel 2001)

    Parameters
    ----------
    S : int
        Total number of species in landscape
    N : int
        Total number of indviduals in landscape

    Returns
    -------
    : ndarray (1D)
        If summary is False, returns array with pmf. If summary is True,
        returns the summed log likelihood of all values in n.

    Notes
    -----
    Multiplying the pmf by S yields the predicted number of species
    with a given abundance

    '''
    assert S < N, "S must be less than N"
    assert S > 1, "S must be greater than 1"
    assert N > 0, "S must be greater than 0"
    start = -2
    stop = 1 - 1e-10
    
    eq = lambda x,S,N: (((N/x) - N) * (-(m.log(1 - x)))) - S
    
    x = scipy.optimize.brentq(eq, start, stop, args=(S,N), disp=True)
    k = np.linspace(1, N, num=N)
    pmf = stats.logser.pmf(k,x)

    #NOTE: If N is very large the sum will be nearly infinite...
    if summary: return -sum(np.log(pmf))
    else:       return pmf

def plognorm_pmf(abundances, summary=False):
    '''
    plognorm_pmf(abundance, summary=False)

    Poisson log-normal pmf (Bulmer 1974)

    Parameters
    ----------
    abundances : A list, np.array, or tuple of the abundance fo each
                 species in the plot.  len(abundance) should equal 
                 the total number of species in the plot and sum(abundances)
                 should equal the total number of individuals in the plot

     Returns
    -------
    : ndarray (1D)
        If summary is False, returns array with pmf. If summary is True,
        returns the summed log likelihood of all values in n.

    Notes
    -----
    Multiplying pmf by len(abundances) gives predicted number of species
    at each n.

    '''
    assert type(abundances) == list or type(abundances) == tuple \
           or type(abundances) == np.ndarray, "Invalid parameter type"
    assert len(np.unique(np.array(abundances))) > 1, "Error: variance is zero"

    abundances = np.array(abundances)
    N = int(np.sum(abundances))
    S = int(len(abundances)) #Might not need S
    log_abund = np.log(abundances)
    mean = np.mean(log_abund)
    var = np.var(log_abund, ddof=1)
    sd = var**0.5

    eq = lambda t, x: np.exp(t * x - np.exp(t) - 0.5*((t - mean) / sd)**2)
    pmf = np.empty(N)
    for i in xrange(1, N + 1):
        if i <= 170:
            integral = integrate.quad(eq, -np.inf, np.inf, args=(i))[0]
            norm = np.exp((-0.5 * m.log(2 * m.pi * var) - m.lgamma(i + 1)))
            pmf[i - 1] = norm * integral
        else:
            z = (m.log(i) - mean) / sd
            pmf[i - 1] = (1 + (z**2 + m.log(i) - mean - 1) / (2 * i * sd**2)) *\
                         np.exp(-0.5 * z**2) / (m.sqrt(2 * m.pi) * sd * i)   
        
    
    if summary: return -sum(np.log(pmf))
    else:       return pmf

        
def mete_lgsr_pmf(S, N, summary=False):
    '''
    mete_lgsr_pmf(S, N, summary=False)

    Truncated log series pmf (Harte 2011)

    Parameters:
    -----------
    S : int
        Total number of species in landscape
    N : int
        Total number of individuals in landscape
    summary: bool (optional)
        see Returns
    
    Returns
    -------
    : ndarray (1D)
        If summary is False, returns array with pmf.  If summary is True,
        returns the summed log likelihood of all values in n.

    Notes
    -----
    This function uses the truncated log series as described in Harte 2011
    eq (7.32).  The equation used in this function to solve for the Lagrange 
    multiplier is equation (7.27) as described in Harte 2011. 
    
    Also note that realistic values of x where x = e^-(beta) (see Harte 2011) are
    in the range (1/e, 1) (exclusive). Therefore, the start and stop parameters 
    for the brentq procedure are close to these values. However, x can
    occasionally be greater than one so the stop value of the brentq optimizer
    is 2.

    '''
    assert S < N, "S must be less than N"
    assert S > 1, "S must be greater than 1"
    assert N > 0, "S must be greater than 0"
    start = 0.3
    stop = 2
    
    #NOTE:  The bound is not N + S - 1
    n = np.linspace(1, N, num=N)

    eq = lambda x, S, N: (sum(x ** n) / sum((x ** n) / n)) - (float(N)/S)
    x = scipy.optimize.brentq(eq, start, stop, args=(S,N), disp=True)
    #Taken from Ethan White's trun_logser_pmf
    nvals = np.linspace(1, N, num=N)
    norm = sum(x ** nvals / nvals)        
    pmf = (x ** nvals/nvals) / norm

    if summary: return -sum(np.log(pmf))
    else:       return pmf


def mete_logsr_approx_pmf(S, N, summary=False, root=2):
    '''
    tflogsApprox(S, N, start=0.30, stop=0.999999, summary=False, root=2)

    Truncated log series using approximation (7.30) and (7.32) in Harte 2011

    Parameters
    ----------
    S : int
        Total number of species in landscape
    N : int
        Total number of individuals in landscape 
    summary: bool (optional)
        see Returns
    root: int (optional)
        1 or 2.  Specifies which root to use for pmf calculations
    
    Returns
    -------
    : ndarray (1D)
        If summary is False, returns array with pmf.  If summary is True,
        returns the summed log likelihood of all values in n.
    
    Notes:
    ------
    This function uses the truncated log series as described in Harte 2011
    eq (7.32).  The equation used in this function to solve for the Lagrange 
    multiplier is equation (7.30) as described in Harte 2011.     
       
    Also note that realistic values of x where x = e^-(beta) (see Harte 2011) are
    in the range (1/e, 1) (exclusive). Therefore, the start and stop parameters
    for the brentq optimizer have been chosen near these values

    
    '''

    #NOTE:  Ethan White has a way to do this as well, but his does not
    #take into account the fact that there can be multiple roots. Based on the
    # way that he implements this the approx method that is not actually
    # a major problem.  However, he does not account for the fact that there are
    # S and N values such that there is no solution to the constraint.  This 
    # method does. 
    assert S < N, "S must be less than N"
    assert S > 1, "S must be greater than 1"
    assert N > 0, "S must be greater than 0"
    start = 0.3
    stop = 1 - 1e-10
    
    eq = lambda x, S, N: ((-m.log(x))*(m.log(-1/(m.log(x))))) - (float(S)/N)
    try:
        x = scipy.optimize.brentq(eq, start, stop, args = (S, N), disp=True)
    except ValueError:
        values = np.linspace(start, stop, num=1000)
        solns = []
        for i in values:
            solns.append(eq(i, S, N))
        ymax = np.max(np.array(solns))
        xmax = values[solns.index(ymax)]
        if ymax > 0:
            print "Warning: More than one solution."
            if root == 1:
                x = scipy.optimize.brentq(eq, start, xmax, args=(S,N), disp=True)
            if root == 2:
                x = scipy.optimize.brentq(eq, xmax, stop, args=(S,N), disp=True)
                       
        if ymax < 0:
            raise RootError('No solution to constraint equation with given ' + \
                            'values of S and N') 
    k = np.linspace(1, N - S + 1, num = N - S + 1)
    g = -1/m.log(x)
    pmf = (1/m.log(g)) * ((x**k)/k) 
    
    if summary: return -sum(np.log(pmf))
    else:       return pmf









    
    


