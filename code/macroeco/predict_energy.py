#!/usr/bin/python

'''
Calculate pmf and likilihood of METE energy distributions.

All distributions have an argument summary, which if False returns the entire 
pmf for the inputted values of n, and if true returns the summed negative 
log-likelihood of the inputted values (useful for likelihood ratio tests or 
AIC).

Distributions
-------------
- 'mete_energy_theta_pdf -- The METE intra-specific energy distribution 
                            (Harte 2011)
- 'mete_energy_psi_pdf' -- The METE community energy distribution

- 'mete_energy_nu_pdf/cdf' -- The METE mean energy per species distribution


References
----------

Harte, J. 2011. Maximum Entropy and Ecology: A Theory of Abundance, Distribution,
and Energetics. Oxford University Press. 

'''
from __future__ import division
import numpy as np
import scipy.optimize 
import scipy.special
import scipy.integrate as integrate
import math as m
import sys

def mete_energy_theta_pdf(S, N, E, n, summary=False):
    '''
    mete_energy_theta_pdf(S, N, E, n, summary=False)

    METE intra-specific energy distribution (Harte 2011)

    Parameters
    ----------
    S : int
        Total number of species in landscape
    N : int
        Total number of individuals in landscape
    E : float
        Total energy ouput of community
    n : int
        Number of individuals in species of interest

    Returns
    -------
    : ndarray (1D)
        If summary is False, returns array with pmf. If summary is True,
        returns the summed log likelihood of all values in n.

    Notes
    -----

    '''
    assert S < N, "S must be less than N"
    assert S > 1, "S must be greater than 1"
    assert N > 0, "N must be greater than 0"
    assert N < E, "N must be less than E"
    assert n < N, "n must be less than N"

    lambda2 = float(S) / (E - N) #Harte (2011) 7.26
    
    epi = np.linspace(1, E, num=E*10) #Ten is arbitrary
    pdf = (n * lambda2 * np.exp(-lambda2 * n * epi)) / (np.exp(-lambda2 * n)\
          - np.exp(-lambda2 * n * E)) #Harte (2011) 7.25

    if summary: return -sum(np.log(pdf))
    else:       return pdf

def mete_energy_theta_cdf():
    pass

def mete_energy_nu_pdf(S, N, E, summary=False):
    '''

    mete_energy_nu_pdf(S, N, E, nmin, nmax summary=False)

    Parameters
    ----------
    S : int
        Total number of species in landscape
    N : int
        Total number of individuals in landscape
    E : float
        Total energy ouput of community

    Returns
    -------
    : ndarray (1D)
        If summary is False, returns array with pmf. If summary is True,
        returns the summed log likelihood of all values in n.
    
    Notes
    -----
    This function uses the nu energy equation found in in Table 7.3 of 
    Harte (2011). Nu is defined as the 'distribution of metabolic rates
    averaged over individuals within species' (Harte 2011). 

    '''

    assert S < N, "S must be less than N"
    assert S > 1, "S must be greater than 1"
    assert N > 0, "N must be greater than 0"
    assert N < E, "N must be less than E"

    
    #Start and stop for root finders
    start = 0.3
    stop = 2

    #Getting beta
    k = np.linspace(1, N, num=N)
    eq = lambda x: sum(x ** k / float(N) * S) -  sum((x ** k) / k)
    x = scipy.optimize.brentq(eq, start, min((sys.float_info[0] / S)**(1/float(N)),\
                              stop), disp=True)
    beta = -m.log(x)
    lambda2 = float(S) / (E - N) #Harte (2011) 7.26

    eta = np.linspace(1 + 1e-10, E, num=1000)
    pdf = (1 / np.log(1 / beta)) * ((np.exp(-(beta / (lambda2 * (eta - 1)))) / \
                                    (eta - 1)))
    
    if summary: return -sum(np.log(pdf))
    else:       return pdf

def mete_energy_nu_cdf(S, N, E, num_samples=100):
    '''
    mete_energy_nu_cdf(S, N, E, emin, emax, num_samples=100)

    CDF for mete_energy_nu_pdf

    Parameters
    ----------
    S : int
        Total number of species in landscape
    N : int
        Total number of individuals in landscape
    E : float
        Total energy ouput of community
    emin : float
        The minumum average energy output of a species
    emax : float
        The maximum average energy output of a species
    num_samples : int (default)
        The number of integrals computed within the interval [emin,emax]

    Returns
    -------
    : structured ndarray (1D)
        Structured array contains fields 'cdf' and 'energy'. The field 'cdf' represents
        the cdf value at a given energy level (Prob (emin <= energy <= upppere)

    Notes:
    ------
    This function uses the equation in Harte (2011) Table 7.3 to calculate
    emin and emax.  Note that because this equation is an approximation
    the integral of the pdf from emin to max will often be less than one.

    '''

    #Start and stop for root finders
    assert S < N, "S must be less than N"
    assert S > 1, "S must be greater than 1"
    assert N > 0, "N must be greater than 0"
    assert N < E, "N must be less than E"

    start = 0.3
    stop = 2

    #Getting beta
    k = np.linspace(1, N, num=N)
    eq = lambda x: sum(x ** k / float(N) * S) -  sum((x ** k) / k)
    x = scipy.optimize.brentq(eq, start, min((sys.float_info[0] / S)**(1/float(N)),\
                              stop), disp=True)
    beta = -m.log(x)
    lambda2 = float(S) / (E - N) #Harte (2011) 7.26

    pdf = lambda eta: (1 / np.log(1 / beta)) * ((np.exp(-(beta / (lambda2 * (eta - 1)))) / \
                                    (eta - 1)))
    cdf = []
    emin = 1 + (1 / (1 * lambda2))#Harte (2011), Table 7.3
    emax = 1 + (1 / (N * lambda2))
    for eta in np.linspace(emin, emax, num=num_samples):
        cdf.append((integrate.quad(pdf, emin, eta)[0], eta))

    return np.array(cdf, dtype=[('cdf', np.float), ('energy', np.float)])


def mete_energy_psi_pdf(S, N, E, summary=False):
    '''
    mete_energy_psi_pdf(S, N, E, summary=False)

    METE community energy distribution (Harte 2011)

    Parameters
    ----------
    S : int
        Total number of species in landscape
    N : int
        Total number of individuals in landscape
    E : float
        Total energy ouput of community

    Returns
    -------
    : ndarray (1D)
        If summary is False, returns array with pmf. If summary is True,
        returns the summed log likelihood of all values in n.


    Notes
    -----
    This is the psi distribution found in Table 7.3 of Harte (2011)

    '''

    assert S < N, "S must be less than N"
    assert S > 1, "S must be greater than 1"
    assert N > 0, "N must be greater than 0"
    assert N < E, "N must be less than E"

    
    #Start and stop for root finders
    start = 0.3
    stop = 2

    #Could make a get beta function like Ethan White...
    #Getting beta
    k = np.linspace(1, N, num=N)
    eq = lambda x: sum(x ** k / float(N) * S) -  sum((x ** k) / k)
    x = scipy.optimize.brentq(eq, start, min((sys.float_info[0] / S)**(1/float(N)),\
                              stop), disp=True)
    beta = -m.log(x)
    lambda2 = float(S) / (E - N) #Harte (2011) 7.26
    lambda1 = beta - lambda2
    sigma = lambda1 + (E * lambda2)

    norm = (float(S) / (lambda2 * N)) * (((np.exp(-beta) - np.exp(-beta*(N + 1))) /\
           (1 - np.exp(-beta))) - ((np.exp(-sigma) - np.exp(-sigma*(N + 1))) / \
           (1 - np.exp(-sigma)))) #Harte (2011) 7.22
    
    epi = np.linspace(1, E, num=E*10) #10 is arbitrary
    exp_neg_gamma = np.exp(-(beta + (epi - 1) * lambda2)) #Notation from E.W.

    pdf1 = (float(S) / (N * norm)) * ((exp_neg_gamma  / (1 - exp_neg_gamma)**2) - \
          ((exp_neg_gamma**N / (1 - exp_neg_gamma)) * (N + (exp_neg_gamma / (1 - \
          exp_neg_gamma))))) # Harte (2011) 7.24

    #EW believes equation 7.24 may have an error
    norm_factor = lambda2 / ((np.exp(-beta) - np.exp(-beta * (N + 1))) / (1 - \
                  np.exp(-beta)) - (np.exp(-sigma) - np.exp(-sigma * (N + 1))) / \
                  (1 - np.exp(-sigma)))

    pdf2 = norm_factor *  exp_neg_gamma * (1 - (N + 1) * exp_neg_gamma \
           ** N + N * exp_neg_gamma ** (N + 1)) / \
           (1 - exp_neg_gamma) ** 2 #EW Code

    #Based on some graphical analysis, there is a slight difference between the 
    #pdf as the values get higher. 


    if summary: return -sum(np.log(pdf1))
    else:       return pdf1

    




    




