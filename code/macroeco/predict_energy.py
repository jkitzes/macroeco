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

References
----------

Harte, J. 2011. Maximum Entropy and Ecology: A Theory of Abundance, Distribution,
and Energetics. Oxford University Press. 

'''

import numpy as np
import scipy as sp
import scipy.stats as stats
import scipy.optimize 
import scipy.special
import math as m
import exceptions
import scipy.integrate as integrate

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
    lambda2 = float(S) / (E - N) #Harte (2011) 7.26

    epi = np.linspace(1, E, num=E*10) #Ten is arbitrary
    pdf = (n * lambda2 * np.exp(-lambda2 * n * epi)) / (np.exp(-lambda2 * n)\
          - np.exp(-lambda2 * n * E)) #Harte (2011) 7.25

    if summary: return -sum(np.log(pdf))
    else:       return pdf

def mete_energy_theta_cdf():
    pass

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

    '''

    start = 0.3
    stop = 2

    #Could make a get beta function like Ethan White...
    k = np.linspace(1, N, num=N) 
    eq = lambda x, S, N: (sum(x ** k) / sum((x ** k) / k)) - (float(N)/S)
    x = scipy.optimize.brentq(eq, start, stop, args=(S,N), disp=True)
    beta = -m.log(x)

    lambda2 = float(S) / (E - N) #Harte (2011) 7.26
    lambda1 = beta - lambda2
    sigma = lambda1 + (E * lambda2)

    norm = (float(S) / (lambda2 * N)) * (((np.exp(-beta) - np.exp(-beta*(N + 1))) /\
           (1 - np.exp(-beta))) - ((np.exp(-sigma) - np.exp(-sigma*(N + 1))) / \
           (1 - np.exp(-sigma)))) #Harte (2011) 7.22
    
    epi = np.linspace(1, E, num=E*10) #10 is arbitrary
    exp_neg_gamma = np.exp(-(beta + (epi - 1) * lambda2)) #Notation from E.W.

    #EW believes equation 7.24 may have an error

    pdf1 = (float(S) / (N * norm)) * ((exp_neg_gamma  / (1 - exp_neg_gamma)**2) - \
          ((exp_neg_gamma**N / (1 - exp_neg_gamma)) * (N + (exp_neg_gamma / (1 - \
          exp_neg_gamma))))) # Harte (2011) 7.24
    
    norm_factor = lambda2 / ((np.exp(-beta) - np.exp(-beta * (N + 1))) / (1 - \
                  np.exp(-beta)) - (np.exp(-sigma) - np.exp(-sigma * (N + 1))) / \
                  (1 - np.exp(-sigma)))

    pdf2 = norm_factor *  exp_neg_gamma * (1 - (N + 1) * exp_neg_gamma \
           ** N + N * exp_neg_gamma ** (N + 1)) / \
           (1 - exp_neg_gamma) ** 2 #EW Code

    #Based on some graphical analysis, there is a slight difference between the 
    #pdf as the values get higher. 

    

    return pdf1, pdf2 




    




