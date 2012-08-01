#!/usr/bin/python

'''
METE energy distributions.

Distributions
-------------
- 'mete_energy_theta_pdf -- The METE intra-specific energy distribution 
                            (Harte 2011)

- 'mete_energy_theta_rank_abund'

- 'mete_energy_psi_pdf' -- The METE community energy distribution

-'mete_energy_psi_rank_abund'

- 'mete_energy_nu_pdf/cdf' -- The METE mean energy per species distribution


References
----------

Harte, J. 2011. Maximum Entropy and Ecology: A Theory of Abundance, 
Distribution, and Energetics. Oxford University Press. 

'''

from __future__ import division
import numpy as np
import scipy.optimize 
import scipy.special
import math as m
import sys

def mete_energy_theta_pdf(epi, S, N, E, n, param_out=True):
    '''
    mete_energy_theta_pdf(epi, S, N, E, n, param_out=True)

    METE intra-specific energy distribution (Harte 2011)

    Parameters
    ----------
    epi : int, float, or array-like object
        The energy values at which to calculate the pdf (1, E]
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
        If param_out is True, returns both the pdf and the estimate for lambda 
        2

    '''
    assert S < N, "S must be less than N"
    assert S > 1, "S must be greater than 1"
    assert N > 0, "N must be greater than 0"
    assert N < E, "N must be less than E"
    assert n < N, "n must be less than N"
    if type(epi) is int or type(epi) is float:
        epi = np.array([epi])
    else:
        epi = np.array(epi)

    lambda2 = float(S) / (E - N) #Harte (2011) 7.26
    
    pdf = (n * lambda2 * np.exp(-lambda2 * n * epi)) / (np.exp(-lambda2 * n)\
          - np.exp(-lambda2 * n * E)) #Harte (2011) 7.25

    if param_out == True:
        return pdf, { 'lambda2': lambda2}
    else:
        return pdf

def mete_energy_nu_pdf(S, N, E, param_out=False):
    '''

    mete_energy_nu_pdf(S, N, E, param_out=False)

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
        If param_out is False, return a pdf.  If param_out is True, returns 
        pdf, eta, and a dictionary with lamda2, beta

    Notes
    -----
    This function uses the nu energy equation found in in Table 7.3 of 
    Harte (2011). Nu is defined as the 'distribution of metabolic rates
    averaged over individuals within species' (Harte 2011). 

    Because of its derivation, this function is discrete!  Only will return
    certain values of e. 

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
    eta_func = lambda n: 1 + (1 / (n * lambda2)) #Harte (2011) pg. 157
    eta = eta_func(np.arange(1, N + 1)[::-1])
    pdf = (1 / np.log(1 / beta)) * ((np.exp(-(beta / (lambda2 * (eta - 1)))) / \
                                    (eta - 1)))
    if param_out == True:
        return pdf, eta, {'lambda2' : lambda2, 'beta' : beta}
    else:
        return pdf
    
def mete_energy_nu_cdf(S, N, E):
    '''
    mete_energy_nu_cdf(S, N, E)

    CDF for mete_energy_nu_pdf

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
    : structured ndarray (1D)
        Structured array contains fields 'cdf' and 'energy'. The field 'cdf' represents
        the cdf value at a given energy level (Prob (emin <= energy <= upppere)

    Notes:
    ------
    '''

    assert S < N, "S must be less than N"
    assert S > 1, "S must be greater than 1"
    assert N > 0, "N must be greater than 0"
    assert N < E, "N must be less than E"

    pdf, epi, params = mete_energy_nu_pdf(S, N, E, param_out=True)
    cdf = []
    for i in xrange(len(pdf)):
        if i < len(pdf) - 1:
            cdf.append(pdf[i] * (epi[i + 1] - epi[i]))
        else:
            cdf.insert(0, 0)
    cdf_struc = np.empty(len(pdf), dtype=[('epi', np.float), ('cdf',\
                         np.float)])
    cdf_struc['cdf'] = np.cumsum(np.array(cdf))
    cdf_struc['epi'] = epi
    return cdf_struc


def mete_energy_psi_pdf(epi, S, N, E, param_out=True):
    '''
    mete_energy_psi_pdf(S, N, E, param_out=True)

    METE community energy distribution (Harte 2011)

    Parameters
    ----------
    epi : int, float, or array-like object
        The energy values at which to calculate the pdf (1, E]
    S : int
        Total number of species in landscape
    N : int
        Total number of individuals in landscape
    E : float
        Total energy output of community

    Returns
    -------
    : ndarray (1D)
        If param_ret is False, returns array with pmf. If param_ret is true,
        returns (pmf, dictionary)
        

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
    x = scipy.optimize.brentq(eq, start, min((sys.float_info[0] / S)**\
                             (1/float(N)), stop), disp=True)
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
    #norm_factor = lambda2 / ((np.exp(-beta) - np.exp(-beta * (N + 1))) / (1 - \
    #              np.exp(-beta)) - (np.exp(-sigma) - np.exp(-sigma * (N + 1))) / \
    #              (1 - np.exp(-sigma)))

    #pdf2 = norm_factor *  exp_neg_gamma * (1 - (N + 1) * exp_neg_gamma \
    #       ** N + N * exp_neg_gamma ** (N + 1)) / \
    #      (1 - exp_neg_gamma) ** 2 #EW Code

    #Based on some graphical analysis, there is a slight difference between the 
    #pdf as the values get higher. 
    if param_out == True: 
        return pdf1, {'lambda2' : lambda2, 'beta' : beta, 'lambda1' : lambda1,\
                      'sigma' : sigma}
    else:       
        return pdf1

def mete_energy_psi_rank_abund(S, N, E):
    '''
    mete_energy_psi_rank_abund(S, N, E)

    Calculates the rank abundance distribution for the METE psi distribution

    Parameters
    ----------
    S : int
        The total number of species in a given landscape
    N : int
        The total number of individuals in a given landscape
    E : int
        The total energy of all individuals in a given landscape with the
        lowest energy being one.

    Returns
    -------
    : structured np.ndarray
        Returns a structured numpy array with dtype = [('rank', np.int),
        ('psi', np.float)] 

    '''

    assert S < N, "S must be less than N"
    assert S > 1, "S must be greater than 1"
    assert N > 0, "N must be greater than 0"
    assert N < E, "N must be less than E"

    #Start and stop for root finders
    start = 0.3
    stop = 2

    #Could make a get beta function like Ethan White...
    k = np.linspace(1, N, num=N)
    eq = lambda x: sum(x ** k / float(N) * S) -  sum((x ** k) / k)
    x = scipy.optimize.brentq(eq, start, min((sys.float_info[0] / S)**\
                             (1/float(N)), stop), disp=True)
    beta = -m.log(x)
    l2 = float(S) / (E - N) #Harte (2011) 7.26
    l1 = beta - l2
    rank = np.arange(1, N + 1)

    def rank_abund_func(r):
        '''
        Uses Harte (2011) eq. 7.37
        '''
        psi = (1 / l2) * np.log(((beta * N) + r - 0.5) / (r - 0.5)) - (l1 / l2)
        return psi

    psi_struc = np.empty(len(rank), dtype=[('rank', np.int), ('psi',\
                                                              np.float)])
    psi_struc['rank'] = rank
    psi_struc['psi'] = rank_abund_func(rank) 
    return psi_struc

def mete_energy_theta_rank_abund(S, N, E, n):
    '''
    mete_energy_theta_rank_abund(S, N, E, n)

    Calculates the rank abundance distribution for the METE theta distribution

    Parameters
    ----------
    S : int
        The total number of species in a given landscape
    N : int
        The total number of individuals in a given landscape
    E : int
        The total energy of all individuals in a given landscape with the
        lowest energy being one.
    n : int
        The total number of inidividuals in a given species

    Returns
    -------
    : structured np.ndarray
        Returns a structured numpy array with dtype = [('rank', np.int),
        ('theta', np.float)] 

    '''

    l2 = float(S) / (E - N)
    rank = np.arange(1, n + 1)
    rank_abund = lambda r : 1 + (1 / (l2 * n)) * np.log( n / (r - 0.5))
    theta_struc = np.empty(len(rank), dtype=[('rank', np.int), ('theta',\
                                                              np.float)])
    theta_struc['rank'] = rank
    theta_struc['theta'] = rank_abund(rank)
    return theta_struc





    
    




    





