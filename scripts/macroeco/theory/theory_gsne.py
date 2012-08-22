#!/usr/bin/python

'''This script contains the GSNE metrics as derived by Harte et al. (2012)'''

from __future__ import division
import numpy as np
import scipy.optimize 
import scipy.special

__author__ = "Mark Wilber"
__copyright__ = "Copyright 2012, Regents of the University of California"
__credits__ = ["John Harte"]
__license__ = None
__version__ = "0.1"
__maintainer__ = "Mark Wilber"
__email__ = "mqw@berkeley.edu"
__status__ = "Development"

def get_lambda1(G, S):
    '''
    This function uses equation 26 to solve for lambda1
    Aprroximation!

    '''
    #TODO: Through assertion 

    assert G < S, 'G must be less than S'
    assert S / G >= 2.8, 'S must be more than 2.8 times greater than G'

    start = 1e-10
    stop = 0.373
    eq = lambda x: 1 / (x * np.log((1 / x))) - (S / G)
    lambda1 = scipy.optimize.brentq(eq, start, stop, disp=True)
    return lambda1

def get_beta(G, N):
    '''
    This function uses equation 28 to solve for beta
    Approximation!

    '''

    assert G < N, 'G must be less than N'
    assert N / G >= 2.8, 'N must be more than 2.8 times greater than G'
    return get_lambda1(G, N)

def lagrange_mult(la, G, S, N, E):
    '''
    '''
    #Assertion Statements needed

    la1 = la[0]
    la2 = la[1]
    la3 = la[2]

    m = np.repeat(np.arange(1, S + 1), N)
    #There is probably a better way to do this
    n = np.arange(1, N + 1)
    for i in xrange(S - 1):
        n = np.concatenate((n, np.arange(1, N + 1)))

    Z = np.exp(-la1 * m) * np.exp(-la2 * m * n) * (-1 / la3) * (1 / m) \
               * (1 / n) * (-np.exp(-la3 * m * n))
    Z_sum = np.sum(Z)

    A = Z * m
    B = Z * m * n
    C = Z * (la3 * m * (n + 1)) / la3
    f1 = np.sum(A) / Z_sum - (S / G)
    f2 = np.sum(B) / Z_sum - (N / G)
    f3 = np.sum(C) / Z_sum - (E / G)
    return (f1, f2, f3)

    """#Constraint S / G
    A = m * np.exp(-la1 * m) * np.exp(-la2 * m * n) * (-1 / la3) * (1 / m) \
               * (1 / n) * (-np.exp(-la3 * m * n))
    constr1 = np.sum(A) / Z_sum - (S / G)

    #Constraint N / G
    B = m * n * np.exp(-la1 * m) * np.exp(-la2 * m * n) * (-1 / la3) * (1 / m)\
               * (1 / n) * (-np.exp(-la3 * m * n))
    constr2 = np.sum(B) / Z_sum - (N / G)

    #Constraint E / G
    C1 = m * n * np.exp(-la1 * m) * np.exp(-la2 * m * n)
    C2 = (-(-la3 * m * n *\
         np.exp(-la3 * m * n)) + np.exp(-la3 * m\
         * n)) / (-la3 * m * n)**2
    C_all = C1 * C2
    constr3 = np.sum(C_all) / Z_sum - (E / G)

    return [constr1, constr2, constr3]"""

def gsne_solver(G, S, N, E):
    '''
    '''
    la1 = G / S
    la2 = G / N
    la3 = G / (E - N)
    la = scipy.optimize.fsolve(lagrange_mult, [5, 5, 5], \
                                args=(G, S, N, E), full_output=True)
    return la





#



#Change name of this distribution later
def gsne_species_per_genera(m, G, S, N):
    '''
    The distribution of the number across all of the genera in the landscape

    m : array-like object
        Array containing desired number of species per genera
    '''

    assert G < S, 'G must be less than S'
    assert S < N, 'S must be less than N'
    if type(m) is int or type(m) is float:
        m = np.array([m])
    else:
        m = np.array(m)

    lambda1 = get_lambda1(G, S)
    pdf = (np.exp(-lambda1 * m)) / (np.log(1 / lambda1) * m)
    return pdf

def gsne_phi_sad(n, G, S, N):
    '''
    The SAD derived from GSNE
    '''
    #NOTE:  This equation does not seem to be properly normalized!
    assert G < S, 'G must be less than S'
    assert S < N, 'S must be less than N'
    if type(n) is int or type(n) is float:
        n = np.array([n])
    else:
        n = np.array(n)

    lambda1 = get_lambda1(G, S)
    beta = get_beta(G, N)
    pdf = (lambda1 * np.exp(-(lambda1 + (beta * n)))) / (n * np.log(1 / beta) \
                                       * (1 - np.exp(-(lambda1 + (beta * n)))))
    return pdf, beta, lambda1

def gsne_K_gad(n, G, S, N):
    '''
    The GAD distribution derived from the GSNE
    '''
    
    assert G < S, 'G must be less than S'
    assert S < N, 'S must be less than N'
    if type(n) is int or type(n) is float:
        n = np.array([n])
    else:
        n = np.array(n)

    lambda1 = get_lambda1(G, S)
    pdf = np.exp((-lambda1 * N) / (n * np.log(1 / lambda1))) / \
                                            (n * np.log(1 / lambda1))
    return pdf








    
