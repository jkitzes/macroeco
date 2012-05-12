#!/usr/bin/python

'''This module provides functions for graphing results of macroeco analyses'''
from __future__ import division
import matplotlib.pyplot as plt
import numpy as np
import logging
from macroeco.workflow import Workflow, loggername
from macroeco import sad_analysis
from macroeco import form_func
from macroeco import predict_sad
from scipy.interpolate import UnivariateSpline
import os

__author__ = "Mark Wilber"
__copyright__ = "Copyright 2012, Regents of the University of California"
__license__ = None
__version__ = "0.5"
__maintainer__ = "Mark Wilber"
__email__ = "mqw@berkeley.edu"
__status__ = "Development"

module_logger = logging.getLogger(loggername)

def sad_cdf_obs_pred_plot(sad, outputID, params={}, interactive=False):
    '''Makes a plot of observed vs predicted cdf values
    based on the given sad

    Parameters
    ----------
    sad : ndarray
        SAD
    outputID : string
        The name of the saved plot
    params : dict
        If empty uses default values, if not incorporates given params
        into params into dict
    interactive : bool
        If True, figure is shown in interactive window.  
        If false, no figure is shown.

    Notes
    -----
    Saves plot to current working directory. Default distribution to graph
    is METE. Need top specify in parameters.xml if another distribution is
    desired.

    '''

    #Allowing the user to set parameters if they pass some in
    prm = {'clr_jit':'grey', 'clr_sct':'black', 'ln_clr':'grey', 'jit_scl':0.007,\
          'ylbl':'Observed cdf', 'xlbl':'Predicted cdf',\
          'title': 'Observed vs. predicted values of SAD cdf', 'distr':'mete'}
    if len(params) != 0:
        inter = set(params.viewkeys()).intersection(set(prm.viewkeys()))
        if len(inter) != 0:
            module_logger.debug('Setting parameters ' + str(inter))
            for key in inter:
                if key is 'jit_scl':
                    prm[key] = eval(params[key])
                else:
                    prm[key] = params[key]
          
    cdf = sad_analysis.get_obs_vs_pred_cdf(sad, prm['distr'])
    jitt_x, jitt_y = jitter(cdf['pred'], cdf['obs'], jitter_scale=prm['jit_scl'])
    plt.plot([0] + list(cdf['obs']),[0] +  list(cdf['obs']),\
                                            color=prm['ln_clr'], linestyle='--')
        
    plt.scatter(jitt_x, jitt_y, color=prm['clr_jit']) 
    plt.scatter(cdf['pred'], cdf['obs'], color=prm['clr_sct'])
    '''for s in xrange(len(cdf)):
        plt.text(cdf['pred'][s], cdf['obs'][s] + 0.01, "n=" + str(cdf['n'][s]),\
                                fontsize=7)'''
    plt.ylabel(prm['ylbl'])
    plt.xlabel(prm['xlbl'])
    plt.title(prm['title'])
    plt.xlim(-0.1,1.1)
    plt.ylim(-0.1,1.1)
    plt.legend(('Slope 1', 'All species cum.\nprob. (jittered)',\
                    'Cum. prob.\nfor given n'),loc='best', prop={'size':12})
    plt.savefig(outputID + '_cdf')
    form_func.output_form(cdf, outputID + '_cdf.csv')
    module_logger.debug('Plot and csv saved to cwd')

    if interactive:
        plt.show()
    else:
        plt.clf()

def sad_rank_abund_plot(sad, outputID, params={}, interactive=False):
    '''Makes a plot of observed vs predicted rank abundance curves based
    on given sad

    Parameters
    ----------
    sad : ndarray
        SAD
    outputID : string
        The name of the saved plot
    params : dict
        If empty uses default values, if not incorporates given params
        into params into dict
    interactive : bool
        If True, figure is shown in interactive window.  
        If false, no figure is shown.

    Notes
    -----
    Saves plot to current working directory. Default distribution to graph
    is METE. Need top specify in parameters.xml if another distribution is
    desired.


    '''
    prm = {'distr':'mete', 'logx':False, 'logy':True}
    if len(params) != 0:
        inter = set(params.viewkeys()).intersection(set(prm.viewkeys()))
        if len(inter) != 0:
            module_logger.debug('Setting parameters ' + str(inter))
            for key in inter:
                if key is 'logx' or key is 'logy':
                    prm[key] = eval(params[key])
                else:
                    prm[key] = params[key]
    #With support [1,N]
    obs_pred_abund = sad_analysis.get_obs_pred_abund(sad, prm['distr'])
    rank = np.arange(1, len(sad) + 1)
    plt.plot(rank, obs_pred_abund['predicted'][::-1], 'o-', color='gray')
    plt.plot(rank, obs_pred_abund['observed'][::-1], 'o-', color='black')
    plt.title('Observed and predicted rank abundance plot')
    plt.xlabel('Rank')
    plt.ylabel('Abundance')
    plt.legend(('Predicted', 'Observed'), loc='best')
    if prm['logx']:
        plt.semilogx()
        plt.xlabel('log(Rank)')
    if prm['logy']:
        plt.semilogy()
        plt.ylabel('log(Abundance)')
    plt.savefig(outputID + '_rank_abund')
    form_func.output_form(obs_pred_abund, outputID + '_rank_abund.csv')
    module_logger.debug('Plot and csv saved to cwd')

    if interactive:
        plt.show()
    else:
        plt.clf()


def obs_pred_rarity(sad, outputID, params={}, interactive=False):
    '''Makes a bar graph of observed and predicted rarity

     Parameters
    ----------
    sad : ndarray
        SAD
    outputID : string
        The name of the saved plot
    params : dict
        If empty uses default values, if not incorporates given params
        into params into dict
    interactive : bool
        If True, figure is shown in interactive window.  
        If False, no figure is shown.

    Notes
    -----
    Saves plot to current working directory. Default distribution to graph
    is METE. Need top specify in parameters.xml if another distribution is
    desired.

    '''
    prm = {'distr':'mete', 'rarity':10}
    if len(params) != 0:
        inter = set(params.viewkeys()).intersection(set(prm.viewkeys()))
        if len(inter) != 0:
            module_logger.debug('Setting parameters ' + str(inter))
            for key in inter:
                if key is 'rarity':
                    prm[key] = eval(params[key])
                else:
                    prm[key] = params[key]
    #TODO: remove ticks from plots!
    obs_pred = sad_analysis.get_obs_and_pred_rarity(sad, prm['distr'], n=prm['rarity'])
    x = np.array([0.5,1.5])
    plt.bar(x, list(obs_pred), width=0.5, color='gray')
    plt.xticks(x + 0.25, ('Obs', 'Pred'))
    plt.ylabel('Number of species less than ' + str(prm['rarity']))
    plt.ylim(0, max(obs_pred) + .25*(max(obs_pred)))
    plt.text(.75 - .05, obs_pred[0] + 0.5, "n = " + str(obs_pred[0]))
    plt.text(1.75 - .05, obs_pred[1] + 0.5, "n = " + str(obs_pred[1]))
    plt.title("Observed vs. predicted rarity")
    plt.savefig(outputID + '_rarity')
    
    if interactive:
        plt.show()
    else:
        plt.clf()

def obs_pred_Nmax(sad, outputID, params={}, interactive=False):
    '''Makes a bar graph of obeserved and predicted N max

     Parameters
    ----------
    sad : ndarray
        SAD
    outputID : string
        The name of the saved plot
    params : dict
        If empty uses default values, if not incorporates given params
        into params into dict
    interactive : bool
        If True, figure is shown in interactive window.  
        If False, no figure is shown.

    Notes
    -----
    Saves plot to current working directory. Default distribution to graph
    is METE. Need top specify in parameters.xml if another distribution is
    desired.

    '''
    prm = {'distr':'mete'}
    if len(params) != 0:
        inter = set(params.viewkeys()).intersection(set(prm.viewkeys()))
        if len(inter) != 0:
            module_logger.debug('Setting parameters ' + str(inter))
            for key in inter:
                prm[key] = params[key]

    #TODO: remove ticks from plots! Add param manipulations
    obs_pred = sad_analysis.get_obs_and_pred_Nmax(sad, prm['distr'])
    x = np.array([0.5,1.5])
    plt.bar(x, list(obs_pred), width=0.5, color='gray')
    plt.xticks(x + 0.25, ('Obs', 'Pred'))
    plt.ylabel('Nmax of SAD')
    plt.ylim(0, max(obs_pred) + .25 * (max(obs_pred)))
    plt.text(.75 - .05, obs_pred[0] + 0.5, "n = " + str(obs_pred[0]))
    plt.text(1.75 - .05, obs_pred[1] + 0.5, "n = " + str(obs_pred[1]))
    plt.title("Observed vs. predicted N max")
    plt.savefig(outputID + '_Nmax')
    
    if interactive:
        plt.show()
    else:
        plt.clf()

def common_plot(cmn_arrays, outputID, params={}, interactive=False):
    '''To be used directly with get_common_arrays() in cmn_analysis
    This function creates a plot of chi-squared as a function of
    A/D**2.

    Parameters
    ----------
    cmn_arrays : list
        list of structured arrays as outputed by get_common_arrays() in 
        cmn_analysis.py
    outputID : string
        The name of the saved plot
    params : dict
        If empty uses default values, if not incorporates given params
        into params into dict
    interactive : bool
        If True, figure is shown in interactive window.  
        If False, no figure is shown.




    '''
    areas = []

    for data in cmn_arrays:
        x = data['dist']
        y = data['cmn']
        s = UnivariateSpline(x, y)
        xs = np.linspace(np.min(x), np.max(x), num=100*len(x))
        ys = s(xs)
        plt.plot(xs**2, ys)
        plt.plot(x**2, y)
        areas.append("Area = " + str(data['area'][0]) + " (spline)")
        areas.append("Area = " + str(data['area'][0]) + " (raw)")

    plt.xlabel('Distance')
    plt.ylabel('Commonality')
    plt.legend(tuple(areas), loc='best', prop={'size':9})
    plt.title('Commonality and distance')
    plt.savefig(outputID + "_commonality")
    for i, data in enumerate(cmn_arrays):
        form_func.output_form(data, outputID + '_ ' + str(eval(params['grid'])[i]))
    module_logger.debug('Commonality plot and csv files saved')
    
    if interactive:
        plt.show()
    else:
        plt.clf()

def common_ADsquared_plot(cmn, outputID, params={}, areas=[], interactive=False):
    '''Generates a plot of commonality as a function of Area over Distance squared


    '''
    #Getting different symbols for scatter plot
    def_col = ['red','green', 'blue', 'orange', 'yellow', 'black', 'grey']
    color = []
    for i in xrange(len(cmn)):
        col = np.random.random_sample(3)
        color.append(col)

    areas = []
    for i, strc in enumerate(cmn): #iterating through each area
        scale = 20 #scale for truncation
        ind = ((strc['dist'] ** 2 / strc['area']) > scale)
        trun_s = strc[ind]
        if i < len(def_col):
            plt.scatter((trun_s['dist'] ** 2 / trun_s['area']), trun_s['cmn'],\
                                color=def_col[i])
        else:
            plt.scatter((trun_s['dist'] ** 2 / trun_s['area']), trun_s['cmn'],\
                                color=color[i])
        areas.append('Area = ' + str(trun_s['area'][0]))
        form_func.output_form(strc, outputID + '_unique_area_' + str(i))
    plt.xlabel("Distance ^ 2 / Area")
    plt.ylabel("Commonality")
    plt.title("Commonality plotted as a function of distance^2 over area (m^2)")
    plt.legend(tuple(areas), loc='best', prop={'size':9})
    plt.savefig(outputID + '_DAsq_linear')
    plt.xlabel("log(Distance ^ 2 / Area)")
    plt.semilogx()
    plt.savefig(outputID + '_DAsq_linearlog')
    module_logger.debug("Plots and csv's saved to " + os.getcwd())
    if interactive:
        plt.show()
    else:
        plt.clf()









    '''#Less than 0.5 because D**2 >> A
    scale = 0.05
    ind = (a_d['A/D**2'] < scale)
    x = 1 / a_d['A/D**2'][ind] #invert to see what happens as D increases
    y = a_d['cmn'][ind]

    bins = np.linspace(1 / scale, np.max(x), num=11)
    xdig = np.digitize(x, bins)
    unq_xdig = np.unique(xdig)
    y_smooth = []
    x_smooth = []
    for i in unq_xdig:
        ind = (i == xdig)
        y_smooth.append(sum(y[ind]) / sum(ind))
        x_smooth.append(sum(x[ind]) / sum(ind))
        
    plt.scatter(x, y, color='grey', alpha=0.5 )
    plt.plot(x_smooth, y_smooth, '-o', color='black')
    plt.xlabel(" log(Distance ^ 2 / Area))")
    plt.ylabel("log(Commonality)")
    if len(areas) == 0:
        plt.title("Commonality plotted as a function of distance^2 over area")
    else:
        plt.title("Commonality plotted as a function of distance^2 over area:\n\
                   Areas (m^2): " + str(areas)) 
    plt.legend(("Binned data points", "All data points"), loc='best', prop={'size':9})
    plt.savefig(outputID + "_ADsquared")
    form_func.output_form(a_d, outputID + "_ADsquared.csv")
    module_logger.debug("Plot and csv saved")

    if interactive:
        plt.show()
    else:
        plt.clf()'''


    



def jitter(x, y, jitter_scale=0.001):
    '''Function takes in x and y values and jitters
    the identical points (x,y) in the x direction
    
    Parameters
    ----------
    x : ndarray
      x points
    y : ndarray
      y points
    jitter_scale : float
        The distance a point is jittered in the x direction

    Returns
    -------
    : ndarrays
        returns jittered x and y as separate ndarrays



    '''

    #TODO: Add option to jitter in y direction
    assert len(x) == len(y), "x and y must be the same length"
    jitt_x = np.copy(x)
    jitt_y = np.copy(y)
    xy_tuple = []
    for i in xrange(len(x)):
        xy_tuple.append((x[i], y[i]))

    xy_tuple = np.array(xy_tuple)
    unq_xy = np.unique(xy_tuple)
    for xy in unq_xy:
        indices = np.where(xy_tuple == xy)[0]
        if len(indices > 1):
            for ind in xrange(len(indices)):
                if ind > 0:
                    jitt_x[indices[ind]] += (jitter_scale * ind)
    return jitt_x, jitt_y



