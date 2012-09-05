#!/usr/bin/python

'''This module provides functions for outputting results of macroeco 
analyses'''


from __future__ import division
import matplotlib.pyplot as plt
import numpy as np
import logging
from macroeco.utils.form_func import output_form, add_field

__author__ = "Mark Wilber"
__copyright__ = "Copyright 2012, Regents of the University of California"
__license__ = None
__version__ = "0.5"
__maintainer__ = "Mark Wilber"
__email__ = "mqw@berkeley.edu"
__status__ = "Development"

class SADOutput(object):
    '''
    This formats and outputs SAD analyses

    '''

    def __init__(self, out_dir):
        '''
        Parameters
        ----------
        out_dir : string
            Output directory of object
        '''
        self.out_dir = out_dir

    def write_summary_table(self, smry, criteria=None):
        '''
        Parameters
        ----------
        smry : dict
            A dictionary as returned by the fucntion compare_summary within
            the CompareDistribution class

        Notes
        -----
        Writes out a formatted txt file to self.out_dir 

        '''
        logging.info('Writing summary table')
        tot_sad = len(smry['obs']['balls'])
        if criteria != None:
            assert len(criteria) == tot_sad, "len(criteria) must  equal" + \
                                      " number of sads under consideration"
        for i in xrange(tot_sad):
            fout = open(self.out_dir + '_summary_table' + str(i) + '.txt', 'w')
            if criteria != None:
                fout.write('CRITERIA: ' + str(criteria[i]) + '\n\n')
            else:
                fout.write('CRITERIA: NONE -- SAD ' + str(i) + '\n\n')
            ob = smry['obs']
            fout.write('EMPIRICAL VALUES:\n' + 'S = ' + str(ob['urns'][i]) + 
                    '\nN = ' + str(ob['balls'][i]) + '\nObserved Nmax= ' + 
                    str(ob['max'][i]) + '\nObserved Rarity = ' +
                    str(ob['tot_min'][i]) + '\n\n')
            for kw in smry.iterkeys():
                if kw != 'obs':
                    dt= smry[kw]
                    fout.write('PREDICTED DISTRIBUTION : ' + kw + '\nS = ' +
                    str(dt['urns'][i]) + '\nN = ' + str(dt['balls'][i]) + 
                    '\nAIC = ' + str(dt['aic'][i]) + '\nDelta_AIC = ' +
                    str(dt['aic_d'][i]) + '\nAIC_weight = ' +
                    str(dt['aic_w'][i]) +'\nPredicted Nmax= ' + 
                    str(dt['max'][i]) + '\nPredicted Rarity = ' + 
                    str(dt['tot_min'][i]) + '\n\n')
            fout.close()

    def plot_rads(self, rads, criteria=None):
        '''
        Plotting the observed and predicted rank abundance distributions

        Parameters
        ----------
        rads : dict
            A dictionary that is returned from the function compare_rads in the
            CompareDistribution class. 

        Notes
        -----
        Saves RAD plots to given out_dir.  Saves as many plots as there are
        SADs.

        '''
        tot_sad = len(rads['obs'])
        recs = make_rec_from_dict(rads, tot_sad)
        if criteria != None:
            assert len(criteria) == tot_sad, "len(criteria) must  equal" + \
                                      " number of sads under consideration"
        for i, data in enumerate(recs):
            names = data.dtype.names
            for nm in names:
                plt.plot(np.arange(1, len(data) + 1), np.sort(data[nm])[::-1],\
                                   '-o')
            plt.legend(names, loc='best')
            if criteria != None:
                plt.title('RAD criteria: ' + str(criteria[i]))
            else:
                plt.title('RAD: plot number ' + str(i))
            plt.semilogy()
            plt.ylabel('log(abundance)')
            plt.xlabel('rank')
            logging.info('Saving figure ' + self.out_dir + '_rank_abundance_' 
                            + str(i))
            plt.savefig(self.out_dir + '_rank_abundance_' + str(i))
            plt.clf()
            output_form(recs[i], self.out_dir + '_rank_abundance_' + str(i))
    
    def plot_cdfs(self, cdfs, obs_sads, criteria=None):
        '''

        Plots observed vs predicted cdfs and returns a csv file with values
        used for plotting.


        Parameters
        ----------
        cdfs : dict
            A dictionary that is returned from the function compare_cdfs in the
            CompareDistribution class. 

        obs_sads : list
            A list of arrays.  The observed sad(s)

        '''
        tot_sad = len(cdfs['obs'])
        recs = make_rec_from_dict(cdfs, tot_sad)
        if criteria != None:
            assert len(criteria) == tot_sad, "len(criteria) must  equal" + \
                                      " number of sads under consideration"
        for i, data in enumerate(recs):
            names = data.dtype.names
            for nm in names:
                plt.plot(np.sort(obs_sads[i]), np.sort(data[nm]), '-o')
            plt.legend(names, loc='best')
            if criteria != None:
                plt.title('CDF criteria: ' + str(criteria[i]))
            else:
                plt.title('CDF: plot number ' + str(i))
            plt.semilogx()
            plt.ylabel('cumulative probability')
            plt.xlabel('abundance')
            logging.info('Saving figure ' + self.out_dir + '_cdf_plot_' 
                            + str(i))
            plt.savefig(self.out_dir + '_cdf_plot_' + str(i))
            plt.clf()
            n_rec = add_field(data, [('n', np.int)])
            n_rec['n'] = obs_sads[i]
            output_form(n_rec, self.out_dir + '_cdf_plot_' + str(i)) 



def make_rec_from_dict(dist_dict, num, dt=np.float):
    '''
    Parameters
    ----------
    dist_dict : dict
        A dictionary with each keyword referencing a list of arrays

    num : int
        Number of rec_arrays to return in list

    '''
    recs = []
    names = list(dist_dict.viewkeys())
    dtype = zip(names, np.repeat(dt, len(names)))
    for i in xrange(num):
        temp = np.empty(len(dist_dict[names[0]][i]), dtype=dtype)
        for kw in dist_dict.iterkeys():
            temp[kw] = dist_dict[kw][i]
        recs.append(temp)
    return recs












           
            
            
                        

            


        



"""
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


"""
