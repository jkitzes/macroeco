#!/usr/bin/python

'''
Contains functions to compare empirical and predicted sars

'''

from __future__ import division
from macroeco import predict_sar
from macroeco import sad_analysis
from macroeco import ssad
from macroeco import empirical
import numpy as np

__author__ = "Mark Wilber"
__copyright__ = "Copyright 2012, Regents of University of California"
__credits__ = "John Harte"
__license__ = None
__version__ = "0.1"
__maintainer__ = "Mark Wilber"
__email__ = "mqw@berkeley.edu"
__status__ = "Development"

def get_obs_pred_sar(datafile, cuts, samples=50, gridded=True):
    '''Gets the observed and predicted sar and returns

    cuts -- either wh_list or divisions.  If sample=False, expects
    divisions, if sample=True expects wh_list

    Maybe I should break this into two functions

    Just doing downscaling for now...

    '''
    sad = sad_analysis.get_gridded_sad_list(datafile, [(1,1)], clean=True)[0][0]
    data = empirical.Patch(datafile)
    anchor_area = data.width * data.height
    if gridded:
        obs_sar = data.sar_grid(cuts)
    else:
        obs_sar = data.sar_sample(cuts, samples)

    a_list = obs_sar[:,1] #areas
    pred_sar = predict_sar.sar_method2(len(sad), np.sum(sad), anchor_area, a_list)
    sar_array = np.empty(len(a_list), dtype=[('obs', np.float), ('pred', np.float),\
                                             ('area', np.float)])
    sar_array['obs'] = obs_sar[:,0]
    sar_array['pred'] = pred_sar['species']
    sar_array['area'] = a_list
    
    return sar_array









    
    



