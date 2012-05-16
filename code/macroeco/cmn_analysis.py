#!/usr/bin/python

'''Script that uses reproducible workflow to test commonality predictions
on data sets.
'''
from __future__ import division
import sys, os
from macroeco import empirical
import numpy as np

__author__ = "Mark Wilber"
__copyright__ = "Copyright 2012, Regents of the University of California"
__license__ = None
__version__ = "0.5"
__maintainer__ = "Mark Wilber"
__email__ = "mqw@berkeley.edu"
__status__ = "Development"

def get_common_arrays(patch, grids):
    '''This function calulates the commonality for the given grids and
    gets the average commonality for a each distance within a grid.
    The length of the returned list is the same length as 
    grids.

    Parameters
    ----------
    patch : Patch object
          a Patch object from macroeco.empirical
    grids : list
        A list of tuples with desired cuts
    
    Returns
    -------
    : list
        A list of structured arrays with 
        dtype=[('dist', np.float), ('cmn', np.float), ('area', np.float)]
        Dist = Distance, cmn = commonality
        Length of list is same length as grids


    '''
    #Check that grids does not have (1,1)
    index = str(grids).find('(1, 1)')
    if index != -1:
        raise ValueError('Grid of (1,1) has no commonality')
    areas = patch.get_div_areas(grids)
    common = patch.QS_grid(grids)
    struc_list = []
    
    

    for i in xrange(len(common)):
        dist = np.round(common[i][:,2], decimals=7)#To account for rounding issues
        cmn = common[i][:,3]
        unq_dist = np.unique(dist)
        cmn_average = []
        for d in unq_dist:
            ind = (d == dist)
            cmn_average.append(sum(cmn[ind]) / sum(ind))

        dist_cmn = np.empty(len(unq_dist), dtype=[('dist', np.float),\
                                                  ('cmn', np.float),\
                                                  ('area', np.float)])
        dist_cmn['dist'] = unq_dist
        dist_cmn['cmn'] = np.array(cmn_average)
        dist_cmn['area'] = areas[i]
        struc_list.append(dist_cmn)

    return struc_list












