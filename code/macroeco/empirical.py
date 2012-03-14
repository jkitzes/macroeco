#!/usr/bin/python

'''
Routines for calculating empirical macroecological metrics.

This module primarily declares a class Patch that represents a single, 
spatially and temporally discrete ecological census (such as the STRI BCI 
forest plot 2005, serpentine grassland plot 1998, etc). Patch has many methods 
that are used to calculate empirical macroecological metrics.

Classes
-------
- `Patch` -- a single censused area

Patch Methods
-------------
- `sad_grid` -- calculate gridded SAD, SAR, or EAR
- `sad_sample` -- calculate sampled SAD, SAR, or EAR
- `qs_grid` -- calculate commonality for all possible pairs of sub-patches
- `get_sub_sad` -- return SAD, SAR, or EAR for rectangular sub-patch
- `get_sub_abund` -- return vector of abundances for given table
- `get_sp_centers` -- 

Misc functions
--------------
- `distance` -- return Euclidean distance between two points
'''

from __future__ import division
import numpy as np
from random import choice
from data import XYTable

__author__ = "Justin Kitzes"
__copyright__ = "Copyright 2012, Regents of the University of California"
__credits__ = ["John Harte"]
__license__ = None
__version__ = "0.1"
__maintainer__ = "Justin Kitzes"
__email__ = "jkitzes@berkeley.edu"
__status__ = "Development"

#
# CLASS PATCH
#
class Patch:
    '''
    A Patch object represents a single spatially and temporally discrete 
    ecological dataset containing information on species abundances.

    Multiple non-contiguous patches or temporal series of data from a single 
    patch should be declared as an object of class Network.

    Parameters
    ----------
    filepath : str
        Path to xytable data file. File must have header row and include cols, 
        at a minimum, named int_code, x, and y.

    Attributes
    ----------
    xy : object of class XYTable
        Object containing patch data and metadata.
    S : int
        Number of species in patch.
    N : int
        Number of total individuals in patch.
    n0_vect : ndarray
        1D array of abundance for each species code.
    x_min, x_max, y_min, y_max : float
        Minimum and maximum values of x and y for patch.
    width, height : float
        Inclusive width and height of patch, in precisions.
    '''

    def __init__(self, datapath):
        ''' Initialize object of class Patch. See class documentation. '''
        self.xy = XYTable(datapath)
        self.set_attributes()


    def set_attributes(self):
        ''' Set attributes other than xy. See class documentation. '''
        xr = self.xy.meta['xrange']
        yr = self.xy.meta['yrange']
        pr = self.xy.meta['precision']

        self.x_min = xr[0]
        self.x_max = rnd(xr[1] + pr)
        self.y_min = yr[0]
        self.y_max = rnd(yr[1] + pr)

        self.width = rnd(xr[1] - xr[0] + pr)
        self.height = rnd(yr[1] - yr[0] + pr)
        
        # S and N are actual num of spp and individs present. n0_vect sums to
        # N  but has max_spp_code + 1 elements, not S elements.
        self.n0_vect = self.get_sub_sad(xr[0], rnd(xr[1] + pr), yr[0], 
                                        rnd(yr[1] + pr), summary = '')
        self.S = sum(self.n0_vect > 0)
        self.N = sum(self.n0_vect)


    def sad_grid(self, div_list, summary = ''):
        '''
        Calculate gridded SAD, SAR, or EAR for patch.

        Divides patch into evenly sized strips along vertical and horizontal 
        dimensions, as specified in div_list, and calculates SAD, SAR, or EAR 
        for all subpatches created by a particular division.
        
        Parameters
        ----------
        div_list : list of tuples
            Number of divisions of patch along (x, y) axes to use for grid.
            Tuple (2, 1), for example, splits the patch with a vertical cut.
        summary : string equal to '', 'SAR', or 'EAR'
            Chooses to summarize results as full SAD, SAR, or EAR. See Returns.

        Returns
        -------
        result : list of ndarrays
            List of same length as div_list, with each element corresponding to 
            a division tuple from div_list. If summary = '', elements are a 2D 
            ndarray with each subpatch in a row and each species in a col. If 
            summary = 'SAR' or 'EAR', elements are a 1D ndarray giving the 
            count of species or endemics in each subpatch.
            
        Notes
        -----
        The values of x and y in the tuples of div_list must be factors of 
        p_width and p_height, respectively, so that every subpatch will have an 
        identical number of survey points.
        '''
        
        result = []

        for div in div_list:  # Loop each division tuple

            # Check that x and y divs are factors of width and height
            x_d = divisible(self.width, self.xy.meta['precision'], div[0])
            y_d = divisible(self.height, self.xy.meta['precision'], div[1])
            if not (x_d and y_d):
                raise IndexError('Patch not evenly divisible by div.')
            
            # Calculate result for this div
            div_result = []
            sub_width = rnd(self.width / div[0])
            sub_height = rnd(self.height / div[1])

            for grid_row in xrange(0, div[0]):  # Loop rows of grid
                x_st = rnd(self.x_min + grid_row * sub_width)
                x_en = rnd(x_st + sub_width)

                for grid_col in xrange(0, div[1]):  # Loop cols of grid
                    y_st = rnd(self.y_min + grid_col * sub_height)
                    y_en = rnd(y_st + sub_height)

                    div_result.append(self.get_sub_sad(x_st, x_en, y_st, y_en, 
                                                       summary))

            result.append(np.array(div_result))

        return result


    def sad_sample(self, wh_list, samples, summary = ''):
        '''
        Calculate sampled SAD, SAR, or EAR for patch.

        Randomly places subpatches of width and height given in wh_list down in 
        patch, and calculates SAD, SAR, or EAR for each of samples subpatches.
        
        Parameters
        ----------
        wh_list : list of tuples
            Width and height of subpatches to be placed in patch. Width and 
            height must be less than patch p_width and p_height.
        samples : int
            Number of randomly sampled subpatches to draw for each width/height 
            combination in wh_list.
        summary : string equal to '', 'SAR', or 'EAR'
            Chooses to summarize results as full SAD, SAR, or EAR. See Returns.

        Returns
        -------
        result : list of ndarrays
            List of same length as wh_list, with each element corresponding to 
            a width/height tuple from wh_list. If summary = '', elements are a 
            2D ndarray with each sampled subpatch in a row and each species in 
            a col. If summary = 'SAR' or 'EAR', elements are a 1D ndarray 
            giving the count of species or endemics in each sampled subpatch.
            
        Notes
        -----
        The values of x and y in the tuples of div_list must be factors of 
        p_width and p_height, respectively, so that every subpatch will have an 
        identical number of survey points.
        '''

        result = []

        for wh in wh_list:  # Loop each width-height tuple
            
            # Calculate result for this width height
            pr = self.xy.meta['precision']
            
            wh_result = []
            sub_w = rnd(wh[0])
            sub_h = rnd(wh[1])
            x_sts = np.arange(self.x_min, rnd(self.x_max - sub_w + pr), pr)
            y_sts = np.arange(self.y_min, rnd(self.y_max - sub_h + pr), pr)

            for s in xrange(0,samples):  # Loop each sample
                x_st = choice(x_sts)
                y_st = choice(y_sts)

                x_en = rnd(x_st + sub_w)
                y_en = rnd(y_st + sub_h)

                wh_result.append(self.get_sub_sad(x_st, x_en, y_st, y_en, 
                                                  summary))

            result.append(np.array(wh_result))

        return result


    def QS_grid(self, div_list):
        '''
        QS_grid(div_list)

        Calculates commonality between pairs of subpatches in a gridded patch. 
        Result is the Sorensen index, equivalent to Chi in Harte et al., also 
        equivalent to 1 - Bray-Curtis dissimilarity.

        Divides patch into evenly sized strips along vertical and horizontal 
        dimensions, as specified in div_list, and calculates commonality for 
        all possible pairs of subpatches created by a particular division.
        
        Parameters
        ----------
        div_list : list of tuples
            Number of divisions of patch along (x, y) axes to use for grid.
            Input of (1, 1) corresponds to no divisions, ie the whole patch.

        Returns
        -------
        result : list of ndarrays
            List of same length as div_list, with each element corresponding to 
            an division tuple from div_list. Elements are a 2D ndarray with the 
            index of the pair of subpatches (counting row wise from the top 
            left) in the first 2 cols, the distance between the center of the 
            subpatches in the 3rd col, and the commonality in the 4th col.
            
        Notes
        -----
        The values of x and y in the tuples of div_list must be factors of 
        p_width and p_height, respectively, so that every subpatch will have an 
        identical number of survey points.
        '''
        # TODO: Make sure that divs divide plot evenly and that divs are of
        # width at least one precision.
        # TODO: Failed with test_dense?
        # TODO: Arg to spit out individ species, not Sorensen

        result = []

        # SAD and sp_cent have the same number of list elements and row
        # structure within arrays
        SAD = self.SAD_grid(div_list, summary = '')
        sp_cent = self._get_sp_centers(div_list)

        for ind_div, div in enumerate(div_list):
            div_result = []
            nsp = div[0] * div[1]

            div_SAD = SAD[ind_div]
            div_sp_cent = sp_cent[ind_div]

            for ind_a in xrange(0, nsp - 1):
                spp_in_a = (div_SAD[ind_a,:] > 0)
                for ind_b in xrange(ind_a + 1, nsp):
                    spp_in_b = (div_SAD[ind_b,:] > 0)
                    dist = distance(div_sp_cent[ind_a,:], div_sp_cent[ind_b,:])
                    QS = sum(spp_in_a * spp_in_b) / (0.5 * (sum(spp_in_a) + 
                                                            sum(spp_in_b)))
                    div_result.append((ind_a, ind_b, dist, QS))

            result.append(np.array(div_result))

        return result


    def get_sub_sad(self, x_st, x_en, y_st, y_en, summary):
        '''
        Calculates a SAD, SAR, or EAR (according to summary, see class 
        docstring) for a subpatch of known starting and ending coordinates.
        '''
        
        sub_table = self.xy.get_sub_table(x_st, x_en, y_st, y_en)
        sub_abund = self.get_sub_abund(sub_table)

        if summary is 'SAR': 
            return sum(sub_abund > 0)
        elif summary is 'EAR':  # Don't count if self.total_abund = 0
            return sum((sub_abund == self.n0_vect) * \
                       (self.n0_vect != 0))
        else:
            return sub_abund


    def get_sub_abund(self, table):
        '''
        Calculate abundance of each species in a xytable.
        
        If the maximum spp code is much higher than S - 1, this will be
        inefficient.
        '''
        abund = np.zeros(self.xy.max_spp_code + 1)
        for row in table:
            abund[row[self.xy.col_spp_code]] += row[self.xy.col_count]
        return abund


    def get_sp_centers(self, div_list):
        '''
        Get coordinate of center of patches in landscape gridded according to 
        divisions in div_list.
        '''
        # TODO: Did not confirm that this works for x_min and y_min > 0
        sp_centers = []
        for div in div_list:
            sp_width = self.p_width / float(div[0])
            sp_height = self.p_height / float(div[1])

            div_sp_cent = []
            for sp_x in xrange(0, div[0]):  # Same sp order as SAD_grid
                x_origin = (self.x_min + sp_x * sp_width)
                x_cent = x_origin + 0.5 * (sp_width - self.precision)
                for sp_y in xrange(0, div[1]):
                    y_origin = (self.y_min + sp_y * sp_height)
                    y_cent = y_origin + 0.5 * (sp_height - self.precision)

                    div_sp_cent.append((x_cent, y_cent))

            sp_centers.append(np.array(div_sp_cent))

        return sp_centers




#
# MISC FUNCTIONS
#

def distance(pt1, pt2):
    ''' Calculate Euclidean distance between two points '''
    return np.sqrt((pt1[0] - pt2[0]) ** 2 + (pt1[1] - pt2[1]) ** 2)

def divisible(dividend, precision, divisor, tol = 1e-9):
    '''
    Check if dividend (here width or height of patch) is evenly divisible by 
    divisor (here a number of patch divs) while accounting for floating point 
    rounding issues.
    '''
    if divisor == 0:
        return False
    if divisor > round(dividend / precision):
        return False

    quot_raw = (dividend / precision) / divisor
    quot_round = round(quot_raw)
    diff = abs(quot_raw - quot_round)

    if diff < tol:
        return True
    else:
        return False

def rnd(num):
    '''
    Round num to number of decimal places in precision. Used to avoid issues 
    with floating points in the patch and subpatch width and height that make 
    subpatches not lie exactly on even divisions of patch.

    WARNING: The current implementation requires that precision be no less than 
    0.01! I don't know a good way to check for this.
    '''
    return round(num, 2)
