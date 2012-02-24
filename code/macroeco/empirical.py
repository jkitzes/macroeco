#!/usr/bin/python

'''
Routines for calculating empirical macroecological metrics.

This module primarily declares a class Patch that represents a single, 
spatially located ecological survey (such as the STRI BCI forest plot, 
serpentine grassland plot, etc). Patch has many methods that are used to 
calculate macroecological metrics.

Classes
-------
- `Patch` -- a single spatially and temporally identified data set

Patch Methods
-------------
- `SAD_grid` -- calculate gridded SAD, SAR, or EAR
- `SAD_sample` -- calculate sampled SAD, SAR, or EAR
- `QS_grid` -- calculate gridded commonality for pairs of subpatches
- `_sub_SAD` -- calculate SAD, spp count, or endemics count for subset of patch
- `_get_nspp` -- count number of species in sparse or dense data
- `_get_sparse_bool` -- return True if patch data is in sparse format
- `_get_total_abund` -- count total abundance of each species in patch

Misc functions
--------------
- `distance` -- return Euclidean distance between two points
'''


import numpy as np
from random import choice
import data
reload(data)


__author__ = "Justin Kitzes"
__copyright__ = "Copyright 2012, Regents of the University of California"
__credits__ = ["John Harte"]
__license__ = None
__version__ = "0.1"
__maintainer__ = "Justin Kitzes"
__email__ = "jkitzes@berkeley.edu"
__status__ = "Development"

# Temporarily define sample sparse and dense plots for testing
test_dense = np.array((0,0,0,1,1,3,0,0,4,1,0,0,0,2,2,1)*2).reshape(4,4,2)
test_dense[1,0,0] = 0
test_sparse = data.dense_to_sparse(test_dense, 1)
test_xy = test_sparse[:,0:3]
test_xy = np.vstack((test_xy,test_xy[0:8,:]))


#
# CLASS PATCH
#
class Patch:
    '''
    Patch(data, x_minmax, y_minmax, precision)

    A patch object represents a single spatially and temporally discrete 
    ecological dataset containing information on species abundances.

    Multiple non-contiguous patches or temporal series of data from a single 
    patch should be declared as class Network.

    Parameters
    ----------
    (for the __new__ method)

    data : sparse or dense ndarray
        See data.py for routines for properly formatting the data array
    x_minmax : tuple or list
        Minimum and maximum survey points along x-axis
    y_minmax : tuple or list
        Minimum and maximum survey points along y-axis
    precision : float
        Precision of census data (ie, minimum mapping distance between points)

    Attributes
    ----------
    data : ndarray
        Array specifying locations and abundances of all species in patch
    precision : float
        Smallest survey precision (ie, minimum mapping distance between points)
    sparse : bool
        True if data is of form sparse, False if data is of form dense
    nspp : int
        Total number of species in patch
    total_abund : ndarray
        1D array of length nspp giving total abundances of each species
    x_min, x_max, y_min, y_max : float
        x and y coordinates of min and max survey point in patch
    p_width, p_height : float
        Inclusive width and height of patch, in precisions
        Equal to number of survey points along x and y axis
    '''
    # TODO: Configure to work with 3D data sets

    def __init__(self, data, x_minmax, y_minmax, precision):
        '''
        Initialize object of class Patch. See class documentation.
        '''
        # TODO: Error checking for correct plot type
        # TODO: Take extracted metadata fields and loaded data as args

        self.data = data
        self.precision = precision
        self.sparse = self._get_sparse_bool()
        self.sppcodes = self._get_sppcodes()
        self.nspp = len(self.sppcodes)
        self.total_abund = self._get_total_abund()

        self.x_min = x_minmax[0]
        self.x_max = x_minmax[1]
        self.p_width = self.x_max - self.x_min + precision
        self.y_min = y_minmax[0]
        self.y_max = y_minmax[1]
        self.p_height = self.y_max - self.y_min + precision

    def SAD_grid(self, div_list, summary = ''):
        '''
        SAD_grid(div_list, summary = '')

        Calculate gridded SAD, SAR, or EAR for patch.

        Divides patch into evenly sized strips along vertical and horizontal 
        dimensions, as specified in div_list, and calculates SAD, SAR, or EAR 
        for all subpatches created by a particular division.
        
        Parameters
        ----------
        div_list : list of tuples
            Number of divisions of patch along (x, y) axes to use for grid.
            Input of (1, 1) corresponds to no divisions, ie the whole patch.
        summary : string equal to '', 'SAR', or 'EAR'
            Chooses to summarize results as full SAD, SAR, or EAR. See Returns.

        Returns
        -------
        result : list of ndarrays
            List of same length as div_list, with each element corresponding to 
            an division tuple from div_list. If summary = '', elements are a 2D 
            ndarray with each subpatch in a row and each species in a col. If 
            summary = 'SAR' or 'EAR', elements are a 1D ndarray giving the 
            count of species or endemics in each subpatch.
            
        Notes
        -----
        The values of x and y in the tuples of div_list must be factors of 
        p_width and p_height, respectively, so that every subpatch will have an 
        identical number of survey points.
        '''
        # TODO: Error check that div must be >= 1

        result = []

        for div in div_list:  # Loop each division tuple
            div_result = []
            sp_width = self.p_width / float(div[0])
            sp_height = self.p_height / float(div[1])
            # TODO: Error check that x and y strips divide dimensions evenly - 
            # use remainder function on *_max + precision and ensure 0.

            for x_div_count in xrange(0,div[0]):  # Loop x_divs
                x_st = self.x_min + x_div_count * sp_width
                x_en = x_st + sp_width

                for y_div_count in xrange(0,div[1]):  # Loop y_divs
                    y_st = self.y_min + y_div_count * sp_height
                    y_en = y_st + sp_height

                    div_result.append(self._sub_SAD(x_st, x_en, y_st, y_en, 
                                                    summary))

            result.append(np.array(div_result))

        return result


    def SAD_sample(self, wh_list, samples, summary = ''):
        '''
        SAD_sample(wh_list, samples, summary = '')

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
        # TODO: Check that width and height < pwidth and pheight
        result = []

        for wh in wh_list:  # Loop each width-height tuple
            wh_result = []
            sp_width = wh[0]
            sp_height = wh[1]
            x_origins = np.arange(self.x_min, self.x_max - sp_width + 
                                  self.precision, self.precision)
            y_origins = np.arange(self.y_min, self.y_max - sp_width + 
                                  self.precision, self.precision)

            for s in xrange(0,samples):  # Loop each sample
                # TODO: Currently fails for sp_width = whole plot
                x_st = choice(x_origins)
                y_st = choice(y_origins)

                x_en = x_st + sp_width
                y_en = y_st + sp_height

                wh_result.append(self._sub_SAD(x_st, x_en, y_st, y_en, 
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


    def _sub_SAD(self, x_st, x_en, y_st, y_en, summary):
        '''
        Calculates a SAD, SAR, or EAR (according to summary, see class 
        docstring) for a subpatch of known starting and ending coordinates. 
        Only works for rectangles.
        '''
        
        if self.sparse:
            # TODO: Use masked array instead of this bool?
            in_sp = np.all([self.data[:,1] >= x_st, self.data[:,1] < x_en,
                           self.data[:,2] >= y_st, self.data[:,2] < y_en],
                           axis = 0)
            sp_data = self.data[in_sp]
            sub_abund = self._get_sparse_abund(sp_data)
                
        else:
            sub_abund = self.data[y_st:y_en, x_st:x_en, 
                                  :].sum(axis=0).sum(axis=0)

        if summary is 'SAR':
            return sum(sub_abund > 0)
        elif summary is 'EAR': # Don't count if self.total_abund = 0
            return sum((sub_abund == self.total_abund) * \
                       (self.total_abund != 0))
        else:
            return sub_abund


    def _get_sp_centers(self, div_list):
        '''
        Get coordinate of center of patches in landscape gridded according to 
        divisions in div_list. Works for both sparse and dense data.
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


    def _get_total_abund(self):
        ''' Calculate total abundance of each species in entire patch '''
        if self.sparse:
            return self._get_sparse_abund(self.data)
        else:
            self.total_abund = self.data.sum(axis=0).sum(axis=0)


    def _get_sppcodes(self):
        ''' Get array of codes of all species in self.data '''
        # TODO: Dense here still counts 3rd dim 'floors' with 0 individs
        if self.sparse:
            return np.unique(self.data[:,0]).astype(int)
        else:
            return np.size(self.data, 2)  # Size of 3rd dim


    def _get_sparse_bool(self):
        ''' Return true if data is sparse, false if dense '''
        return len(self.data.shape) == 2  # If data is 2D, is sparse


    def _get_sparse_abund(self, data):
        ''' Calculate abundance of each species in a sparse patch '''
        # If the maximum spp code is much higher than S - 1, this will be
        # inefficient.
        abund = np.zeros(max(self.sppcodes) + 1)
        for row in data:
            abund[row[0]] += row[3]
        return abund[self.sppcodes] # Only return codes in sppcodes


#
# MISC FUNCTIONS
#

def distance(pt1, pt2):
    ''' Calculate Euclidean distance between two points '''
    return np.sqrt((pt1[0] - pt2[0]) ** 2 + (pt1[1] - pt2[1]) ** 2)
