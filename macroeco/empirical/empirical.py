#!/usr/bin/python

'''
Calculating empirical macroecological metrics on censused patch.

Declares a class Patch that represents a spatially and temporally discrete 
ecological census (such as the STRI BCI forest plot 2005, serpentine grassland 
plot 1998, etc).

Classes
-------
- `Patch` -- a censused area

Patch Methods
-------------
- `sad_grid` -- calculate gridded SAD, SAR, or EAR
- `sad_sample` -- calculate sampled SAD, SAR, or EAR
- `qs_grid` -- calculate commonality for all possible pairs of sub-patches
- `get_sub_sad` -- return SAD, SAR, or EAR for rectangular sub-patch
- `get_sub_abund` -- return vector of abundances for given table
- `get_sp_centers` -- 
- 'get_div_areas' -- return list of areas made by div_list

Misc functions
--------------
- `distance` -- return Euclidean distance between two points
'''

from __future__ import division
import numpy as np
import itertools
from random import choice
from macroeco.data import DataTable

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
    ecological dataset.

    Parameters
    ----------
    datapath : str
        Path to csv file containing census data.
    params : dict
        Dictionary of parameters.

    Attributes
    ----------
    datatable : object of class DataTable
        Object containing patch data and metadata.
    meta_dict : dict
        Dictionary of metadata, containing minimum, maximum, and precision for 
        each column of datatable.
    S : int
        Number of species in patch.
    N : int
        Number of total individuals in patch.
    n0_array : recarray
        Recarray of spp name and abundance.
    width, height : float
        Inclusive width and height of patch, in precisions.
    '''

    def __init__(self, datapath, params):
        '''Initialize object of class Patch. See class documentation.'''

        self.datatable = DataTable(datapath, params)
        self.params = params

        if self.params['count_col'] == None:
            count_col = None
        else:
            count_col = self.datatable.table[params['count_col']]
        self.sad = self.sum_by_spp(self.spp_list, 
                                   self.datatable.table[params['spp_col']],
                                   count_col)

        self.spp_list = (np.unique(self.table[params['spp_col']])).astype(str)
        self.S = self.sad.shape[0]
        self.N = sum(self.sad)

    def grid(self):
        '''
        Divides div_cols into div_list divisions, and sum across sum_col. If no 
        sum_col is given in params, counts individuals of each species.
        
        Parameters
        ----------
        None (but see elements of self.params).

        Returns
        -------
        grid_dict : dict of ndarrays
            Dict with one key, value pair for each element of div_list, where 
            key is the div_list tuple and value is a 2D ndarray with row for 
            each subtable created by division and col for each spp.          
            
        Notes
        -----
        The values of in the tuples of div_list must be factors of the total 
        range of the associated div_col (maximum - minimum + precision) so that 
        every subpatch will have an identical number of survey points.
        '''
        
        grid_dict = {}

        # Get range and number of unique elements of all div_cols
        div_col_ranges = []
        div_col_mins = []
        div_col_maxs = []
        div_unique_vals = []
        for div_col in self.params['div_cols']:
            try:
                div_col_mins.append(self.datatable.meta[(div_col, 'minimum')])
                div_col_maxs.append(self.datatable.meta[(div_col, 'minimum')])

                div_col_ranges.append(\
                    (self.datatable.meta[(div_col, 'maximum')] - 
                     self.datatable.meta[(div_col, 'minimum')] + 
                     self.datatable.meta[(div_col, 'precision')]))
                div_unique_vals.append(None)
            except:
                div_col_ranges.append(None)
                div_unique_vals.append(\
                    len(np.unique(self.datatable.table[div_col])))

        # Loop each division tuple
        for div_tuple in self.params['div_list']:

#            # Loop through each element of the division tuple. Get the
#            # associated col and check whether it has a range. If it does not,
#            # make sure that division is either 1 (or 'none') or
#            # number of unique values (or 'all'). If div_col_ranges, make sure
#            # range is divisible by division.
#            for n, this_elem in enumerate(div_tuple):
#                if not div_col_ranges[n]:  # If no range
#                    if this_elem == 1:
#                        pass
#                    elif this_elem == div_unique_vals[n]:
#                        pass
#                    elif this_elem == 'none':
#                        new_tuple = list(div_tuple)
#                        new_tuple[n] = 1
#                        div_tuple = tuple(new_tuple)
#                    elif this_elem == 'all':
#                        new_tuple = list(div_tuple)
#                        new_tuple[n] = div_unique_vals[n]
#                        div_tuple = tuple(new_tuple)
#                    else:
#                        raise ValueError('Field %s cannot be divided into\
#                                         %s divisions.' % 
#                                         (self.params['div_cols'][n],
#                                          this_elem))
#                else:
#                    assert (divisible(div_col_ranges[n],
#                        self.datatable.meta[(self.params['div_cols'][n], 
#                                             'precision')], this_elem),
#                            'Field %s cannot be divided into %s divisions.' % 
#                            (self.params['div_cols'][n], this_elem))
            
            # Declare array to hold results
            result = np.empty((np.prod(div_tuple), self.S))
            row_counter = 0

            # Loop through each row  - FIX FOR UNIQUES
            for row_tuple in itertools.product(*[range(0, x) for x in 
                                                 div_tuple]):

                # Get conditions for each col
                cond_list = []
                for n, elem in enumerate(row_tuple):
                    this_col = self.params['div_cols'][n]
                    if div_col_ranges[n]:  # If there is a range, do <, >
                        width = rnd(div_col_ranges[n] / div_tuple[n])
                        st = rnd(div_col_mins[n] + elem * width)
                        en = rnd(st + width)
                        cond_list.append((this_col, ">" + str(st)))
                        cond_list.append((this_col, "<" + str(en)))
                    else:  # If no range, do unique subset bool
                        if elem == 'none':
                            pass
                        else:  # DON'T KNOW IF THIS IS A STRING OR NOT
                            cond_list.append((this_col, "=='" + str(elem) + 
                                              "'"))


                # Get subtable

                # Get sub spp count and place in results array

                # Increment row_counter
                row_counter += 1


            for grid_col in xrange(0, div[0]):  # Loop columns of grid
                x_st = rnd(self.x_min + grid_col * sub_width)
                x_en = rnd(x_st + sub_width)

                for grid_row in xrange(0, div[1]):  # Loop rows of grid
                    y_st = rnd(self.y_min + grid_row * sub_height)
                    y_en = rnd(y_st + sub_height)

                    div_result.append(self.get_sub_sad(x_st, x_en, y_st, y_en, 
                                                       summary))

            grid_dict[div_tuple] = result

        return grid_dict


    def sample(self, wh_list, samples, summary = ''):
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
            if wh[0] > self.width or wh[1] > self.height:
                raise ValueError('Width and height of sample patch must be ' +
                                 'less than width and height of full patch') 
            
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
            a division tuple from div_list. Elements are a 2D ndarray with the 
            index of the pair of subpatches (counting row wise from the top 
            left) in the first 2 cols, the distance between the center of the 
            subpatches in the 3rd col, and the commonality in the 4th col.
            
        Notes
        -----
        The values of x and y in the tuples of div_list must be factors of 
        width and height, respectively, so that every subpatch will have an 
        identical number of survey points.
        '''
        # TODO: Make sure that divs divide plot evenly and that divs are of
        # width at least one precision.
        # TODO: Failed with test_dense?
        # TODO: Arg to spit out individ species, not Sorensen

        result = []

        # SAD and sp_cent have the same number of list elements and row
        # structure within arrays
        SAD = self.sad_grid(div_list, summary = '')
        sp_cent = self.get_sp_centers(div_list)

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
                    denom = (0.5 * (sum(spp_in_a) + sum(spp_in_b)))
                    if denom == 0:
                        QS = 0
                    else:
                        QS = sum(spp_in_a * spp_in_b) / denom
                    div_result.append((ind_a, ind_b, dist, QS))

            result.append(np.array(div_result))

        return result

    def sar_grid(self, div_list):
        '''
        Calulates the average SAR for the given div_list

        Parameters
        ----------
        div_list : list of tuples
            Number of divisions of patch along (x, y) axes to use for grid.
            Input of (1, 1) corresponds to no divisions, ie the whole patch.

        Returns
        -------
        : 2D ndarray
            Column one contains the average number of species in a cell of an
            area given in column 2.  Units of column 2 depend on units in 
            metadata.
        '''

        sar_list = self.sad_grid(div_list, summary='SAR')
        sar_means = []
        for sar in sar_list:
            sar_means.append(np.mean(sar))
        sar_array = np.empty((len(div_list), 2), dtype=np.float)
        sar_array[:,0] = np.array(sar_means)
        sar_array[:,1] = np.array(self.get_div_areas(div_list))
        return sar_array

    def sar_sample(self, wh_list, samples):
        '''Calculates the average SAR for the given number of samples
        of each grid in wh_list.

        Parameters
        ----------
        wh_list : list of tuples
            Width and height of subpatches to be placed in patch. Width and 
            height must be less than patch p_width and p_height.
        samples : int
            Number of randomly sampled subpatches to draw for each width/height 
            combination in wh_list.

        Returns
        -------
        : 2D ndarray
            Column one contains the average number of species in a cell of an 
            area given in column 2.  Units of column 2 depend on units in 
            metadata.
        '''

        sar_samples = self.sad_sample(wh_list, samples, summary='SAR')
        sar_means = []
        for sar in sar_samples:
            sar_means.append(np.mean(sar))
        sar_array = np.empty((len(wh_list), 2), dtype=np.float)
        sar_array[:,0] = np.array(sar_means)
        sar_array[:,1] = np.array([x[0] * x[1] for x in wh_list])
        return sar_array

    def sum_by_spp(self, spp_list, spp_col, sum_col=None):
        '''
        Returns a 1D ndarray of length spp counting the number of instances of 
        each spp in spp_col or summing values for each spp in sum_col (if 
        sum_cols is not None).
        '''

        n_spp = len(spp_list)
        result = np.zeros(n_spp)

        if sum_col == None:
            for n in xrange(0, n_spp):
                result[n] = np.sum(spp_col == spp_list[n])
        else:
            for n in xrange(0, n_spp):
                result[n] = np.sum(sum_col[spp_col == spp_list[n]])

        return result

    def get_sp_centers(self, div_list):
        '''
        Get coordinate of center of patches in landscape gridded according to 
        divisions in div_list.
        '''
        # TODO: Did not confirm that this works for x_min and y_min > 0
        sp_centers = []
        for div in div_list:
            sp_width = self.width / float(div[0])
            sp_height = self.height / float(div[1])

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

    def get_div_areas(self, div_list):
        '''
        Get the areas of the cells made by the given div_list
        '''

        area_list = []
        for div in div_list:
            if divisible(self.width, self.xy.meta['precision'], div[0])\
                and divisible(self.height, self.xy.meta['precision'], div[1]):
                
                width = self.width / div[0]
                height = self.height / div[1]
                area = width * height
                area_list.append(area)

            else:
                raise IndexError('Patch not evenly divisible by div.')
        return area_list


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
    '''
    return round(num, 6)
