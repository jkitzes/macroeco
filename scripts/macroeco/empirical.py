#!/usr/bin/python

'''
Calculating macroecological metrics for empirical or theoretical patch. Patch 
is interpreted broadly as any temporally and spatially defined census.

Classes
-------
- `EPatch` -- empirical metrics for census data
- `TPatch` -- parallel predicted/theoretical/fit metrics

Patch Methods
-------------
- `sad` -- calculate species abundance distribution (grid or sample)
- `sar` -- calculate species-area relationship (grid or sample)
- `ear` -- calculate endemics-area relationship (grid or sample)
- `qs` -- calculate commonality between sub-patches (grid)
- `ssad` -- calculate species-level spatial abundance distrib (grid or sample)

- `get_sp_centers` --
- 'get_div_areas' -- return list of areas made by div_list

Misc functions
--------------
- `distance` -- return Euclidean distance between two points
'''

from __future__ import division
import numpy as np
import itertools
from copy import deepcopy
from random import choice
from data import DataTable

__author__ = "Justin Kitzes"
__copyright__ = "Copyright 2012, Regents of the University of California"
__credits__ = ["John Harte"]
__license__ = None
__version__ = "0.1"
__maintainer__ = "Justin Kitzes"
__email__ = "jkitzes@berkeley.edu"
__status__ = "Development"


class EPatch:
    '''
    An object representing an empirical census.

    Parameters
    ----------
    data_path : str
        Path to csv file containing census data.
    subset : dict
        Dictionary of permanent subset to data, {'column_name': 'condition'}, 
        which will limit all analysis to records in which column_name meets the 
        condition, ie, {'year': '==2005', 'x': ('>20', '<40')} restricts 
        analysis to year 2005 and x values between 20 and 40. These conditions 
        can also be passed to the invididual methods, but subsetting the data 
        table up front may save analysis time.

    Attributes
    ----------
    data_table : object of class DataTable
        Object containing patch data and metadata.

    '''

    def __init__(self, datapath, subset = {}):
        '''Initialize object of class Patch. See class documentation.'''

        self.data_table = DataTable(datapath)
        self.data_table.table = self.data_table.get_subtable(subset)


    def sad(self, criteria):
        '''
        Calculates an empirical species abundance distribution given criteria.

        Parameters
        ----------
        criteria : dict
            Dictionary of form {column_name: value}. Must contain a key with a 
            value of 'species' indicating the column with species identifiers 
            (this column must be type categorical in metadata). If a column 
            giving the counts of species found at a point is also in the data, 
            a key with the value 'counts' should also be given.

            Value has a different meaning depending on column type:
            - metric - number of divisions of data along this axis, int/float
            - categorical - 'split' calculates each category separately
        
        Returns
        -------
        spp_list : ndarray
            1D array listing identifiers for species in the same order as they 
            appear in arrays found in result.
        result : dict
            List of tuples containing results, where first element is 
            dictionary of criteria for this calculation and second element is a 
            1D ndarray of length species containing the abundance for each 
            species.
        '''

        spp_list, spp_col, count_col, combinations = \
            self.parse_criteria(criteria)

        result = []
        for comb in combinations:
            
            subtable = self.data_table.get_subtable(comb)

            sad_list = []
            for species in spp_list:
                spp_subtable = subtable[subtable[spp_col] == species]
                if count_col:
                    count = np.sum(spp_subtable[count_col])
                else:
                    count = len(spp_subtable)
                sad_list.append(count)
                    
            result.append((comb, np.array(sad_list)))

        # TODO: Consider changing array to structured with spp name as field

        return spp_list, result


    def parse_criteria(self, criteria):
        '''
        Parses criteria list to get all possible column combinations.

        Parameters
        ----------
        criteria : dict
            (See docstring for EPatch.sad)

        Returns
        -------
        spp_list : ndarray
            1D array listing identifiers for species in the same order as they 
            appear in arrays found in result.
        spp_col : str
            Name of column containing species identifiers.
        count_col : str
            Name of column containing counts, if any.
        combinations : list of dicts
            List of dictionaries giving all possible combinations of criteria. 
            Columns not mentioned in criteria are ignored and will be averaged 
            over in later analyses.

        '''

        spp_list = None
        spp_col = None
        count_col = None
        combinations = []

        # Calculate all possible combinations of columns based on criteria
        # TODO: Add error checking
        for key, value in criteria.items():
            
            # Look for two special values indicating species and count cols
            if value == 'species':
                spp_list = np.unique(self.data_table.table[key])
                spp_col = key
                continue
            if value == 'counts':
                count_col = key
                continue

            # Get levels of categorial or metric data
            if value == 'split':  # Categorial
                levels = np.unique(self.data_table.table[key])
                levels_str = ['==' + str(x) for x in levels]

            else:  # Metric

                # TODO: Make sure divisions equal if count (like LBRI), if true 
                # continuous (like BCI) let divisions be anything. Currently no 
                # checking.
                # TODO: Add support for sampling, if 'sample' in args, randomly 
                # choose n start coords and store in levels

                min = self.data_table.meta[(key, 'minimum')]
                max = self.data_table.meta[(key, 'maximum')]
                prec = self.data_table.meta[(key, 'precision')]

                # TODO: Error if step < prec
                step = (max + prec - min) / value
                starts = np.arange(min, max + prec, step)
                ends = starts + step

                starts_str = ['>=' + str(x) for x in starts]
                ends_str = ['<' + str(x) for x in ends]
                levels_str = zip(starts_str, ends_str)

            # Add these levels to combinations dictionary
            if len(combinations) == 0:  # If first criteria
                for i, level in enumerate(levels_str):
                    combinations.append({key: level})
            else:
                temp_comb = []
                for i, level in enumerate(levels_str):
                    exist_recs = deepcopy(combinations)
                    for rec in exist_recs:
                        rec[key] = level
                    temp_comb += exist_recs
                combinations = temp_comb

        return spp_list, spp_col, count_col, combinations


    def sar(self, div_cols, div_list, criteria, type='sar'):
        '''
        Calulate an empirical species-area relationship given criteria.

        Parameters
        ----------
        div_cols : tuple
            Column names to divide, eg, ('x', 'y'). Must be metric.
        div_list : list of tuples
            List of division pairs in same order as div_cols, eg, [(2,2), 
            (2,4), (4,4)]. Values are number of divisions of div_col.
        criteria : dict
            See docstring for EPatch.sad. Here, criteria SHOULD NOT include 
            items referring to div_cols (if there are any, they are ignored).
        type : string
            'sar' or 'ear' for species or endemics area relationship. EAR is 
            relative to the subtable selected after criteria is applied.

        Returns
        -------
        areas : ndarray
            1D array of areas associated with each element of div_list.
        mean_result : ndarray
            1D array of same length as areas giving mean number of species or 
            endemics at this area.
        full_result : list of ndarrays
            List of same length as areas containing arrays with element for 
            count of species or endemics in each subpatch at corresponding 
            area.
        '''

        # If any element in div_cols in criteria, remove from criteria
        criteria = {k: v for k, v in criteria.items() if k not in div_cols}

        # Loop through div combinations (ie, areas), calc sad, and summarize
        areas = []
        mean_result = []
        full_result = []

        for div in div_list:

            # Add divs to criteria dict
            this_criteria = deepcopy(criteria)
            for i, col in enumerate(div_cols):
                this_criteria[col] = div[i]

            # Get flattened sad for all criteria and this div
            spp_list, sad = self.sad(this_criteria)
            flat_sad = flatten_sad(sad)

            # Store results
            if type == 'sar':
                this_full = np.sum((flat_sad > 0), axis=0)
                this_mean = np.mean(this_full)
            elif type == 'ear':
                totcnt = np.sum(flat_sad, axis=1)
                totcnt_arr = \
                    np.array([list(totcnt),]*np.shape(flat_sad)[1]).transpose()

                this_full = np.sum(np.equal(flat_sad, totcnt_arr), axis=0)
                this_mean = np.mean(this_full)
            else:
                raise NotImplementedError('No SAR of type %s available' % type)

            full_result.append(this_full)
            mean_result.append(this_mean)

            # Store area
            area = 1
            for i, col in enumerate(div_cols):
                min = self.data_table.meta[(col, 'minimum')]
                max = self.data_table.meta[(col, 'maximum')]
                prec = self.data_table.meta[(col, 'precision')]
                length = (max + prec - min)

                area *= length / div[i]

            areas.append(area)

        # Return
        return areas, mean_result, full_result


    def comm(self, div_list):
        '''
        Calculates commonality between pairs of subpatches.
        
        Result is the Sorensen index, equivalent to Chi in Harte et al., also 
        equivalent to 1 - Bray-Curtis dissimilarity.

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


def flatten_sad(sad):
    '''
    Takes a list of tuples, like sad output, ignores keys, and converts values 
    into a 2D array with each value as a column (ie, species in rows, samples 
    in columns.
    '''

    result = np.zeros((len(sad[0][1]), len(sad)))

    for i, tup in enumerate(sad):
        result[:,i] = tup[1]

    return result


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
