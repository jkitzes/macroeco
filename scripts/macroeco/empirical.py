#!/usr/bin/python

'''
Calculating macroecological metrics for empirical or theoretical patch. Patch 
is interpreted broadly as any temporally and spatially defined census.

Classes
-------
- `Patch` -- empirical metrics for census data

Patch Methods
-------------
- `sad` -- calculate species abundance distribution (grid or sample)
- `sar` -- calculate species-area relationship (grid or sample)
- `ear` -- calculate endemics-area relationship (grid or sample)
- `comm` -- calculate commonality between sub-patches (grid)
- `ssad` -- calculate species-level spatial abundance distrib (grid or sample)
- `sed` -- calculate species energy distribution (grid or sample)
- `ied` -- calculate the community energy distribution

- `get_sp_centers` --
- 'get_div_areas' -- return list of areas made by div_list

Misc functions
--------------
- `distance` -- return Euclidean distance between two points
'''

from __future__ import division
import numpy as np
from copy import deepcopy
from data import DataTable

__author__ = "Justin Kitzes"
__copyright__ = "Copyright 2012, Regents of the University of California"
__credits__ = ["John Harte"]
__license__ = None
__version__ = "0.1"
__maintainer__ = "Justin Kitzes"
__email__ = "jkitzes@berkeley.edu"
__status__ = "Development"


class Patch:
    '''
    An object representing an empirical census.

    Parameters
    ----------
    data_path : str
        Path to csv file containing census data.
    subset : dict or str
        Dictionary of permanent subset to data, {'column_name': 'condition'},
        which will limit all analysis to records in which column_name meets the
        condition, ie, {'year': ('==', 2005), 'x': [('>', 20), ('<', 40)]}
        restricts analysis to year 2005 and x values between 20 and 40. These
        conditions can also be passed to the individual methods, but subsetting
        the data table up front may save analysis time.  Subsetting on a string
        would look something like {'name' : [('==', 'John'), ('==', 'Harry')]}.
        In addition, subset can be a query string for a SQL database.

    Attributes
    ----------
    data_table : object of class DataTable
        Object containing patch data and metadata.

    '''

    def __init__(self, datapath, subset = {}):
        '''Initialize object of class Patch. See class documentation.'''
        
        # Handle csv 
        self.data_table = DataTable(datapath, subset=subset)
        
        # If datapath is sql or db the subsetting is already done.
        if type(subset) == type({}):
            self.data_table.table = self.data_table.get_subtable(subset)

    
    def sad(self, criteria, clean=False):
        '''
        Calculates an empirical species abundance distribution given criteria.

        Parameters
        ----------
        criteria : dict
            Dictionary of form {column_name: value}. Must contain a key with a 
            value of 'species' indicating the column with species identifiers 
            (this column must be type categorical in metadata). If a column 
            giving the counts of species found at a point is also in the data, 
            a key with the value 'count' should also be given.

            Value has a different meaning depending on column type:
            - metric - number of divisions of data along this axis, int/float
            - categorical - 'split' calculates each category separately,
              'whole' takes the entire column.
        clean : bool
            If True, all the zeros are removed from the sads.  If False, sads
            are left as is.
        
        Returns
        -------
        result : list 
            List of tuples containing results, where the first element is a
            dictionary of criteria for this calculation and second element is a 
            1D ndarray of length species containing the abundance for each 
            species. The third element is 1D array listing identifiers for 
            species in the same order as they appear in the second element of 
            result. 
        '''
        
        spp_list, spp_col, count_col, engy_col, mass, combinations = \
            self.parse_criteria(criteria)

        if spp_col == None:
            raise TypeError('No species column specified in "criteria" ' +
                                                                   'parameter')
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

            sad_list = np.array(sad_list)

            if clean:
                ind = np.where(sad_list != 0)[0]
                sad_list = sad_list[ind]
                temp_spp_list = spp_list[ind]
            else:
                temp_spp_list = spp_list

                    
            result.append((comb, sad_list, temp_spp_list))

        return result

    def ssad(self, criteria):
        '''
        Calculates empirical species-level spatial abundance distributions
        given criteria.

        Parameters
        ----------
        criteria : dict
            See Patch.sad docstring

        Returns
        -------
        : tuple
            Returns a tuple with two objects.  The first object is an array of
            dicts that correspond to the criteria used to generate each cell.
            The length of the first object in equal to the number of divisions
            specified.  The second object is a dictionary that has length
            species and each keyword is a species.  Each species keyword looks
            up an array with the ssad for the given species.  The array that
            each keyword looks up is the same length as criteria. 


        '''
        sad_return = self.sad(criteria, clean=False)
        spp_list = sad_return[0][2]
        combs, array_res = flatten_sad(sad_return)
        ssad = {}
        
        for i, spp in enumerate(spp_list):
            ssad[spp] = array_res[i,:]

        return combs, ssad

    def parse_criteria(self, criteria):
        '''
        Parses criteria list to get all possible column combinations.

        Parameters
        ----------
        criteria : dict
            (See docstring for Patch.sad)
        energy : bool
            If False, does not return an energy column, if True, returns an
            energy column.
            
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
        engy_col = None
        mass_col = None
        combinations = []

        # Calculate all possible combinations of columns based on criteria
        # TODO: Add error checking
        for key, value in criteria.items():
            
            # Look for two special values indicating species and count cols
            if value == 'species':
                spp_list = np.unique(self.data_table.table[key])
                spp_col = key
                continue
            if value == 'count':
                count_col = key
                continue
            if value == 'energy':
                engy_col = key
                continue
            if value == 'mass':
                mass_col = key
                continue

            # Get levels of categorial or metric data
            if value == 'split':  # Categorial
                levels = np.unique(self.data_table.table[key])
                levels_str = [('==' , x.astype(levels.dtype)) for x in levels]
            elif value == 'whole':
                # Random string to minimize chance of overlap?
                levels_str = [('==','whole')]
            else:  # Metric

                # TODO: Throw a warning if the data is not divisible by the
                # divisions specified.
                try:
                    dmin = self.data_table.meta[(key, 'minimum')]
                    dmax = self.data_table.meta[(key, 'maximum')]
                    dprec = self.data_table.meta[(key, 'precision')]

                    # TODO: Error if step < prec
                    step = (dmax + dprec - dmin) / value
                    starts = np.arange(dmin, dmax + dprec, step)
                    ends = starts + step
                except TypeError:
                    raise TypeError('Unable to proceed to with values ' +
                                    'obtained from metadata.  Please check ' + 
                                    'the metadata file')
                

                starts_str = [('>=', x) for x in starts]
                ends_str = [('<', x) for x in ends]
                levels_str = [list(lvl) for lvl in zip(starts_str, ends_str)]


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
        
        if len(combinations) == 0:
            combinations.append({})
        
        return spp_list, spp_col, count_col, engy_col, mass_col, combinations



    def sar(self, div_cols, div_list, criteria, form='sar'):
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
        form : string
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
            sad_return = self.sad(this_criteria)
            flat_sad = flatten_sad(sad_return)[1]

            # Store results
            if form == 'sar':
                this_full = np.sum((flat_sad > 0), axis=0)
                this_mean = np.mean(this_full)
            elif form == 'ear':
                totcnt = np.sum(flat_sad, axis=1)
                totcnt_arr = \
                    np.array([list(totcnt),]*np.shape(flat_sad)[1]).transpose()

                this_full = np.sum(np.equal(flat_sad, totcnt_arr), axis=0)
                this_mean = np.mean(this_full)
            else:
                raise NotImplementedError('No SAR of form %s available' % form)

            full_result.append(this_full)
            mean_result.append(this_mean)

            # Store area
            area = 1
            for i, col in enumerate(div_cols):
                dmin = self.data_table.meta[(col, 'minimum')]
                dmax = self.data_table.meta[(col, 'maximum')]
                dprec = self.data_table.meta[(col, 'precision')]
                length = (dmax + dprec - dmin)

                area *= length / div[i]

            areas.append(area)

        # Return
        rec_sar = np.array(zip(mean_result, areas), dtype=[('items', np.float),
                                                           ('area', np.float)])
        return rec_sar, full_result

    def ied(self, criteria, normalize=True, exponent=0.75):
        '''
        Calculates the individual energy distribution for the entire community
        given the criteria

        Parameters
        ----------
        criteria : dict
            Dictionary must have contain a key with the value 'energy'.  See
            sad method for further requirements.
        normalize : bool
            If True, this distribution is normalized by dividing by the lowest
            energy value within each element of criteria. If False, returns raw
            energy values.
        exponent : float
            The exponent of the allometric scaling relationship if energy is
            calculated from mass.

        Returns
        -------
        result : list
            List of tuples containing results, where first element is 
            dictionary of criteria for this calculation and second element is a 
            1D ndarray containing the energy measurement of each individual in
            the subset.  The third element is the full (not unique) species
            list for the given criteria. 

        Notes
        -----
        If count_col is None or is all ones, the entire energy column for each
        subtable is returned.  Else, the average energy per individual,
        repeated for each individual is returned. This is equivalent to the psi
        distribution from Harte (2011).


        '''
        
        spp_list, spp_col, count_col, engy_col, mass_col, combinations = \
            self.parse_criteria(criteria)

        if engy_col == None and mass_col == None:
            raise ValueError("No energy or mass column given")
        elif engy_col == None and mass_col != None:
            mass = True
            this_engy = mass_col
        else:
            mass = False
            this_engy = engy_col

        result = []
        for comb in combinations:

            subtable = self.data_table.get_subtable(comb)
            
            # If all counts are not 1
            if count_col and (not np.all(subtable[count_col] == 1)):
                
                # Remove any zero counts
                subtable = subtable[subtable[count_col] != 0]
                # Convert counts to ints
                temp_counts = subtable[count_col].astype(int)

                energy = np.repeat((subtable[this_engy] /
                        subtable[count_col]), temp_counts)
                species = np.repeat(subtable[spp_col], temp_counts)
            else:
                energy = subtable[this_engy] 
                species = subtable[spp_col]

            # Convert mass to energy if mass is True
            if mass:
                energy = (energy ** exponent)
                
            # Normalizing energy
            if normalize:
                energy = energy / np.min(energy)
            result.append((comb, energy, species))

        return result

    def sed(self, criteria, normalize=True, exponent=0.75):
        '''
        Calculates the species-level energy distribution for each given species
        in the community.

        Parameters
        ----------
        criteria : dict
            Dictionary must have contain a key with the value 'energy' or
            'mass'.  See sad method for further requirements.

        Returns
        -------
        result : list of tuples
            Each tuple contains two objects.  The first object is a dict with
            the division specifications that generated the given species energy
            distributions.  The second object is a dict with a keyword
            corresponding to each species in the spp_list.  Each species
            keyword looks up a np.array that contains the given species
            energy distribution.

        Note
        ----
        This is equivalent to the theta distribution from Harte (2011).

        '''
        spp_list, spp_col, count_col, engy_col, mass_col, combinations = \
            self.parse_criteria(criteria)

        ied = self.ied(criteria, normalize=normalize, exponent=exponent)

        result = []
        for this_ied in ied:
            this_criteria_sed = {}

            for spp in spp_list:
                spp_ind = (spp == this_ied[2])
                this_spp_sed = this_ied[1][spp_ind]
                this_criteria_sed[spp] = this_spp_sed

            result.append((this_ied[0], this_criteria_sed))
        
        return result
    
    def ased(self, criteria, normalize=True, exponent=0.75):
        '''
        Calculates the average species energy distribution for each given
        species in a subset. 
        
        Parameters
        ----------
        criteria : dict
            Dictionary must have contain a key with the value 'energy' or
            'mass'.  See sad method for further requirements.
        
        Returns
        -------
        result : list 
            List of tuples containing results, where the first element is a
            dictionary of criteria for this calculation and second element is a 
            1D ndarray of length species containing the average energy for each 
            species. The third element is 1D array listing identifiers for 
            species in the same order as they appear in the second element of 
            result.         

        Notes
        -----
        This is equivalent to the nu distribution from Harte 2011

        '''

        sed = self.sed(criteria, normalize=normalize, exponent=exponent)

        result = []
        for this_sed in sed:
            spp_list = list(this_sed[1].viewkeys())
            spp_list.sort()

            # Take the mean energy for each species
            nu = [np.mean(this_sed[1][spp]) for spp in spp_list if
                                                    len(this_sed[1][spp]) != 0]
            # Truncated spp_list if necessary
            spp_list = [spp for spp in spp_list if len(this_sed[1][spp]) != 0]
            
            result.append((this_sed[0], np.array(nu), np.array(spp_list)))

        return result

def flatten_sad(sad):
    '''
    Takes a list of tuples, like sad output, ignores keys, and converts values 
    into a 2D array with each value as a column (ie, species in rows, samples 
    in columns.
    '''

    combs = [cmb[0] for cmb in sad]
    result = np.zeros((len(sad[0][1]), len(sad)))

    for i, tup in enumerate(sad):
        result[:,i] = tup[1]

    return combs, result


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
