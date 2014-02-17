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
- `universal_sar` -- calculates the universal sar curve
- `ear` -- calculate endemics-area relationship (grid or sample)
- `comm` -- calculate commonality between sub-patches (grid)
- `ssad` -- calculate species-level spatial abundance distrib (grid or sample)
- `sed` -- calculate species energy distribution (grid or sample)
- `ied` -- calculate the community (individual) energy distribution
- `ased` -- calculate the average species energy distribution
- `tsed` -- calculate the total species energy distribution

- `get_sp_centers` --
- 'get_div_areas' -- return list of areas made by div_list

Misc functions
--------------
- `distance` -- return Euclidean distance between two points
'''

from __future__ import division
import numpy as np
from math import radians, cos, sin, asin, sqrt
import itertools
from copy import deepcopy
from data import DataTable
import scipy.spatial.distance as dist
import shapely.geometry as geo


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
            appear in arrays found in result. Contains all species in table,
            although all species may not appear in subtables that are defined
            by combinations.
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
                                    'the metadata file and/or parameters file')


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


    def sar(self, div_cols, div_list, criteria, form='sar', output_N=False):
        '''
        Calculate an empirical species-area relationship given criteria.

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
        output_N : bool
            Adds the column N to the output rec array which contains the
            average N for a given area.

        Returns
        -------
        rec_sar: structured array
            Returns a structured array with fields 'items' and 'area' that
            contains the average items/species for each given area specified by
            critieria.
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
        N_result = []

        for div in div_list:

            # Add divs to criteria dict
            this_criteria = deepcopy(criteria)
            for i, col in enumerate(div_cols):
                this_criteria[col] = div[i]

            # Get flattened sad for all criteria and this div
            sad_return = self.sad(this_criteria)

            if output_N:
                N_result.append(np.mean([sum(sad[1]) for sad in sad_return]))

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
        if not output_N:
            rec_sar = np.array(zip(mean_result, areas), dtype=[('items',
                                                np.float), ('area', np.float)])
        else:
            rec_sar = np.array(zip(mean_result, N_result, areas),
              dtype=[('items', np.float), ('N', np.float), ('area', np.float)])

        return rec_sar, full_result


    def universal_sar(self, div_cols, div_list, criteria, include_full=False):
        '''
        Calculates the empirical universal sar given criteria. The universal
        sar calculates the slope of the SAR and the ratio of N / S at all
        the areas in div_cols (where N is the total number of species and S is
        the total number of species).

        This function assumes that the div_list contains halvings.  If they are
        not, the function will still work but the results will be meaningless.
        An example a of div_list with halvings is:

        [(1,1), (1,2), (2,2), (2,4), (4,4)]

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
        include_full : bool
            If include_full = True, the division (1,1) will be included if it
            was now already included. Else it will not be included.  (1,1) is
            equivalent to the full plot


        Returns
        -------
        z_array : a structured array
            Has the columns names:
            'z' : slope of the SAR at the given area
            'S' : Number of species at the given division
            'N' : Number of individuals at the given division
            'N/S' : The ratio of N/S at the given division


        Notes
        -----
        If you give it n divisions in div_list you will get a structured array
        back that has length n - 2.  Therefore, if you only have one
        '''

        # If (1,1) is not included, include it
        if include_full:
            try:
                div_list.index((1,1))
            except ValueError:
                div_list.insert(0, (1,1))

        # Run sar with the div_cols
        sar = self.sar(div_cols, div_list, criteria, output_N=True)[0]

        # sort by area
        sar = np.sort(sar, order=['area'])[::-1]

        # Calculate z's
        if len(sar) >= 3: # Check the length of sar
            z_list = [z(sar['items'][i - 1], sar['items'][i + 1]) for i in
                 np.arange(1, len(sar)) if sar['items'][i] != sar['items'][-1]]
        else:
            return np.empty(0, dtype=[('z', np.float), ('S', np.float), ('N',
                                                 np.float), ('N/S', np.float)])

        N_over_S = sar['N'][1:len(sar) - 1] / sar['items'][1:len(sar) - 1]

        z_array = np.array(zip(z_list, sar['items'][1:len(sar) - 1],
            sar['N'][1:len(sar) - 1], N_over_S), dtype=[('z', np.float), ('S',
            np.float), ('N', np.float), ('N/S', np.float)])

        return z_array

    def comm_sep(self, plot_locs, criteria, loc_unit=None):
        '''
        Calculates commonality (Sorensen and Jaccard) between pairs of plots.

        Parameters
        ----------
        plot_locs : dict
            Dictionary with keys equal to each plot name, which must be
            represented by a column in the data table, and values equal to a
            tuple of the x and y coordinate of each plot
        criteria : dict
            See docstring for Patch.sad.
        loc_unit : str
            Unit of plot locations. Special cases include 'decdeg' (decimal
            degrees), returns result in km. Otherwise ignored.

        Returns
        -------
        result: structured array
            Returns a structured array with fields plot-a and plot-b (names of
            two plots), dist (distance between plots), and sorensen and jaccard
            (similarity indices). Has row for each unique pair of plots.
        '''

        # Set up sad_dict with key=plot and val=clean sad for that plot
        sad_dict = {}

        # Loop through all plot cols, updating criteria, and getting spp_list
        for plot in plot_locs.keys():

            # Find current count col and remove it from criteria
            for crit_key in criteria.keys():
                if criteria[crit_key] == 'count':
                    criteria.pop(crit_key, None)

            # Add this plot as col with counts
            criteria[plot] = 'count'

            # Get SAD for existing criteria with this plot as count col
            sad_return = self.sad(criteria, clean=True)

            # Check that sad_return only has one element, or throw error
            if len(sad_return) > 1:
                raise NotImplementedError('Too many criteria for comm_sep')

            # Get unique species list for this plot and store in sad_dict
            sad_dict[plot] = sad_return[0][2]

        # Set up recarray to hold Sorensen index for all pairs of plots
        n_pairs = np.sum(np.arange(len(plot_locs.keys())))
        result = np.recarray((n_pairs,), dtype=[('plot-a','S32'),
                                                ('plot-b', 'S32'),
                                                ('spp-a', int),
                                                ('spp-b', int),
                                                ('dist', float),
                                                ('sorensen', float),
                                                ('jaccard', float)])

        # Loop through all combinations of plots and fill in result table
        row = 0
        for pair in itertools.combinations(plot_locs.keys(), 2):

            # Names of plots
            plota = pair[0]
            plotb = pair[1]

            result[row]['plot-a'] = plota
            result[row]['plot-b'] = plotb

            # Calculate inter-plot distance
            if loc_unit == 'decdeg':
                result[row]['dist'] = decdeg_distance(plot_locs[plota],
                                                      plot_locs[plotb])
            else:
                result[row]['dist'] = distance(plot_locs[plota],
                                               plot_locs[plotb])

            # Get similarity indices
            spp_a = len(sad_dict[plota])
            spp_b = len(sad_dict[plotb])

            result[row]['spp-a'] = spp_a
            result[row]['spp-b'] = spp_b

            intersect = set(sad_dict[plota]).intersection(sad_dict[plotb])
            union = set(sad_dict[plota]).union(sad_dict[plotb])

            # Fill in zero if denom is zero
            if spp_a + spp_b == 0:
                result[row]['sorensen'] = 0
            else:
                result[row]['sorensen'] = (2*len(intersect)) / (spp_a+spp_b)

            if len(union) == 0:
                result[row]['jaccard'] = 0
            else:
                result[row]['jaccard'] = len(intersect) / len(union)

            # Increment row counter
            row += 1

        return result

    def o_ring(self, div_cols, bin_edges, criteria, n0_min_max=None,
                  edge_correct=False, density=False):
        '''
        Calculates univariate O-ring for a species.

        Parameters
        ----------
        div_cols : tuple
            Column names containing x and y coordinates of individuals
        bin_edges : iterable
            List of edges of distance classes to bin histogram of distances
        criteria : dict
            See docstring for Patch.sad. Count column must be used.
        n0_min_max : tuple
            Optional min and max abundance for species to consider. Useful for
            ignoring rare species with few samples and abundant species for
            which calculation would take a long time.
        edge_correct : bool
            Correct histograms by replacing count of individuals at distance
            bin with expected count if entire ring at that distance was
            available (part of ring may fall outside of plot). Default False.
        density : bool
            If True, return densities (counts divided by area of torus defined
            by bin edges) instead of counts. Default False.

        Returns
        -------
        result : tuple
            List of tuples with three elements each. First is combination used
            to generate results, second is spp_list for that combination
            (includes all species in entire landscape), and third is list of
            length spp_list giving histogram of pairwise distances for each
            species.

        Notes
        -----
        Pairwise distances are directional, giving n(n-1) total distances, as
        edge correction is directional.

        If there are no records in a combination, histogram will be None. If
        there are records but a species has only one individual, histogram
        will be all zeros.

        When using edge_correct or density, the maximum distance used for edge
        correction, given by the mean of the last two bin_edge values, should
        be no greater than one half the longer dimension of the plot. This
        ensures that it is not possible for an entire edge correction buffer
        to be outside of the plot, which could lead to divide by zero errors.

        '''

        spp_list, spp_col, count_col, engy_col, mass, combinations = \
            self.parse_criteria(criteria)

        bin_edges = np.array(bin_edges)

        result_list = []

        for comb in combinations:

            # If comb includes division, cannot also use edge correction
            # This would require better parsing of plot boundaries for division
            if (not comb.keys() == []) and edge_correct:
                raise NotImplementedError("Edge correction cannot be used "
                                          "with combinations.")

            # Get appropriate subtable for this combination
            subtable = self.data_table.get_subtable(comb)

            # Declare empty list for all histograms for all species
            spp_hist_list = []

            # If density is True, set edge_correct to True
            if density:
                edge_correct = True

            # Set up plot polygon for edge correction
            if edge_correct:
                xmin = self.data_table.meta[(div_cols[0], 'minimum')]
                xmax = self.data_table.meta[(div_cols[0], 'maximum')]
                ymin = self.data_table.meta[(div_cols[1], 'minimum')]
                ymax = self.data_table.meta[(div_cols[1], 'maximum')]

                plot = geo.box(xmin, ymin, xmax, ymax)

                all_r = (bin_edges[:-1] + bin_edges[1:]) / 2

            # Calculate areas of all toruses
            if density:
                ring_areas = []
                for i in range(len(bin_edges) - 1):
                    ring_areas.append(np.pi*(bin_edges[i+1]**2 -
                                             bin_edges[i]**2))
                ring_areas = np.array(ring_areas)

            # Loop all species
            for spp in spp_list:

                spp_subtable = subtable[subtable[spp_col] == spp]

                # If spp not present or singleton, continue
                # Ensure that if single record but count > 1, do analysis
                if len(spp_subtable) == 0:
                    spp_hist_list.append(None)
                    continue

                # Get n0, accounting for count col
                if count_col:
                    count = np.sum(spp_subtable[count_col])
                else:
                    count = len(spp_subtable)

                # Skip this spp if there is a min_max set and n0 out of range
                if n0_min_max and (count < n0_min_max[0] or count >
                                   n0_min_max[1]):
                    spp_hist_list.append(None)
                    continue

                # Get list of all points and all counts
                x = spp_subtable[div_cols[0]]
                y = spp_subtable[div_cols[1]]
                all_points = zip(x,y)
                all_counts = list(spp_subtable[count_col])

                # Declare array to hold histogram of pairwise distances
                all_hist = np.zeros(len(bin_edges) - 1)

                # Declare array to hold all sampled areas per bin
                if density:
                    all_areas = np.zeros(len(ring_areas))

                # Go through all_points
                for i, this_point in enumerate(all_points):

                    # Get this point and remove from list of all points
                    this_count = all_counts[i]

                    # Create list of all other points and counts except this
                    all_other_points = all_points[0:i] + all_points[i+1:]
                    all_other_counts = all_counts[0:i] + all_counts[i+1:]

                    # Get dist from this point to all other points
                    # If no other points, other_dist is empty
                    # May still be other individs at this point
                    if all_other_points:
                        other_dist = dist.cdist(np.array([this_point]),
                                                np.array(all_other_points))
                    else:
                        other_dist = np.array(())

                    # Repeat other point distances to acccount for their counts
                    other_dist = np.repeat(other_dist, all_other_counts)

                    # Repeat entire other_dist to account for count here
                    other_dist = np.tile(other_dist, this_count)

                    # Add 0 distances between individs at this point
                    # Multiplied by two to get directional pairwise dists
                    n_this_dists = this_count - 1
                    if n_this_dists > 0:
                        other_dist = np.concatenate((other_dist,
                                                np.zeros(n_this_dists*2)))

                    # Calculate histogram of distances to other points
                    hist, _ = np.histogram(other_dist, bin_edges)

                    # Edge correct distance
                    if edge_correct:
                        corr_fact = np.zeros(len(all_r))
                        for i, r in enumerate(all_r):
                            x, y = this_point
                            circ = geo.Point(x,y).buffer(r,resolution=64)
                            out_len = circ.boundary.difference(plot).length
                            in_frac = ((circ.boundary.length - out_len) /
                                       circ.boundary.length)
                            corr_fact[i] = in_frac
                        hist = hist / corr_fact

                    # Store sampled area at each dist for density calculation
                    if density:
                        all_areas += (ring_areas * corr_fact)

                    # Add this point results to main histogram
                    all_hist += hist

                # If density, divide all values by summed sampled torus areas
                if density:
                    all_hist = all_hist / all_areas

                # Append final hist for this species to running list
                spp_hist_list.append(all_hist)

            # For this comb, create and append tuple to result list
            result_list.append((comb, spp_list, spp_hist_list))

        return result_list


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

    def sed(self, criteria, normalize=True, exponent=0.75, clean=False):
        '''
        Calculates the species-level energy distribution for each given species
        in the community.

        Parameters
        ----------
        criteria : dict
            Dictionary must have contain a key with the value 'energy' or
            'mass'.  See sad method for further requirements.
        normalize : bool
            If True, this distribution is normalized by dividing by the lowest
            energy value within each element of criteria. If False, returns raw
            energy values.
        exponent : float
            The exponent of the allometric scaling relationship if energy is
            calculated from mass
        clean : bool
            If False, sed dictionary contains all species.  If True, species
            with no individuals are removed.  This is useful when subsetting.

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
        The theta distribution from Harte (2011) is a an sed.

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

                if clean: # If True, don't add empty species lists
                    if len(this_spp_sed) > 0:
                        this_criteria_sed[spp] = this_spp_sed
                else:
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

    def tsed(self, criteria, normalize=True, exponent=0.75):
        '''
        Calculates the total species energy distribution for each given
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

        '''

        sed = self.sed(criteria, normalize=normalize, exponent=exponent)

        result = []
        for this_sed in sed:
            spp_list = list(this_sed[1].viewkeys())
            spp_list.sort()

            # Take the mean energy for each species
            omega = [np.sum(this_sed[1][spp]) for spp in spp_list if
                                                    len(this_sed[1][spp]) != 0]
            # Truncated spp_list if necessary
            spp_list = [spp for spp in spp_list if len(this_sed[1][spp]) != 0]
            
            result.append((this_sed[0], np.array(omega), np.array(spp_list)))

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


def decdeg_distance(pt1, pt2):
    ''' Calculate Earth surface distance (in km) between decimal latlong points
    using Haversine approximation.

    http://stackoverflow.com/questions/15736995/how-can-i-quickly-estimate-the-distance-between-two-latitude-longitude-points
    '''
    lat1, lon1 = pt1
    lat2, lon2 = pt2

    # Convert decimal degrees to radians
    lon1, lat1, lon2, lat2 = map(radians, [lon1, lat1, lon2, lat2])

    # haversine formula
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    a = sin(dlat/2)**2 + cos(lat1) * cos(lat2) * sin(dlon/2)**2
    c = 2 * asin(sqrt(a))
    km = 6367 * c

    return km

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

def z(doubleS, halfS):
    '''Calculates the z for a double S value and a half S value'''

    return np.log(doubleS / halfS) / (2 * np.log(2))
