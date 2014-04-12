from __future__ import division
import os
import re
import copy
from configparser import ConfigParser
import itertools
from copy import deepcopy
from twiggy import log
log = log.name('emp ')

import numpy as np
import pandas as pd
import scipy.spatial.distance as dist
try:
    import shapely.geometry as geo
except:
    pass
# TODO: Make shapely import work with pyinstaller

from ..misc import doc_sub, log_start_end

metric_params = \
    """patch : Patch obj
        Patch object containing data for analysis
    cols : str
        Indicates which column names in patch data table are associated with
        species identifiers, counts, energy, and mass. See Notes.
    splits : str
        If multiple analyses for subsets of patch data table are desired,
        specifies how columns should be split. See Notes."""

metric_return = \
    """list
        List of tuples containing results, where the first element of each
        tuple is a string indicating the split values used for that result and
        second element is a dataframe giving the result."""

cols_note = \
    """The parameter ``cols`` is a dictionary with keys for four special
    columns and values giving the column name in the patch data table
    associated with each special column.

    - spp_col - Unique species identifiers
    - count_col - Number of individuals at a location
    - energy_col - Energy of individuals
    - mass_cal - Mass of individuals

    Only spp_col is always mandatory. Note that the value of spp_col may be
    set to a columm in the data table giving the genus, family, functional
    group, etc., which allows for analysis of this metric by those groups.
    count_col is used when multiple individuals of a species may be found at
    a single recorded location, as is the case in gridded censuses where all
    individuals in a quadrat are "assigned" to a single point. energy_col
    and mass_col are used for energy-based metrics."""

splits_note = \
    """The parameter ``splits`` is a semicolon-separated string in the form of
    "column: value", where column is a name of a column in the patch data
    table and value is either (a) an integer giving the number of
    equally-spaced divisions of a column, or (b) the special keyword
    'split', which evaluates all unique levels of a column.

    For example, presume a data table has columns for x and y spatial
    coordinates and a column for year, of which there are three. The string
    "x:2; y:2; year:split" will perform the analysis separately for each of
    four subplots of the patch (created by dividing the x and y coordinates
    each into two equally sized divisions) within each of the three years,
    for a total of 12 separate analyses."""

division_note = \
    """The parameter divisions describes how to successively divide the patch
    along the x_col and y_col dimensions. For
    example, the string '1,2; 2,2; 2,4' will produce an output table with three
    rows, giving the result across two subplots when the patch is split
    along y_col, across four subplots when the patch is split into a 2x2 grid,
    and across eight subplots when the patch is split into 2 parts along x_col
    and 4 parts along y_col."""


class Patch(object):
    """
    An object representing an empirical census

    Parameters
    ----------
    metadata_path : str
        Path to metadata file describing census data
    subset : str
        String describing subset of data to use for Patch analysis. See Notes.

    Attributes
    ----------
    table : dataframe
        Table of census data recorded in patch
    meta : ConfigParser obj
        Object similar to dict describing data table, loaded from metadata
        file at metadata_path and processed by subset
    subset : str
        Subset string passed as parameter

    Notes
    -----
    The table file described by the metadata must contain column names
    consisting only of letters and numbers, with no spaces or other special
    characters.

    The parameter subset takes different forms depending on whether the data
    file described by the metadata is a csv or a sql/db file.

    For csv data files, subset is a semicolon-separated string describing
    subset operations. For example, the string "year==2005; x>20; x<40;
    spp=='cabr'" loads a data table containing only records for which the year
    is 2005, x values are between 20 and 40, and species is 'cabr'. Note that
    for categorical columns, the value of the column must be enclosed in single
    quotes.

    For sql/db files, subset is a SQL query string that selects the data from
    the data file.

    The meta attribute of this object is processed to reflect the value of
    subset. If columns with a min and a max are included in the subset string,
    the min and max values for that column in meta will be updated to reflect
    the specified limits.

    An empty Patch object can be created with a metadata_path of None.

    """

    def __init__(self, metadata_path, subset=''):

        if not metadata_path:  # Allow for creation of empty patch
            self.meta = None
            self.subset = ''
            self.table = None
        else:
            self.meta = ConfigParser()
            self.meta.read(metadata_path)
            self.subset = subset
            self.table = self._load_table(metadata_path,
                                          self.meta['Description']['datapath'])


    def _load_table(self, metadata_path, data_path):
        """
        Load data table, taking subset if needed

        Parameters
        ----------
        metadata_path : str
            Path to metadata file
        data_path : str
            Path to data file, absolute or relative to metadata file

        Returns
        -------
        dataframe
            Table for analysis

        """

        metadata_dir = os.path.dirname(metadata_path)
        data_path = os.path.normpath(os.path.join(metadata_dir, data_path))

        extension = data_path.split('.')[-1]

        if extension == 'csv':
            full_table = pd.read_csv(data_path, index_col=False)
            table = _subset_table(full_table, self.subset)
            self.meta = _subset_meta(self.meta, self.subset)
        elif extension in ['db', 'sql']:
            table = self._get_db_table(data_path, extension)
        else:
            raise TypeError('Cannot process file of type %s' % extension)

        return table

    def _get_db_table(self, data_path, extension):
        """
        Query a database and return query result as a recarray

        Parameters
        ----------
        data_path : str
            Path to the database file
        extension : str
            Type of database, either sql or db

        Returns
        -------
        table : recarray
            The database query as a recarray

        """
        # TODO: This is probably broken
        raise NotImplementedError, "SQL and db file formats not yet supported"

        # Load table
        if extension == 'sql':
            con = lite.connect(':memory:')
            con.row_factory = lite.Row
            cur = con.cursor()

            with open(data_path, 'r') as f:
                sql = f.read()

            cur.executescript(sql)

        else:
            con = lite.connect(data_path)
            con.row_factory = lite.Row
            cur = con.cursor()

        cur.execute(self.subset)

        # Check that table is not empty
        db_info = cur.fetchall()
        try:
            col_names = db_info[0].keys()
        except IndexError:
            raise lite.OperationalError("Query %s to database %s is empty" %
                                        (query_str, data_path))

        # Convert objects to tuples
        converted_info = [tuple(x) for x in db_info]

        # NOTE: Using default value for Unicode: Seems better than checking
        # lengths.  Should we keep the type as unicode?
        dtypes=[type(x) if type(x) != unicode else 'S150' for x in db_info[0]]

        table = np.array(converted_info, dtype=zip(col_names, dtypes))
        con.commit()
        con.close()

        # Return a recarray for consistency
        # TODO: This should now be a pd.dataframe
        return table.view(np.recarray)


def _subset_table(full_table, subset):
    """
    Return subtable matching all conditions in subset

    Parameters
    ----------
    full_table : dataframe
        Entire data table
    subset : str
        String describing subset of data to use for analysis

    Returns
    -------
    dataframe
        Subtable with records from table meeting requirements in subset

    """
    if not subset:
        return full_table

    # TODO: Figure out syntax for logical or
    conditions = subset.replace(' ','').split(';')

    valid = np.ones(len(full_table), dtype=bool)
    for condition in conditions:
        this_valid = eval('full_table.' + condition)
        valid = np.logical_and(valid, this_valid)

    return full_table[valid]

def _subset_meta(full_meta, subset):
    """
    Return metadata reflecting all conditions in subset

    Parameters
    ----------
    full_meta : ConfigParser obj
        Metadata object
    subset : str
        String describing subset of data to use for analysis

    Returns
    -------
    Configparser object or dict
        Updated version of full_meta accounting for subset string

    """
    if not subset:
        return full_meta

    meta = {}  # Make deepcopy of entire meta (all section dicts in meta dict)
    for key, val in full_meta.iteritems():
        meta[key] = copy.deepcopy(dict(val))

    conditions = subset.replace(' ','').split(';')

    for condition in conditions:
        condition_list = re.split('[<>=]', condition)
        col = condition_list[0]
        val = condition_list[-1]

        try:
            col_step = meta[col]['step']
        except:  # If there's no metadata for this col, do nothing
            continue

        operator = re.sub('[^<>=]', '', condition)

        if operator == '==':
            meta[col]['min'] = val
            meta[col]['max'] = val
        elif operator == '>=':
            meta[col]['min'] = val
        elif operator == '>':
            meta[col]['min'] = str(eval(val) + eval(col_step))
        elif operator == '<=':
            meta[col]['max'] = val
        elif operator == '<':
            meta[col]['max'] = str(eval(val) - eval(col_step))
        else:
            raise ValueError, "Subset %s not valid" % condition

    return meta


@log_start_end
@doc_sub(metric_params, metric_return, cols_note, splits_note)
def sad(patch, cols, splits, clean=True):
    """
    Calculates an empirical species abundance distribution

    Parameters
    ----------
    {0}
    clean : bool
        If True, all species with zero abundance are removed from SAD results.
        Default False.

    Returns
    -------
    {1} Result has two columns: spp (species identifier) and y (individuals of
    that species).

    Notes
    -----
    {2}

    {3}

    """

    (spp_col, count_col), patch = \
        _get_cols(['spp_col', 'count_col'], cols, patch)

    full_spp_list = np.unique(patch.table[spp_col])

    # Loop through each split
    result_list = []
    for substring, subpatch in _yield_subpatches(patch, splits):

        # Get abundance for each species
        sad_list = []
        for spp in full_spp_list:
            this_spp = (subpatch.table[spp_col] == spp)
            count = np.sum(subpatch.table[count_col][this_spp])
            sad_list.append(count)

        # Create dataframe of spp names and abundances
        subdf = pd.DataFrame({'spp': full_spp_list, 'y': sad_list})

        # Remove zero abundance rows if requested
        if clean:
            subdf = subdf[subdf['y'] > 0]

        # Append subset result
        result_list.append((substring, subdf))

    # Return all results
    return result_list


@log_start_end
@doc_sub(metric_params, metric_return, cols_note, splits_note)
def ssad(patch, cols, splits):
    """
    Calculates an empirical intra-specific spatial abundance distribution

    Parameters
    ----------
    {0}

    Returns
    -------
    {1} Result has one column giving the individuals of species in each
    subplot.

    Notes
    -----
    {2}

    {3}

    """

    # Get and check SAD
    sad_results = sad(patch, cols, splits, clean=False)

    # Create dataframe with col for spp name and numbered col for each split
    for i, sad_result in enumerate(sad_results):
        if i == 0:  # For first result, create dataframe
            fulldf = sad_result[1]
            fulldf.columns = ['spp', '0']  # Renames y col to 0
        else:  # For other results, append col to dataframe, named by num
            fulldf[str(i)] = sad_result[1]['y']

    # Get each spp SSAD (row of fulldf) and append as tuple in result_list
    result_list = []
    for _, row in fulldf.iterrows():
        row_values_array = np.array(row[1:], dtype=float)
        result_list.append((row[0], pd.DataFrame({'y': row_values_array})))

    # Return all results
    return result_list


@log_start_end
@doc_sub(metric_params, metric_return, cols_note, splits_note, division_note)
def sar(patch, cols, splits, divs, ear=False):
    """
    Calculates an empirical species area or endemics area relationship

    Parameters
    ----------
    {0}
    divs : str
        Description of how to divide x_col and y_col. See notes.
    ear : bool
        If True, calculates an endemics area relationship

    Returns
    -------
    {1} Result has three columns, div, x, and y, that give the ID for the
    division given as an argument, fractional area, and the mean species
    richness at that division.

    Notes
    -----
    {2}

    For the SAR and EAR, cols must also contain x_col and y_col, giving the x
    and y dimensions along which to grid the patch.

    {3}

    {4}

    """

    def sar_y_func(spatial_table, all_spp):
        return np.mean(spatial_table['n_spp'])

    def ear_y_func(spatial_table, all_spp):
        endemic_counter = 0
        for spp in all_spp:
            spp_in_cell = [spp in x for x in spatial_table['spp_set']]
            spp_n_cells = np.sum(spp_in_cell)
            if spp_n_cells == 1:  # If a spp is in only 1 cell, endemic
                endemic_counter += 1
        n_cells = len(spatial_table)
        return endemic_counter / n_cells # mean endemics / cell

    if ear:
        y_func = ear_y_func
    else:
        y_func = sar_y_func

    return _sar_ear_inner(patch, cols, splits, divs, y_func)


def _sar_ear_inner(patch, cols, splits, divs, y_func):
    """
    y_func is function calculating the mean number of species or endemics,
    respectively, for the SAR or EAR
    """

    (spp_col, count_col, x_col, y_col), patch = \
        _get_cols(['spp_col', 'count_col', 'x_col', 'y_col'], cols, patch)

    # Loop through each split
    result_list = []
    for substring, subpatch in _yield_subpatches(patch, splits):

        # Get A0
        A0 = _patch_area(subpatch, x_col, y_col)

        # Loop through all divisions within this split
        all_spp = np.unique(subpatch.table[spp_col])
        subresultx = []
        subresulty = []
        subresultnspp = []
        subresultnindivids = []
        subdivlist = _split_divs(divs)
        for subdiv in subdivlist:
            spatial_table = _yield_spatial_table(subpatch, subdiv, spp_col,
                                                 count_col, x_col, y_col)
            subresulty.append(y_func(spatial_table, all_spp))
            subresultx.append(A0 / eval(subdiv.replace(',', '*')))
            subresultnspp.append(np.mean(spatial_table['n_spp']))
            subresultnindivids.append(np.mean(spatial_table['n_individs']))

        # Append subset result
        subresult = pd.DataFrame({'div': subdivlist, 'x': subresultx,
                                  'y': subresulty, 'n_spp': subresultnspp,
                                  'n_individs': subresultnindivids})
        result_list.append((substring, subresult))

    return result_list


def _split_divs(divs):
    if type(divs) == type((1,1)):  # Tuple (occurs when main evals single div)
        subdivlist = [str(divs)[1:-1]]
    else: # String
        subdivlist = divs.split(';')
    return subdivlist


@log_start_end
@doc_sub(metric_params, metric_return, cols_note, splits_note)
def comm_grid(patch, cols, splits, divs, metric='Sorensen'):
    """
    Calculates commonality as a function of distance for a gridded patch

    Parameters
    ----------
    {0}
    divs : str
        Description of how to divide x_col and y_col. Unlike SAR and EAR, only
        one division can be given at a time. See notes.
    metric : str
        One of Sorensen or Jaccard, giving the metric to use for commonality
        calculation

    Returns
    -------
    {1} Result has three columns, pair, x, and y, that give the locations of
    the pair of patches for which commonality is calculated, the distance
    between those cells, and the Sorensen or Jaccard result.

    Notes
    -----
    {2}

    For gridded commonality, cols must also contain x_col and y_col, giving the
    x and y dimensions along which to grid the patch.

    {3}

    """

    (spp_col, count_col, x_col, y_col), patch = \
        _get_cols(['spp_col', 'count_col', 'x_col', 'y_col'], cols, patch)

    # Loop through each split
    result_list = []
    for substring, subpatch in _yield_subpatches(patch, splits):

        # Get spatial table and break out columns
        spatial_table = _yield_spatial_table(subpatch, divs, spp_col,
                                             count_col, x_col, y_col)
        spp_set = spatial_table['spp_set']
        cell_loc = spatial_table['cell_loc']
        n_spp = spatial_table['n_spp']

        # Get all possible pairwise combinations of cells
        pair_list = []
        dist_list = []
        comm_list = []
        for i in range(len(spatial_table)):
            for j in range(i+1, len(spatial_table)):

                iloc = np.round(cell_loc[i], 6)
                jloc = np.round(cell_loc[j], 6)
                pair_list.append('('+str(iloc[0])+' '+str(iloc[1])+') - '+
                                 '('+str(jloc[0])+' '+str(jloc[1])+')')

                dist_list.append(_distance(cell_loc[i], cell_loc[j]))

                ij_intersect = spp_set[i] & spp_set[j]
                if metric.lower() == 'sorensen':
                    comm = 2*len(ij_intersect) / (n_spp[i] + n_spp[j])
                elif metric.lower() == 'jaccard':
                    comm = len(ij_intersect) / len(spp_set[i] | spp_set[j])
                else:
                    raise ValueError, ("Only Sorensen and Jaccard metrics are "
                                      "available for gridded commonality")
                comm_list.append(comm)

        # Append subset result
        subresult = pd.DataFrame({'pair': pair_list, 'x': dist_list,
                                  'y': comm_list})
        result_list.append((substring, subresult))

    # Return all results
    return result_list


def _yield_spatial_table(patch, div, spp_col, count_col, x_col, y_col):
    """
    Calculates an empirical spatial table

    Yields
    -------
    DataFrame
        Spatial table for each division. See Notes.

    Notes
    -----
    The spatial table is the precursor to the SAR, EAR, and grid-based
    commonality metrics. Each row in the table corresponds to a cell created by
    a given division. Columns are cell_loc (within the grid defined by the
    division), spp_set, n_spp, and n_individs.

    """

    div_split_list = div.replace(';','').split(',')
    div_split = (x_col + ':' + div_split_list[0] + ';' +
                 y_col + ':' + div_split_list[1])

    # Get cell_locs
    # Requires _parse_splits and _product functions to go y inside of x
    x_starts, x_ends = _col_starts_ends(patch, x_col, div_split_list[0])
    x_offset = (x_ends[0] - x_starts[0]) / 2
    x_locs = x_starts + x_offset

    y_starts, y_ends = _col_starts_ends(patch, y_col, div_split_list[1])
    y_offset = (y_ends[0] - y_starts[0]) / 2
    y_locs = y_starts + y_offset

    cell_locs = _product(x_locs, y_locs)

    # Get spp set and count for all cells
    n_spp_list = [] # Number of species in cell
    n_individs_list = []
    spp_set_list = []   # Set object giving unique species IDs in cell
    for cellstring, cellpatch in _yield_subpatches(patch,div_split,name='div'):
        spp_set = set(np.unique(cellpatch.table[spp_col]))
        spp_set_list.append(spp_set)
        n_spp_list.append(len(spp_set))
        n_individs_list.append(np.sum(cellpatch.table[count_col]))

    # Create and return dataframe
    df = pd.DataFrame({'cell_loc': cell_locs, 'spp_set': spp_set_list,
                       'n_spp': n_spp_list, 'n_individs': n_individs_list})

    return df




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
            result[row]['dist'] = _decdeg_distance(plot_locs[plota],
                                                   plot_locs[plotb])
        else:
            result[row]['dist'] = _distance(plot_locs[plota],
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

    Notes
    -----
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


def _distance(pt1, pt2):
    """Euclidean distance between two points"""
    return np.sqrt((pt1[0] - pt2[0]) ** 2 + (pt1[1] - pt2[1]) ** 2)


def decdeg_distance(pt1, pt2):
    ''' Calculate Earth surface distance (in km) between decimal latlong points
    using Haversine approximation.

    http://stackoverflow.com/questions/15736995/how-can-i-quickly-estimate-the-distance-between-two-latitude-longitude-points
    '''
    lat1, lon1 = pt1
    lat2, lon2 = pt2

    # Convert decimal degrees to radians
    lon1, lat1, lon2, lat2 = map(np.radians, [lon1, lat1, lon2, lat2])

    # haversine formula
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    a = np.sin(dlat/2)**2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon/2)**2
    c = 2 * np.asin(np.sqrt(a))
    km = 6367 * c

    return km


def _get_cols(special_col_names, cols, patch):
    """
    Retrieve values of special_cols from cols string or patch metadata
    """

    # If cols not given, try to fall back on cols from metadata
    if not cols:
        if 'cols' in patch.meta['Description'].keys():
            cols = patch.meta['Description']['cols']
        else:
            raise NameError, ("cols argument not given, spp_col at a minimum "
                              "must be specified")

    # Parse cols string into dict
    cols = cols.replace(' ', '')
    col_list = cols.split(';')
    col_dict = {x.split(':')[0]: x.split(':')[1] for x in col_list}

    # Get special_col_names from dict
    result = []
    for special_col_name in special_col_names:
        col_name = col_dict.get(special_col_name, None)

        # Create a count col if its requested and doesn't exist
        if special_col_name is 'count_col' and col_name is None:
            col_name = 'count'
            patch.table['count'] = np.ones(len(patch.table))

        # All special cols must be specified (count must exist by now)
        if col_name is None:
            raise ValueError, ("Required column %s not specified" %
                               special_col_name)

        result.append(col_name)

    return tuple(result), patch


@doc_sub(splits_note)
def _yield_subpatches(patch, splits, name='split'):
    """
    Iterator for subtables defined by a splits string

    Parameters
    ----------
    patch : obj
        Patch object containing data to subset
    splits : str
        Specifies how a column of a dataset should be split. See Notes.

    Yields
    ------
    tuple
        First element is subset string, second is subtable dataframe

    Notes
    -----
    {0}

    """

    if splits:
        subset_list = _parse_splits(patch, splits)
        for subset in subset_list:
            log.info('Analyzing subset %s: %s' % (name, subset))
            subpatch = copy.copy(patch)
            subpatch.table = _subset_table(patch.table, subset)
            subpatch.meta = _subset_meta(patch.meta, subset)
            yield subset, subpatch
    else:
        yield '', patch


@doc_sub(splits_note)
def _parse_splits(patch, splits):
    """
    Parse splits string to get list of all associated subset strings.

    Parameters
    ----------
    patch : obj
        Patch object containing data to subset
    splits : str
        Specifies how a column of a dataset should be split. See Notes.

    Returns
    -------
    list
        List of subset strings derived from splits string

    Notes
    -----
    {0}

    """

    split_list = splits.replace(' ','').split(';')
    subset_list = []  # List of all subset strings

    for split in split_list:
        col, val = split.split(':')

        if val == 'split':
            uniques = []
            for level in patch.table[col]:
                if level not in uniques:
                    uniques.append(level)
            level_list = [col + '==' + str(x) + '; ' for x in uniques]
        else:
            starts, ends = _col_starts_ends(patch, col, val)
            level_list = [col + '>=' + str(x) + '; ' + col + '<' + str(y)+'; '
                          for x, y in zip(starts, ends)]

        subset_list.append(level_list)

    # Get product of all string levels as list, conv to string, drop final ;
    return [''.join(x)[:-2] for x in _product(*subset_list)]


def _patch_area(patch, x_col, y_col):

    lengths = []
    for col in [x_col, y_col]:
        col_step = eval(patch.meta[col]['step'])
        col_min = eval(patch.meta[col]['min'])
        col_max = eval(patch.meta[col]['max'])
        lengths.append(col_max - col_min + col_step)

    return lengths[0] * lengths[1]

def _col_starts_ends(patch, col, slices):

    col_step = eval(patch.meta[col]['step'])
    col_min = eval(patch.meta[col]['min'])
    col_max = eval(patch.meta[col]['max'])
    edges = np.linspace(col_min-col_step/2, col_max+col_step/2, eval(slices)+1)
    starts = edges[:-1]
    ends = edges[1:]

    return starts, ends


def _product(*args, **kwds):
    """
    Generates cartesian product of lists given as arguments

    From itertools.product documentation
    """

    pools = map(tuple, args) * kwds.get('repeat', 1)
    result = [[]]
    for pool in pools:
        result = [x+[y] for x in result for y in pool]
    return result


def _distance(pt1, pt2):
    """Euclidean distance between two points"""
    return np.sqrt((pt1[0] - pt2[0]) ** 2 + (pt1[1] - pt2[1]) ** 2)


def _decdeg_distance(pt1, pt2):
    """
    Earth surface distance (in km) between decimal latlong points using
    Haversine approximation.

    http://stackoverflow.com/questions/15736995/
    how-can-i-quickly-estimate-the-distance-between-two-latitude-longitude-
    points
    """

    lat1, lon1 = pt1
    lat2, lon2 = pt2

    # Convert decimal degrees to radians
    lon1, lat1, lon2, lat2 = map(np.radians, [lon1, lat1, lon2, lat2])

    # haversine formula
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    a = np.sin(dlat/2)**2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon/2)**2
    c = 2 * np.arcsin(np.sqrt(a))
    km = 6367 * c

    return km


def empirical_cdf(data):
    """
    Generates an empirical cdf from data

    Parameters
    ----------
    data : iterable
        Empirical data

    Returns
    --------
    DataFrame
        Columns 'data' and 'ecdf'. 'data' contains ordered data and 'ecdf'
        contains the corresponding ecdf values for the data.

    """

    vals = pd.Series(data).value_counts()
    ecdf = pd.DataFrame(data).set_index(keys=0)
    probs = pd.DataFrame(vals.sort_index().cumsum() / np.float(len(data)))
    ecdf = ecdf.join(probs, how="right")
    ecdf = ecdf.reset_index()
    ecdf.columns = ['data', 'ecdf']

    return ecdf
