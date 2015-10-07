from __future__ import division
import os
import re
import copy
from configparser import ConfigParser
import itertools
from copy import deepcopy
import logging

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
    """The parameter ``cols`` is a string describing which column in the data
    table should be used for which "special columns" in analysis. The five
    possible special columns are

    - spp_col - Unique species identifiers
    - count_col - Number of individuals at a location
    - x_col - x coordinate of location
    - y_col - y coordinate of location
    - energy_col - Energetic requirements of individual(s) at a location

    For example, setting ``cols`` to ``spp_col: spp: count_col: number`` will
    use the column named "spp" in the data table to represent the unique
    species identifiers, and the column "number" in the data table to represent
    the count of individuals at a point.

    Different special columns are required for different analyses. count_col is
    used when multiple individuals of a species may be found at a single
    recorded location, as is the case in gridded censuses where all individuals
    in a quadrat are "assigned" to a single point. If count_col is not
    specified, each record in the data table will be presumed to represent a
    single individual (i.e., a count of 1).

    Note that the value of spp_col may be set to a columm in the data table
    giving the genus, family, functional group, etc., which allows for analysis
    of this metric by those groups. """

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
    for a total of 12 separate analyses.  Note that if you pass in the x
    split you MUST also pass in a y split (even if it is just "y:1") or vice
    versa.  Otherwise, the computed areas will be incorrect."""

division_note = \
    """The parameter divisions describes how to successively divide the patch
    along the x_col and y_col dimensions. For
    example, the string '1,2; 2,2; 2,4' will produce an output table with three
    rows, giving the result across two subplots when the patch is split
    along y_col, across four subplots when the patch is split into a 2x2 grid,
    and across eight subplots when the patch is split into 2 parts along x_col
    and 4 parts along y_col."""

start_emp_example = \
    """>>> # Using the ANBO data provided in demo_files_ANBO.zip found at
    >>> # https://github.com/jkitzes/macroeco/releases/

    >>> import macroeco as meco

    >>> # Pass in path to metadata in order to make patch object
    >>> pat = meco.empirical.Patch('~/Desktop/ANBO.txt')"""

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
        Dict-like metadata, loaded from metadata_path and processed by subset
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

    Examples
    --------

    >>> # Using the ANBO data provided in demo_files_ANBO.zip found at
    >>> # https://github.com/jkitzes/macroeco/releases/

    >>> import macroeco as meco

    >>> # Pass in path to metadata in order to make patch object
    >>> pat = meco.empirical.Patch('~/Desktop/ANBO.txt')

    >>> # Subset data upon loading using subset string.
    >>> # Should be a conditional statements separated by a semicolon
    >>> pat = meco.empirical.Patch('~/Desktop/ANBO.txt',
                                    subset="year==2010; row>2")

    """

    def __init__(self, metadata_path, subset=''):

        if not metadata_path:  # Allow for creation of empty patch
            self.meta = None
            self.subset = ''
            self.table = None
        else:
            self.meta = ConfigParser()
            self.meta.read(os.path.expanduser(metadata_path))
            self.subset = subset
            self.table = self._load_table(metadata_path,
                                          self.meta['Description']['datapath'])

        self.incremented = False

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

        metadata_dir = os.path.dirname(os.path.expanduser(metadata_path))
        data_path = os.path.normpath(os.path.join(metadata_dir, data_path))

        extension = data_path.split('.')[-1]

        if extension == 'csv':
            full_table = pd.read_csv(data_path, index_col=False)
            table = _subset_table(full_table, self.subset)
            self.meta, _ = _subset_meta(self.meta, self.subset)
        elif extension in ['db', 'sql']:

            # TODO: deal with incrementing in DB table
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
        col = re.split("[<>=!]", condition)[0]
        comp = condition.replace(col, "")
        try:
            this_valid = eval("full_table['{0}']{1}".format(col, comp))
        except KeyError as e: # catch error and redisplay for twiggy
            raise KeyError("Column '%s' not found" % e.message)
        valid = np.logical_and(valid, this_valid)

    return full_table[valid]

def _subset_meta(full_meta, subset, incremented=False):
    """
    Return metadata reflecting all conditions in subset

    Parameters
    ----------
    full_meta : ConfigParser obj
        Metadata object
    subset : str
        String describing subset of data to use for analysis
    incremented : bool
        If True, the metadata has already been incremented

    Returns
    -------
    Configparser object or dict
        Updated version of full_meta accounting for subset string

    """
    if not subset:
        return full_meta, False

    meta = {}  # Make deepcopy of entire meta (all section dicts in meta dict)
    for key, val in full_meta.iteritems():
        meta[key] = copy.deepcopy(dict(val))

    conditions = subset.replace(' ','').split(';')

    inc = False
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
            if incremented:
                meta[col]['min'] = val
            else:
                meta[col]['min'] = str(eval(val) + eval(col_step))
            inc = True
        elif operator == '<=':
            meta[col]['max'] = val
        elif operator == '<':
            if incremented:
                meta[col]['max'] = val
            else:
                meta[col]['max'] = str(eval(val) - eval(col_step))
            inc = True
        else:
            raise ValueError, "Subset %s not valid" % condition

    return meta, inc


@log_start_end
@doc_sub(metric_params, metric_return, cols_note, splits_note, start_emp_example)
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

    Examples
    --------

    {4}

    >>> # Get the SAD of the full plot
    >>> sad = meco.empirical.sad(pat, 'spp_col:spp; count_col:count', '')

    >>> # Extract the SAD
    >>> sad_df = sad[0][1]
    >>> sad_df
           spp     y
    0    arsp1     2
    1     cabr    31
    2   caspi1    58
    3     chst     1
    4    comp1     5
    5     cran     4
    6     crcr    65
    7    crsp2    79
    8     enfa     1
    9     gnwe    41
    10   grass  1110
    11   lesp1     1
    12    magl     1
    13    mesp     6
    14    mobe     4
    15    phdi   210
    16   plsp1     1
    17    pypo    73
    18    sasp     2
    19    ticr   729
    20   unsh1     1
    21   unsp1    18
    22   unsp3     1
    23   unsp4     1

    >>> # Get SAD for 4 subplots within the full plot and keep absent species
    >>> # using clean = False
    >>> sad_subplots = meco.empirical.sad(pat, 'spp_col:spp; count_col:count', splits = "row:2; column:2", clean=False)
    >>> len(sad_subplots)
    4

    >>> # Look at SAD in one of the 4 cells
    >>> sad_subplots[0]
    ('row>=-0.5; row<1.5; column>=-0.5; column<1.5',
           spp    y
    0    arsp1    0
    1     cabr    7
    2   caspi1    0
    3     chst    1
    4    comp1    1
    5     cran    3
    6     crcr   21
    7    crsp2   16
    8     enfa    0
    9     gnwe    8
    10   grass  236
    11   lesp1    0
    12    magl    0
    13    mesp    4
    14    mobe    0
    15    phdi   33
    16   plsp1    1
    17    pypo    8
    18    sasp    2
    19    ticr  317
    20   unsh1    1
    21   unsp1    0
    22   unsp3    1
    23   unsp4    1)

    See http://www.macroeco.org/tutorial_macroeco.html for additional
    examples and explanation

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
@doc_sub(metric_params, metric_return, cols_note, splits_note, start_emp_example)
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

    Examples
    --------

    {4}

    >>> # Get the spatial abundance distribution for all species for each of
    >>> # the cells in the ANBO plot
    >>> all_spp_ssads = meco.empirical.ssad(pat, cols='spp_col:spp; count_col:count', splits='row:4; column:4')

    >>> # Convert to dict for easy searching
    >>> all_ssads_dict = dict(all_spp_ssads)

    >>> # Look up the spatial abundance distribution for 'grass'
    >>> all_ssads_dict['grass']
         y
    0    42
    1    20
    2    60
    3    60
    4    88
    5    86
    6    20
    7     0
    8   110
    9    12
    10  115
    11  180
    12  160
    13  120
    14   26
    15   11

    >>> # Each value in 'y' gives the abundance of grass in one of the 16 cells

    See http://www.macroeco.org/tutorial_macroeco.html for additional
    examples and explanation

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
@doc_sub(metric_params, metric_return, cols_note, splits_note, division_note,
            start_emp_example)
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
    {1} Result has 5 columns; div, x, and y; that give the ID for the
    division given as an argument, fractional area, and the mean species
    richness at that division.

    Notes
    -----
    {2}

    For the SAR and EAR, cols must also contain x_col and y_col, giving the x
    and y dimensions along which to grid the patch.

    {3}

    {4}

    Examples
    --------

    {5}

    >>> # Get the SAR at the full area (1,1), 1 x 2 division,
    >>> # 2 x 1 division, 2 x 2 division, 2 x 4 division, 4 x 2 division, and
    >>> # 4 x 4 division
    >>> sar = meco.empirical.sar(pat,
                cols='spp_col:spp; count_col:count; x_col:row; y_col:column',
                splits="",
                divs="1,1; 1,2; 2,1; 2,2; 2,4; 4,2; 4,4")

    >>> sar[0][1]
       div  n_individs    n_spp   x        y
    0  1,1   2445.0000  24.0000  16  24.0000
    1  1,2   1222.5000  18.5000   8  18.5000
    2  2,1   1222.5000  17.0000   8  17.0000
    3  2,2    611.2500  13.5000   4  13.5000
    4  2,4    305.6250  10.1250   2  10.1250
    5  4,2    305.6250  10.5000   2  10.5000
    6  4,4    152.8125   7.5625   1   7.5625

    The column div gives the divisions specified in the function call. The
    column n_individs specifies the average number of individuals across the
    cells made from the given division. n_spp gives the average species across
    the cells made from the given division. x gives the absolute area of a
    cell for the given division. y gives the same information as n_spp and is
    included for easy plotting.

    See http://www.macroeco.org/tutorial_macroeco.html for additional
    examples and explanation

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

    # Clean divlist
    return ["".join(s.strip().split(" ")) for s in subdivlist]


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

    # Catch error if you don't use ; after divs in comm_grid in MacroecoDesktop
    try:
        div_split_list = div.replace(';','').split(',')
    except AttributeError:
        div_split_list = str(div).strip("()").split(',')

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


@log_start_end
@doc_sub(metric_params, metric_return, cols_note, splits_note)
def o_ring(patch, cols, splits, spp, bin_edges, density=True, full=False):
    """
    Calculates univariate O-ring for a species

    Parameters
    ----------
    {0}
    bin_edges : iterable
        List of edges of distance classes to bin histogram of distances
    spp : str
        String corresponding to focal species code
    density : bool
        If True, return densities (counts divided by area of torus defined
        by bin edges) instead of counts. Default True.
    full : bool
        If True, return a separate column giving density at distance x for
        every individual, rather than mean density. Default False.

    Returns
    -------
    {1} Result has two columns, x and y, that give the distance to the center
    of a torus and the number or density of individuals found in that torus.

    Notes
    -----
    If density is False, counts are raw counts, non-edge corrected, within
    rings.

    Pairwise distances are directional, giving n(n-1) total distances for a
    species with n individuals, as edge correction is inherently directional.

    Bins include the lower edge and exclude the upper edge, except for the
    final bin which includes both the lower and upper edge. Floating point
    arithmetic may cause points located "exactly" on edges to be allocated
    contrary to this rule, however.

    If there are no records for a species, result table will be a dataframe
    with no records. If there are records but a species has only one
    individual, dataframe will have zero count at all torus areas.

    When using density, the maximum distance used for edge correction, given by
    the mean of the last two bin_edge values, should ideally be set to no
    greater than one half the diagonal distance across the plot. This ensures
    that it is not possible for an entire edge correction buffer to be outside
    of the plot.

    {2}

    For the 0-ring analysis, cols must also contain x_col and y_col, giving the
    x and y dimensions along which to analyze spatial pattern.

    {3}

    """

    try:
        geo.box(0, 0, 1, 1)
    except:
        raise ImportError, "O-ring analysis requires shapely package"

    (spp_col, count_col, x_col, y_col), patch = \
        _get_cols(['spp_col', 'count_col', 'x_col', 'y_col'], cols, patch)

    # Loop through each split
    result_list = []
    for substring, subpatch in _yield_subpatches(patch, splits):

        # Get table for just this species
        spp_table = subpatch.table[subpatch.table[spp_col] == spp]

        # If spp not present, continue
        if (len(spp_table) == 0):
            result_list.append((substring, pd.DataFrame(columns=['x','y'])))
            continue

        # Set up plot geometry
        plot_poly, radii, torus_areas = \
            _get_plot_geometry(subpatch, bin_edges, x_col, y_col)

        # Get lists of all points and counts in spp_table
        x = spp_table[x_col]
        y = spp_table[y_col]
        points = zip(x,y)
        counts = list(spp_table[count_col])

        # Arrays to hold summed areas and distance histograms for all points

        if full:
            hists = []  # Vectors of len(radii) appended for each point
            areas = []
        else:
            hists = np.zeros(len(radii))
            areas = np.zeros(len(radii))

        # Go through each point and associated count
        for i, (point, count) in enumerate(zip(points, counts)):

            # Create list of all other points and counts except this
            other_points = points[0:i] + points[i+1:]
            other_counts = counts[0:i] + counts[i+1:]

            # Get dist from this point to all other points
            if other_points:
                other_dists = dist.cdist(np.array([point]),
                                         np.array(other_points))
            else:
                other_dists = np.array(())

            # Repeat other point distances to acccount for their counts
            other_dists = np.repeat(other_dists, other_counts)

            # Repeat entire other_dist array to account for count here
            other_dists = np.tile(other_dists, count)

            # Add 0's for count at this point to account for count here
            if count > 1:
                other_dists = np.concatenate((other_dists,
                                              np.zeros(count*(count-1))))

            # Calculate histogram of distances to other points
            hist, _ = np.histogram(other_dists, bin_edges)

            # Convert histogram to density if desired
            corr_factor = np.ones(len(radii))  # Frac length in plot
            for j, r in enumerate(radii):
                circ = geo.Point(*point).buffer(r, resolution=64)
                outside_len = circ.boundary.difference(plot_poly).length
                corr_factor[j] = ((circ.boundary.length - outside_len) /
                                   circ.boundary.length)

            # Add hist and corrected area for this point to running totals
            if full:
                hists.append(hist)
                areas.append(torus_areas * corr_factor * count)
            else:
                hists += hist
                areas += torus_areas * corr_factor * count

        # If density, divide summed torus counts by summed areas
        if density:
            hists = np.array(hists) / np.array(areas)

        # Append subset result
        subresult = pd.DataFrame({'x': radii})
        if full:
            for i in range(len(hists)):
                subresult[i] = hists[i]
        else:
            subresult['y'] = hists
        result_list.append((substring, subresult))

    # Return all results
    return result_list


def _get_plot_geometry(subpatch, bin_edges, x_col, y_col):

    # Plot polygon
    xmin = eval(subpatch.meta[x_col]['min'])
    xmax = eval(subpatch.meta[x_col]['max'])
    ymin = eval(subpatch.meta[y_col]['min'])
    ymax = eval(subpatch.meta[y_col]['max'])
    plot_poly = geo.box(xmin, ymin, xmax, ymax)

    # Radii of toruses
    bin_edges = np.array(bin_edges)
    radii = (bin_edges[:-1] + bin_edges[1:]) / 2

    # Areas of all toruses
    torus_areas = []
    for i in range(len(bin_edges) - 1):
        torus_areas.append(np.pi * (bin_edges[i+1]**2 - bin_edges[i]**2))

    return plot_poly, radii, np.array(torus_areas)



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
            logging.info('Analyzing subset %s: %s' % (name, subset))
            subpatch = copy.copy(patch)
            subpatch.table = _subset_table(patch.table, subset)
            subpatch.meta, subpatch.incremented = _subset_meta(patch.meta,
                                                subset, incremented=True)

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

        if patch.incremented:
            lengths.append(col_max - col_min)
        else:
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
    ecdf = ecdf.join(probs)
    ecdf = ecdf.reset_index()
    ecdf.columns = ['data', 'ecdf']

    return ecdf
