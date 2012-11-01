#!/usr/bin/python

'''This module contains 4 separate classes, each built to handle a
canonical data type

This module provides the user with some formatting functions but does provide
the user with all formatting functions that may be required.  This module is
not a substitute for thorough examination of ones data to remove irrelevant
data'''

# NOTE: For some reason, ipython notebook did not like the import * command.
# Need to investigate
import numpy as np
from matplotlib.mlab import csv2rec
from form_func import *
from numpy.lib.recfunctions import drop_fields

__author__ = "Mark Wilber"
__copyright__ = "Copyright 2012, Regents of University of California"
__credits__ = "John Harte"
__license__ = None
__version__ = "0.1"
__maintainer__ = "Mark Wilber"
__email__ = "mqw@berkeley.edu"
__status__ = "Development"

class Columnar_Data:
    '''
    This is the data form that the macroeco software package wants the data
    file in.  All other canonical data sets are converted to columnar data and
    then turned into Columnar_Data objects.

    Examples of columnar data include BCIS, LUQU, and COCO

    Multiple data files must have same format if they are to be merged

    '''

    def __init__(self, datalist, delimiter=',', missingd=None,\
                delete_missing=True, archival=True):
        '''
        This __init__ method takes in data and stores it in rec_arrays.
        If specified,  it will located missing data points and remove them
        from the data set.

        Parameters
        ----------
        datalist : string, list of strings, or list of ndarrays.
            Data filenames or list of data arrays

        delimiter : string
            The file delimiter. Default is ','

        missingd : dict
            Dictionary mapping munged column names to field values which 
            signify that the field does not contain actual data and should be
            masked, e.g. '0000-00-00' or 'unused'
        
        delete_missing : bool
            If True, deletes all of the missing values. If False, only deletes
            the NaNs from the data.

        archival : bool
            If True, a copy of self.columnar_data is made and stored in
            self.columnar_archival. If dataset is very large, set to False.

        Note
        ----
        If column type is integer, missing values are set to -1.  If column
        type is float, missing values are set to NaN.  If column type is
        string, missing values are set to ''.  If column type is object,
        missing values are set to None.

        '''
        if type(datalist) == str:
            datalist = [datalist]

        if np.all(np.array([type(x) == str for x in datalist])):
            self.columnar_data = []
            self.data_names = []
            for file_name in datalist:
                self.columnar_data.append(csv2rec(file_name, delimiter=delimiter,\
                              missingd=missingd))
                self.data_names.append(file_name)
            if missingd != None:
                if delete_missing:
                    trun_data = []
                    for data in self.columnar_data:
                        for key in missingd.iterkeys():
                            try:
                                # Missing float
                                notNaN = (False == np.isnan(data[key]))
                            except:
                                notNaN = np.ones(len(data[key]), dtype=bool)
                            notBlank = np.array([it != '' for it in data[key]])
                            notMinusOne = (data[key] != -1)# Missing int
                            # Missing other
                            notNone = np.array([i != None for i in data[key]])
                            ind = np.bitwise_and(notNaN, notBlank)
                            ind = np.bitwise_and(ind, notMinusOne)
                            ind = np.bitwise_and(ind, notNone)
                            data = data[ind]
                        trun_data.append(data)
                    self.columnar_data = trun_data
                else:
                    trun_data = []
                    for data in self.columnar_data:
                        for key in missingd.iterkeys():
                            try:
                                notNaN = (False == np.isnan(data[key]))
                            except:
                                notNaN = np.ones(len(data[key]), dtype=bool)
                            data = data[notNaN]
                        trun_data.append(data)
                    self.columnar_data = trun_data
        elif np.all(np.array([type(x) == np.ndarray for x in datalist])):
            self.columnar_data = datalist

        if archival:
            self.columnar_archival = [np.copy(data) for data in 
                                                            self.columnar_data]
        else:
            self.columnar_archival = []

    def reset_columnar_data(self):
        '''
        Resets self.columnar_data to self.columnar_archival
        
        Need to be careful about excessive memory usage!
        '''
        if len(self.columnar_archival) == 0:
            raise ValueError("The self.columnar_archival attribute of this %s"
                             % (self.__class__.__name__) + " object has not" 
                             + " been initialized")
        else:
            self.columnar_data = [np.copy(data) for data in 
                                                        self.columnar_archival]

    def split_up_data_by_field(self, split_columns=None):
        '''
        This function will take in the split-columns list and and split the
        data into separate arrays based on the list.  For example, if one were
        to pass in dbh1, dbh2,  dbh3 three copies of the data would be
        made, each being identical except that each would only contain one of
        the instances of dbh. One could also pass [(dbh1, recr1), (dbh2, recr2),
        (dbh3, recr3)].  All other fields in split_columns will be excluded
        other than the fields within the tuple under consideration.

        Parameters
        ----------
        split_columns : list
            a list of tuples specifying the columns by which to split the array
        
        Notes
        -----
        Saves the split array as self.columnar_data.
        
        '''
        #Note: If they enter the wrong column name nothing will be removed
        #Should I error check for this?
        if split_columns != None:
            split_data = []
            given_col_names = []
            for tup in split_columns:
                for name in tup:
                    given_col_names.append(name)
            given_col_names = np.array(given_col_names)


            for data in self.columnar_data:
                for tup in split_columns:
                    ind = np.ones(len(given_col_names), dtype=bool)
                    for name in tup:
                        ind = np.bitwise_and((name != given_col_names), ind)
                    remove_names = given_col_names[ind]
                    split_data.append(drop_fields(data, list(remove_names)))
            self.columnar_data = split_data
    
    def change_column_names(self, change=None, changed_to=None):
        '''
        This function takes a list of column names to be changed and a name
        that they should be changed to

        Parameters
        ----------
        change : list
            List of column names.  Columns names are strings
        changed_to : string
            Name to be changed to

        Notes
        -----
        This function useful if you would like to merge self.columnar_data but
        the dtype.names are different.

        '''
        if change != None and changed_to != None: 
            for data in self.columnar_data:
                column_names = np.array(data.dtype.names)
                for name in change:
                    find = np.where((name == column_names))[0]
                    if len(find) != 0:
                        column_names[find[0]] = changed_to
                        data.dtype.names = tuple(column_names)
        
    def add_fields_to_data_list(self, fields=None, values=None):
        '''
        This functions adds given fields and values to the data list. The
        length of values should be the same length as fields and the length of
        each tuple within each element of values should be the same length as
        the self.columnar_data

        Parameters
        ----------
        fields : list
            A list of strings specifying field names
        values : list of tuples
            A list of tuples with the length of each tuple equalling the length
            of self.columnar_data

        '''
        #NOTE: Should probably make a single dictionary for field/values
        if fields != None and values != None:
            self.columnar_data = add_data_fields(self.columnar_data, fields,\
                                                                        values)

    def remove_columns(self, col_names=None):
        '''
        This function will remove the all the columns within with names in
        col_names from all the datasets in self.columnar_data.

        Parameters
        ----------
        col_names : string or list
            The name or names of columns to be removed

        '''
        
        if col_names != None:
            if type(col_names) == str:
                col_names = [col_names]
            else:
                col_names = list(col_names)
            removed_data = []
            for data in self.columnar_data:
                removed_data.append(drop_fields(data, col_names))
            self.columnar_data = removed_data

    def fractionate_data(self, wid_len=None, step=None, col_names=None):
        '''
        This function converts grid numbers to length measurements in
        self.columnar_data

        Parameters
        ----------
        wid_len : tuple
            A tuple containing the the absolute length of the columns being
            converted
        step : tuple
            The precision (step or stride length) of each grid.  The first
            element in the step tuple corresponds with the first element in the
            wid_len tuple and so on.
        col_names : array-like object
            An array-like object of strings giving the names of the columns
            that will be fractionated

        '''
        if wid_len != None and step != None and col_names != None:
            self.columnar_data = fractionate(self.columnar_data, wid_len, step,
                                                                     col_names)


    def merge_data(self):
        '''
        This function concatenates the data files in data_list.  The dtypes of
        the data in data_list must be identical or this function will fail.
        '''

        self.merged_data = merge_formatted(self.columnar_data)

    def output_merged_data(self, filename):
        '''
        This function merges self.columnar_data and outputs the merged data.

        Parameters
        ----------
        filename : string
            The filename to be output

        '''
        #Merge data in case it has not been done
        self.merge_data()
        output_form(self.merged_data, filename)

    def output_columnar_data(self, filenames):
        '''
        This function outputs the self.columnar_data

        Parameters
        ----------
        filenames : list
            A list of filenames

        '''
        assert len(filenames) == len(self.columnar_data), "Number of filenames\
                                 must be the same as the number of datasets"
        for i, name in enumerate(filenames):
            output_form(self.columnar_data[i], name)

class Grid_Data:
    '''This class handles data should look like the EarthFlow data after a 
    census.  It is a grid with species abundance data in each cell. 
    ex.
    ARTDRA - 6
    GERTYR - 8

    NOTE: Need to consider things that might break this class
    '''

    def __init__(self, filenames, num_cols, archival=True):
        '''
        Pass in the file name(s) of the grid data that you want converted and
        the number of columns in each grid.

        Parameters
        ----------

        filenames : str or list of strings
            A filename or list of filenames

        num_cols : int or list of ints
            If an int or list of length one, this number specifies the number
            of columns for all grids with in filenames.  If a list of
            len(filename) this list specifies the number of columns for each
            individual grid.

        archival : bool
            If True, a copy of self.grid_data is made and stored in
            self.grid_archival. If dataset is very large, set to False.

        '''
        #NOTE: Handle missing data!!!!

        if type(filenames) == str:
            filenames = [filenames]

        if type(num_cols) == int:
            num_cols = [num_cols]

        assert np.all(np.array([name.split('.')[-1] for name in filenames]) ==\
                      'csv'), "Files must be csv"
        assert len(num_cols) == len(filenames) or len(num_cols) == 1, 'Length\
                       of num_cols must be 1 or equal len(filenames)'

        self.grid_data = []
        self.cols = []
        self.rows =[]

        for i, name in enumerate(filenames):
            if len(num_cols) == 1:
                self.cols.append(num_cols[0])
            else:
                self.cols.append(num_cols[i])
            self.grid_data.append(csv2rec(name, names=list(np.arange(0,\
                                            self.cols[i]).astype('S10'))))
            self.rows.append(len(self.grid_data[i]))

        #Remove all '\n' from the end of each cell in grid
        #Not technically necessary but just being clean
        self.grid_data = remove_char(self.grid_data)
        self.grid_data = remove_white_spaces(self.grid_data)

        if archival == True:
            self.grid_archival = [np.copy(data) for data in self.grid_data]
        else:
            self.grid_archival = []

    def reset_grid_data(self):
        '''
        Resets self.grid_data to self.archival_data
        
        Need to be careful about excessive memory usage!
        '''

        if len(self.grid_archival) == 0:
            raise ValueError("The self.grid_archival attribute of this %s"
                             % (self.__class__.__name__) + " object has not" 
                             + " been initialized")
        else:
            self.grid_data = [np.copy(data) for data in self.grid_archival]

    def truncate_grid_cells(self, symbol=None):
        '''
        This function will look at each cell in grid list and truncated the
        string within the cell at AND after the first instance of a given
        symbol.

        Parameters
        ----------
        symbol : char (string of length one)
            The symbol at which to being truncation

        Notes
        -----
        symbol is a keyword argument because format_grid_data script gives the
        option to run every method.

        '''
        if symbol != None: 
            for i in xrange(len(self.grid_data)):
                for nm in self.grid_data[i].dtype.names:
                    for j in xrange(len(self.grid_data[i][nm])):
                        ind = self.grid_data[i][nm][j].find(symbol)
                        if ind != -1:
                            self.grid_data[i][nm][j] = \
                                                 self.grid_data[i][nm][j][:ind]

            self.grid_data = remove_char(self.grid_data)

    def remove_and_replace(self, remove=None, replace=None):
        '''
        Removes a string from a grid cell and replaces it with another one

        Paramters
        ---------
        remove : string
            String to be removed
        replace : string
            String to replace removed string

        '''
        
        if remove != None and replace != None:
            for i in xrange(len(self.grid_data)):
                for nm in self.grid_data[i].dtype.names:
                    for j in xrange(len(self.grid_data[i][nm])):
                        ind = self.grid_data[i][nm][j].find(remove)
                        if ind != -1:
                            self.grid_data[i][nm][j] =\
                              self.grid_data[i][nm][j].replace(remove, replace)

    def find_unique_spp_in_grid(self, spacer='-', spp_sep='\n'):
        '''
        This function finds all of the unique species in the grid.
        It assumes that your grid data is in the proper format.

        Parameters
        ----------
        spacer : str
            The character separating the species code from the species count.
            Default value is '-' (n-dash)

        spp_sep : str
            The character that separates a speces/count combination from
            another species/count combination.  Default value is '\n'

        '''
        self.unq_spp_lists = []
        for num, data in enumerate(self.grid_data):
            spp_names = []
            for col in data.dtype.names:
                for row in xrange(self.rows[num]):
                    if data[col][row].find(spacer) != -1:
                        nam_lst = data[col][row].split(spacer)
                        if len(nam_lst) == 2:
                            spp_names.append(nam_lst[0].strip())
                        else:
                            spp_names.append(nam_lst[0].strip())
                            for i in xrange(1, len(nam_lst) - 1):
                                spp_names.append(nam_lst[i].split(spp_sep)[1].\
                                                                    strip())
            self.unq_spp_lists.append(np.unique(np.array(spp_names)))

    def grid_to_dense(self, spacer='-', spp_sep='\n', archival=True):
        '''
        This function converts a the list of gridded data sets into dense 
        data sets and stores them in dense_data.  In addition, it
        makes a Dense_Data object out of the newly converted data.

        Parameters
        ----------
        spacer : str
            The character separating the species code from the species count.
            Default value is '-' (n-slash)

        spp_sep : str
            The character that separates a speces/count combination from
            another species/count combination.  Default value is '\n'


        '''

        self.find_unique_spp_in_grid(spacer=spacer, spp_sep=spp_sep)
        dense_data = []
        for i, data in enumerate(self.grid_data):
            dtype_list = [('cell', np.int), ('row', np.int), ('column', np.int)]
            for name in self.unq_spp_lists[i]:
                tuple_type = (name, np.int)
                dtype_list.append(tuple_type)
            matrix = np.empty(self.rows[i] * self.cols[i], dtype=dtype_list)
            #Iterate through the plot
            count = 0
            for col in data.dtype.names:
                for row in xrange(self.rows[i]):
                    matrix['cell'][count] = count
                    matrix['row'][count] = row
                    matrix['column'][count] = int(col)
                    for spp_name in self.unq_spp_lists[i]:

                        # Check if cell has species
                        start = data[col][row].find(spp_name)
                        if start != -1:
                            # could be nested in another word
                            if data[col][row][start + len(spp_name)] == spacer:
                               # The nesting could be at the end of the word
                               if start == 0 or data[col][row][start - 1] ==\
                                                                       spp_sep:

                                    raw = data[col][row][start:].split(spacer)[1]
                                    if raw.find(spp_sep) != -1:
                                        tot_spp = raw.split(spp_sep)[0].strip()
                                        matrix[spp_name][count] = int(tot_spp)
                                    else:
                                        tot_spp = raw.split()[0].strip()
                                        matrix[spp_name][count] = int(tot_spp)
                               else:
                                    matrix[spp_name][count] = 0
                            else:
                                matrix[spp_name][count] = 0
                        else:
                            matrix[spp_name][count] = 0
                    count += 1
            dense_data.append(matrix)
        self.Dense_Object = Dense_Data(dense_data, archival=archival)

                    
    def output_grid_data(self, filenames):
        '''
        This function prints the data within self.grid_data with the given
        filenames.

        Parameters
        -----------
        filenames : list
            A list of filnames to which the data will be saved

        '''

        assert len(filenames) == len(self.grid_data), "Number of filenames\
                                 must be the same as the number of datasets"
        for i, data in enumerate(self.grid_data):
            output_form(data, filenames[i]) 

    
class Dense_Data:
    '''This class handles data that are in the dense format.  An example of the
    dense format is a csv file that has columns named 'row' and 'column' and
    the remainder of columns named after each species in the plot.  The values
    within each species column are the counts within the cell specified by the
    columns names 'row' and 'column'.

    Note: Need to consider how I might break this class
    '''

    def __init__(self, datalist, delim=',', replace=None, archival=True):
        '''

        Parameters
        -----------
        datalist : string, list of strings or list of arrays
            List of filenames to be loaded or list of arrays to be set to
            self.dense_data
        delim : string
            The file delimiter
        missingd : dict
            A dictionary mapping column names to values which signfy the that
            the value is missing in the field. If dtype is string value is set
            to '', if dtype is int value is set to -1, if dtype is float value
            is set to NaN, if dtype is object value is set to None.
        replace : tuple
            A tuple of length 2.  The first element is a string which
            represents the missing values that you would like to replace.  The
            second element is the value with which you would like to replace
            the missing values.
        archival : bool
            If True, a copy of self.dense_data is made and stored in
            self.dense_archival. If dataset is very large, set to False.

        '''
        #TODO: What kind of files could break this
        if type(datalist) == str:
            datalist = [datalist]

        if np.all(np.array([type(x) == str for x in datalist])):
            self.dense_data = []
            if replace != None:
                assert len(replace) == 2, "Replace must contain 2 elements"

                for name in datalist:
                    data = csv2rec(name, delimiter=delim, missing=replace[0])
                    for nm in data.dtype.names:
                        try:
                            # Missing float
                            isNaN = (np.isnan(data[nm]))
                        except:
                            isNaN = np.zeros(len(data[nm]), dtype=bool)
                        isBlank = np.array([it == '' for it in data[nm]])
                        isMinusOne = (data[nm] == -1)# Missing int
                        # Missing other
                        isNone = np.array([i == None for i in data[nm]])
                        ind = np.bitwise_or(isNaN, isBlank)
                        ind = np.bitwise_or(ind, isMinusOne)
                        ind = np.bitwise_or(ind, isNone)
                        data[nm][ind] = replace[1]
                    self.dense_data.append(data)
            else:
                for name in datalist:
                    data = csv2rec(name, delimiter=delim)
                    self.dense_data.append(data)

        elif np.all(np.array([type(x) == np.ndarray for x in datalist])):
            self.dense_data = datalist

        if archival:
            self.dense_archival = [np.copy(data) for data in 
                                                            self.dense_data]
        else:
            self.dense_archival = []

    def reset_grid_data(self):
        '''
        Resets self.grid_data to self.archival_data
        
        Need to be careful about excessive memory usage!
        '''

        if len(self.dense_archival) == 0:
            raise ValueError("The self.dense_archival attribute of this %s"
                             % (self.__class__.__name__) + " object has not" 
                             + " been initialized")
        else:
            self.dense_data = [np.copy(data) for data in self.dense_archival]


    def dense_to_columnar(self, spp_col_num, num_spp, archival=True):
        '''
        This function uses a function in form_func to convert dense data into
        columnar data. Stores the columnar data as a Columnar Object.

        Parameters
        ----------
        spp_col_num : int
            The column number in the dense array where the spp_names begin

        num_spp : tuple
            Number of species in each dataset in self.dense_data

        '''

        columnar_data = format_dense(self.dense_data, spp_col_num,\
                                                                      num_spp)
        self.Columnar_Object = Columnar_Data(columnar_data, archival=archival)

    def output_dense_data(self, filenames):
        '''
        This function prints the data within self.dense_data with the given
        filenames.  If self.dense_data has not been filled, error is thrown.

        Parameters
        ----------
        filenames : list
            A list of filenames to which the data will be saved
    
        '''

        assert len(filenames) == len(self.dense_data), "Number of filenames\
                                 must be the same as the number of datasets"
        for i, data in enumerate(self.dense_data):
            output_form(data, filenames[i])

class Transect_Data:
    '''
    This class handles data that are similar to the Breeding Bird survey data.
    One column has the species ID, one column has stop and all the other
    columns have transects.  This class can handle data with "n" nestings, not
    just two.  For example, the data could have location, transect and stop.

    The "stop" data should all be in consecutive columns

    '''
    
    def __init__(self, filenames, delimiter=',', archival=True):
        '''

        Parameters
        ----------
        filenames : list
            A list of filenames
        delimiter : string
            The file delimiter
        archival : bool
            If True, a copy of self.transect_data is made and stored in
            self.transect_archival. If dataset is very large, set to False.


        '''
        self.columnar_data = []
        self.transect_data = []
        if type(filenames) == str:
            filenames = [filenames]

        for name in filenames:
            data = csv2rec(name, delimiter=delimiter)
            self.transect_data.append(data)

        if archival:
            self.transect_archival = [np.copy(data) for data in 
                                                            self.transect_data]
        else:
            self.transect_archival = []

    def reset_transect_data(self):
        '''
        Resets self.transect_data to self.transect_archival
        
        Need to be careful about excessive memory usage!
        '''
        if len(self.transect_archival) == 0:
            raise ValueError("The self.transect_archival attribute of this %s"
                             % (self.__class__.__name__) + " object has not" 
                             + " been initialized")
        else:
            self.transect_data = [np.copy(data) for data in 
                                                        self.transect_archival]

    def transect_to_columnar(self, stop_col_num, stop_name, tot_stops,\
                                                    count_name='count'):
        '''
        This function takes transect data and convertes it into columnar data.
        In addition it saves the columnar data as a Columnar_Data object. 
        

        Parameters
        ----------
        stop_col_num : int
            The column number where the stop counts begin (1 is the first
            column)
        
        stop_name : str
            The name of the new stop column in the formatted data

        tot_stops : int
            The number of columns with stops

        count_name : str
            The name of the count column. Default is "count"


        Notes
        -----
        This function assumes that all data in self.transect_data are formatted
        the same way.  For example, the column that contains species names or
        codes has the same name throughout all data sets.

        '''
        for data in self.transect_data:
            nstops = tot_stops
            dtypes = data.dtype.descr[ : stop_col_num - 1]
            if (len(dtypes) + nstops) != len(data.dtype.names):
                #Accounting for data fields after stops
                end_dtypes = data.dtype.descr[(len(dtypes) + nstops) : ]
                for x in end_dtypes:
                    dtypes.append(x)
            dtypes.append((stop_name, 'S20'))
            dtypes.append((count_name, np.int))
            column_data = np.empty(len(data) * nstops, dtype=dtypes)
            for i in xrange(len(data)):
                for name in column_data.dtype.names:
                    if name is stop_name:
                        column_data[name][i * nstops:(i + 1) * nstops] = \
                                                           np.arange(0, nstops)
                    elif name is count_name:
                        column_data[name][i * nstops:(i + 1) * nstops] = \
                                   np.array(list(data[i]))[stop_col_num - 1 : \
                                                              -len(end_dtypes)]
                    else:
                        column_data[name][i * nstops:(i + 1) * nstops] = \
                                                                  data[name][i]
            self.columnar_data.append(column_data)
        self.Columnar_Object = Columnar_Data(self.columnar_data)

    def output_transect_data(self, filenames):
        '''
        This function prints the data within self.columnar_data with the given
        filenames.  If self.columnar_data has not been filled, an error is 
        thrown.

        Parameters
        ----------
        filenames : list
            A list of filenames to which the data will be saved. Must be the
            same length as self.columnar_data

        '''

        assert len(filenames) == len(self.transect_data), "Number of filenames\
                                 must be the same as the number of datasets"
        for i, data in self.transect_data:
            output_form(data, filenames[i]) 


def remove_char(grid_list, char='\n'):
    '''
    Removes the given char from the end of each cell in grid list
    '''

    for grid in grid_list:
        for name in grid.dtype.names:
            for i in xrange(len(grid[name])): 
                while grid[name][i][::-1].find('\n') == 0:
                    grid[name][i] = grid[name][i][:-1]
    
    return grid_list

def remove_white_spaces(grid_list):
    '''
    Removes all of the white spaces from strings.
    '''
    for grid in grid_list:
        for name in grid.dtype.names:
            for i in xrange(len(grid[name])): 
                grid[name][i] = ''.join(grid[name][i].strip(' '))

    return grid_list



                     

            



            



















       









            





