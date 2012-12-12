#!/usr/bin/python

'''This module contains 4 separate classes, each built to handle a
canonical data type

This module provides the user with some formatting functions but does provide
the user with all formatting functions that may be required.  This module is
not a substitute for thorough examination of ones data to remove irrelevant
data'''

import numpy as np
from matplotlib.mlab import csv2rec
import form_func as ff
from numpy.lib.recfunctions import drop_fields
import csv

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
                delete_missing=False, archival=True):
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

    def subset_data(self, subset={}):
        '''
        Subset any given column of the data
        
        Parameters
        ----------
        subset : dict
            Dictionary of permanent subset to data, {'column_name':
            'condition'}, which will limit all analysis to records in which
            column_name meets the condition, ie, {'year': ('==', 2005), 'x':
            [('>', 20), ('<', 40)]} restricts analysis to year 2005 and x
            values between 20 and 40. These conditions can also be passed to
            the individual methods, but subsetting the data table up front may
            save analysis time.  Subsetting on a string would look something
            like {'name' : [('==', 'John'), ('==', 'Harry')]}
        '''

        if subset != {}:
            sub_data = []
            for data in self.columnar_data:
                valid = np.ones(len(data), dtype=bool)

                for key, value in subset.iteritems():
                    if type(value) is not type(['a']):  # Make all iterables
                        value = [value]

                    # Merge tuples into a string
                    merged_values = []
                    for val in value:
                        try: # check if val[1] is a string
                            eval(str(val[1]))
                            merged_values.append(val[0] + str(val[1]))
                        except:
                            merged_values.append(val[0]  + "'" +  val[1] + "'")

                    for this_value in merged_values:
                        try:
                            this_valid = eval("data[key]" + this_value)
                            valid = np.logical_and(valid, this_valid)
                        except ValueError: #If key can't be found do nothing
                            pass
                                        
                sub_data.append(data[valid])

            self.columnar_data = sub_data

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
            # Check if split_columns is a list of strings. If so, change it
            # into a list of tuples
            split_columns = [(s,) if type(s) == str else tuple(s) for s in 
                                                                 split_columns]

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
        change : list of tuples or strings
            Each tuple or string contains column names. All the column names in
            the first tuple will be changed to the first element in the
            changed_to list and so on.
        changed_to : list
            A list of strings that contain the names that the columns in change
            will be changed to. 

        Notes
        -----
        This function is useful if you would like to merge self.columnar_data
        but the dtype.names are different.

        '''
        if change != None and changed_to != None: 
            if len(change) != len(changed_to):
                raise ValueError('Length of params change and changed_to must'
                                + ' be equal')
            # Convert to tuples if just received strings
            change = [(x,) if type(x) == str else tuple(x) for x in change]

            for data in self.columnar_data:
                column_names = np.array(data.dtype.names)
                for i, name_tup in enumerate(change):
                    for name in name_tup:
                        find = np.where((name == column_names))[0]
                        if len(find) != 0:
                            max_len = np.max([len(x) for x in column_names])
                            if max_len < len(changed_to[i]):
                                column_names = column_names.astype('S' +
                                                       str(len(changed_to[i])))
                            column_names[find[0]] = changed_to[i]
                            data.dtype.names = tuple(column_names)
        
    def add_fields_to_data_list(self, fields_values=None, descr='S20'):
        '''
        This functions adds given fields and values to the data list. The
        length of values should be the same length as fields and the length of
        each tuple within each element of values should be the same length as
        the self.columnar_data

        Parameters
        ----------
        fields_values : dict
            dictionary with keyword being the the field name to be added and
            the value being a tuple with length self.columnar_data specifying
            the values to be added to each field in each data set.
         descr : list of data types or single data type
        descr : a single data type or a dictionary
            A single value will be broadcast to appropriate length.  The
            dictionary must have the same keywords as fields_values and must be
            the same length.  Each keyword should lookup a dtype.
        '''
        if fields_values != None:
            self.columnar_data = ff.add_data_fields(self.columnar_data,
                                                fields_values, descr=descr)

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

    def fractionate_data(self, wid_len=None, step=None, col_names=None,
                                wid_len_old=None, min_old=None, step_old=None):
        '''
        This function converts grid numbers to length measurements in
        self.columnar_data

        Parameters
        ----------
        wid_len : tuple
            A tuple containing the the absolute length of the columns being
            converted
        step : tuple
            The desierd precision (step or stride length) of each grid.  The
            first element in the step tuple corresponds with the first element
            in the wid_len tuple and so on.
        col_names : array-like object
            An array-like object of strings giving the names of the columns
            that will be fractionated
        wid_len_old : tuple or None
            If None, it assumes that a np.unique on datayears[col_name[i]]
            gives a array that is the same length as np.arange(0,
            wid_len_new[i], step=step_new[i]).  If it doesn't, an error will be
            thrown.  If not None, expects the old maximum length for the given
            columns. 
        min_old : tuple or None
            Same as wid_len_old but the old minimum value for each given column
        step_old : tuple or None
            Same as wid_len_old but the old step (or stride length/spacing) for
            each given column.

        '''
        if wid_len != None and step != None and col_names != None:
            self.columnar_data = ff.fractionate(self.columnar_data, wid_len, step,
                                 col_names, wid_len_old=wid_len_old,
                                 min_old=min_old, step_old=step_old)


    def merge_data(self):
        '''
        This function concatenates the data files in data_list.  The dtypes of
        the data in data_list must be identical or this function will fail.
        '''

        self.merged_data = ff.merge_formatted(self.columnar_data)

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
        ff.output_form(self.merged_data, filename)

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
            ff.output_form(self.columnar_data[i], name)

class Grid_Data:
    '''This class handles data should look like the EarthFlow data after a 
    census.  It is a grid with species abundance data in each cell. 
    ex.
    ARTDRA - 6
    GERTYR - 8

    '''

    def __init__(self, filenames, archival=True, spp_sep='\n'):
        '''
        Pass in the file name(s) of the grid data that you want converted and
        the number of columns in each grid.

        Parameters
        ----------

        filenames : str or list of strings
            A filename or list of filenames

        archival : bool
            If True, a copy of self.grid_data is made and stored in
            self.grid_archival. If dataset is very large, set to False.

        '''
        #NOTE: Handle missing data!!!!

        if type(filenames) == str:
            filenames = [filenames]

        assert np.all(np.array([name.split('.')[-1] for name in filenames]) ==\
                      'csv'), "Files must be csv"

        self.grid_data = []
        self.cols = []
        self.rows =[]

        for i, name in enumerate(filenames):
            # Sometimes csv.reader reads an extra column so you have to read to
            # whole file. Seems stupid to read in the file twice but oh well...
            with open(name, 'rb') as csvreader:
                reader = csv.reader(csvreader)
                rows = [row for row in reader]
            min_len = np.min([len(row) for row in rows])
            self.cols.append(min_len)
                
            self.grid_data.append(csv2rec(name, names=list(np.arange(0,\
                                            self.cols[i]).astype('S10'))))
            self.rows.append(len(self.grid_data[i]))

        #Remove all '\n' from the end of each cell in grid
        #Not technically necessary but just being clean
        self.grid_data = remove_char(self.grid_data, char=spp_sep)
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
        symbol : string or list of strings
            The symbol at which to being truncation

        Notes
        -----
        symbol is a keyword argument because format_grid_data script gives the
        option to run every method.

        '''
        if symbol != None: 

            if type(symbol) == str:
                symbol = [symbol]
            else:
                symbol = list(symbol)

            for i in xrange(len(self.grid_data)):
                for nm in self.grid_data[i].dtype.names:
                    for j in xrange(len(self.grid_data[i][nm])):
                        for sym in symbol:
                            ind = self.grid_data[i][nm][j].find(sym)
                            if ind != -1:
                                self.grid_data[i][nm][j] = \
                                                 self.grid_data[i][nm][j][:ind]

            self.grid_data = remove_char(self.grid_data)
    
    # List of remove replace tuples?
    def remove_and_replace(self, remove=None, replace=''):
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
                tuple_type = (name, np.float)
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

                        # Check if cell has species. May be nested occurence!
                        matrix[spp_name][count] = 0 # Set base to 0
                        start = data[col][row].find(spp_name)
                        if start == -1: # Nothing is there
                            pass # Count already set to zero

                        else: # Something is there, but is it nested?
                            found = start
                            while found != -1:
                                # If this is true, it is nested
                                if (data[col][row][start + len(spp_name)] !=\
                                    spacer) or not(start == 0 or \
                                    data[col][row][start - 1] == spp_sep):

                                    pass

                                else: # Actually a species, so add some
                                      # abundance

                                    raw = data[col][row][start:].split(spacer)[1]
                                    if raw.find(spp_sep) != -1:
                                        tot_spp = raw.split(spp_sep)[0].strip()
                                    else:
                                        tot_spp = raw.split()[0].strip()
                                    matrix[spp_name][count] += float(tot_spp)
                                found = data[col][row][start + 1
                                                              :].find(spp_name)
                                start += found + 1
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
            ff.output_form(data, filenames[i]) 

    
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
        replace : tuple
            A tuple of length 2.  The first element is a string that
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
                    self.dense_data.append(replace_vals(name, replace,
                                                                  delim=delim))
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


    def dense_to_columnar(self, spp_col_num, num_spp, count_col='count',\
                                                                archival=True):
        '''
        This function uses a function in form_func to convert dense data into
        columnar data. Stores the columnar data as a Columnar Object.

        Parameters
        ----------
        spp_col_num : int
            The column number in the dense array where the spp_names begin

        num_spp : tuple or int
            Number of species in each dataset in self.dense_data. If it is an
            int, it will be broadcasted to the length of self.dense_data

        count_col : str
            This string specifies the name of the count column.  The default is
            'count'.

        '''
        columnar_data = ff.format_dense(self.dense_data, spp_col_num,\
                                            num_spp, count_col=count_col)
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
            ff.output_form(data, filenames[i])

class Transect_Data:
    '''
    This class handles data that are similar to the Breeding Bird survey data.
    One column has the species ID, one column has stop and all the other
    columns have transects.  This class can handle data with "n" nestings, not
    just two.  For example, the data could have location, transect and stop.

    The "stop" data should all be in consecutive columns

    '''
    
    def __init__(self, filenames, delim=',', replace=None, archival=True):
        '''

        Parameters
        ----------
        filenames : list
            A list of filenames
        delim : string
            The file delimiter
        replace : tuple
            A tuple of length 2.  The first element is a string which
            represents the missing values that you would like to replace.  The
            second element is the value with which you would like to replace
            the missing values.
        archival : bool
            If True, a copy of self.transect_data is made and stored in
            self.transect_archival. If dataset is very large, set to False.


        '''
        self.transect_data = []
        if type(filenames) == str:
            filenames = [filenames]

        if replace != None:
                
            assert len(replace) == 2, "Replace must contain 2 elements"
            replace = (str(replace[0]), replace[1])

            for name in filenames:
                self.transect_data.append(replace_vals(name, replace,
                                                                  delim=delim))
        else:
            for name in filenames:
                data = csv2rec(name, delimiter=delim)
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

    def transect_to_columnar(self, stop_col_num, tot_stops, stop_name='stop',
                                     count_name='count', archival=True):
        '''
        This function takes transect data and convertes it into columnar data.
        In addition it saves the columnar data as a Columnar_Data object. 
        

        Parameters
        ----------
        stop_col_num : iterable or int
            The column number where the stop counts begin (0 is the first
            column). Can be len(transect_data) or length == 1. Broadcast if
            length equals 1.

        tot_stops : iterable or int
            The number of columns with stops. Can be len(transect_data) or 
            length == 1. Broadcast if length equals 1.
        
        stop_name : str
            The name of the new stop column in the formatted data

        count_name : str
            The name of the count column. Default is "count"


        Notes
        -----
        This function assumes that all data in self.transect_data are formatted
        the same way.  For example, the column that contains species names or
        codes has the same name throughout all data sets.

        '''
        # Broadcast stop_col_num
        stop_col_num = ff.broadcast(len(self.transect_data), stop_col_num)
        tot_stops = ff.broadcast(len(self.transect_data), tot_stops)

        columnar_data = []
        for j, data in enumerate(self.transect_data):
            nstops = tot_stops[j]
            dtypes = data.dtype.descr[ : stop_col_num[j] ]
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
                                   np.array(list(data[i]))[stop_col_num[j] : \
                                                              -len(end_dtypes)]
                    else:
                        column_data[name][i * nstops:(i + 1) * nstops] = \
                                                                  data[name][i]
            # Remove all zeros
            column_data = column_data[column_data[count_name] != 0]
            columnar_data.append(column_data)
        self.Columnar_Object = Columnar_Data(columnar_data, archival=archival)

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
            ff.output_form(data, filenames[i]) 


def remove_char(grid_list, char='\n'):
    '''
    Removes the given char from the end of each cell in grid list
    '''

    for grid in grid_list:
        for name in grid.dtype.names:
            for i in xrange(len(grid[name])): 
                while grid[name][i][::-1].find(char) == 0:
                    grid[name][i] = grid[name][i][:-1]
    
    return grid_list

def remove_white_spaces(grid_list):
    '''
    Removes all of the white spaces from strings.
    '''
    for grid in grid_list:
        for name in grid.dtype.names:
            for i in xrange(len(grid[name])):
                grid[name][i] = ''.join(grid[name][i].split(' '))

    return grid_list

def replace_vals(filename, replace, delim=','):
    '''
    Replace the values in filename with specified values in replace_values
    
    Parameters
    ----------
    filename : string 
        Will be read into a rec array

    replace_values : tuple
        First object is value to replace and second object is what to replace
        it with

    
    '''
    data = csv2rec(filename, delimiter=delim, missing=replace[0])
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
    return data






                     

            



            



















       









            





