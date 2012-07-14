#!/usr/bin/python

'''This module contains 4 separate classes, each built to handle a
canonical data type

This module provides the user with some formatting functions but does provide
the user with all formatting functions that may be required.  This module is
not a substitute for thorough examination of ones data to remove irrelevant
data'''


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
    then turned into Columnar_Data obejcts.

    Examples of columnar data include BCIS, LUQU, and COCO

    Multiple data files must have same format if they are to be merged

    '''

    def __init__(self, datalist, delimiter=',', missingd=None,\
                delete_missing=True):
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

        '''
        if type(datalist) == str:
            datalist = [datalist]

        if np.all(np.array([type(x) == str for x in datalist])):
            self.data_list = []
            self.data_names = []
            for file_name in datalist:
                self.data_list.append(csv2rec(file_name, delimiter=delimiter,\
                              missingd=missingd))
                self.data_names.append(file_name)
            if missingd != None:
                if delete_missing:
                    for key in missingd.iterkeys():
                        for data in self.data_list:
                            notNaN = (False == np.isnan(data[key]))
                            notBlank = (data[key] != '')
                            ind = np.bitwise_or(notNaN, notBlank)
                            data = data[ind]
                else:
                    for key in missingd.iterkeys():
                        for data in self.data_list:
                            notNaN = (False == np.isnan(data[key]))
                            data = data[notNaN]
        elif np.all(np.array([type(x) == np.ndarray for x in datalist])):
            self.data_list = datalist
        self.archival_data = [np.copy(data) for data in self.data_list]

    def reset_data_list(self):
        '''
        Resets self.data_list to self.archival_data
        
        Need to be careful about excessive memory usage!
        '''
        self.data_list = [np.copy(data) for data in self.archival_data]

    def split_up_data_by_field(self, split_columns):
        '''
        This function will take in the split-columns list and and split the
        data into separate arrays based on the list.  For example, if one were
        to pass in dbh1, dbh2,  dbh3 the data three copies of the data would be
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
        Saves the split array as self.data_list
        
        
        '''
        #Note: If they enter the wrong column name nothing will be removed
        #Should I error check for this?
        split_data = []
        given_col_names = []
        for tup in split_columns:
            for name in tup:
                given_col_names.append(name)
        given_col_names = np.array(given_col_names)


        for data in self.data_list:
            for tup in split_columns:
                ind = np.ones(len(given_col_names), dtype=bool)
                for name in tup:
                    ind = np.bitwise_and((name != given_col_names), ind)
                remove_names = given_col_names[ind]
                split_data.append(drop_fields(data, list(remove_names)))
        self.data_list = split_data
    
    def change_column_names(self, change, changed_to):
        '''
        This function takes a list of column names to be changed and a name
        that they should be changed to

        Parameters
        ----------
        change : list
            List of column names.  Columns names are strings
        changed_to : string
            Name to be changed to

        '''
        
        for data in self.data_list:
            column_names = np.array(data.dtype.names)
            for name in change:
                find = np.where((name == column_names))[0]
                if len(find) != 0:
                    column_names[find[0]] = changed_to
                    data.dtype.names = tuple(column_names)
        
    def add_fields_to_data_list(self, fields, values):
        '''
        This functions adds given fields and values to the data list. The
        length of values should be the same length as fields and the length of
        each tuple within each element of values should be the same length as
        the self.data_list

        Parameters
        ----------
        fields : list
            A list of strings specifying field names
        values : list of tuples
            A list of tuples with the length of each tuple equalling the length
            of self.data_list

        '''
        #NOTE: Should probably make a single dictionary for field/values
        self.data_list = add_data_fields(self.data_list, fields, values)
    
    def fractionate_data(self, wid_len, step, col_names):
        '''
        fractionate_data(self, wid_len, step, col_names)

        This function converts grid numbers to length measurements in
        self.data_list

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
        self.data_list = fractionate(self.data_list, wid_len, step, col_names)


    def merge_data(self):
        '''
        This function concatenates the data files in data_list.  The dtypes of
        the data in data_list must be identical or this function will fail.
        '''

        self.merged_data = merge_formatted(self.data_list)

    def output_merged_data(self, filename):
        '''
        This function merges self.data_list and outputs the merged data.

        Parameters
        ----------
        filename : string
            The filename to be output

        '''
        #Merge data in case it has not been done
        self.merge_data()
        output_form(self.merged_data, filename)

    def output_data_list(self, filenames):
        '''
        This function outputs the self.data_list

        Parameters
        ----------
        filenames : list
            A list of filenames

        '''
        assert len(filenames) == len(self.data_list), "Number of filenames\
                                 must be the same as the number of datasets"
        for i, name in enumerate(filenames):
            output_form(self.data_list[i], name)

class Grid_Data:
    '''This class handles data should look like the EarthFlow data after a 
    census.  It is a grid with species abundance data in each cell. 
    ex.
    ARTDRA - 6
    GERTYR - 8

    NOTE: Need to consider things that might break this class
    '''

    def __init__(self, filenames, num_cols):
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

        '''
        #NOTE: Handle missing data!!!!

        self.dense_data = []

        if type(filenames) == str:
            filenames = [filenames]

        if type(num_cols) == int:
            num_cols = [num_cols]

        assert np.all(np.array([name.split('.')[-1] for name in filenames]) ==\
                      'csv'), "Files must be csv"
        assert len(num_cols) == len(filenames) or len(num_cols) == 1, 'Length\
                       of num_cols must be 1 or equal len(filenames)'

        self.grids = []
        self.cols = []
        self.rows =[]

        for i, name in enumerate(filenames):
            if len(num_cols) == 1:
                self.cols.append(num_cols[0])
            else:
                self.cols.append(num_cols[i])
            self.grids.append(csv2rec(name, names=list(np.arange(0,\
                                            self.cols[i]).astype('S10'))))
            self.rows.append(len(self.grids[i]))

    def find_unique_spp_in_grid(self, spacer='-', spp_sep='\n'):
        '''
        This function finds all of the unique species in the grid.
        It assumes that your grid data is in the proper format.

        Parameters
        ----------
        spacer : str
            The character separating the species code from the species count.
            Default value is '-' (n-slash)

        spp_sep : str
            The character that separates a speces/count combination from
            another species/count combination.  Default value is '\n'

        '''
        self.unq_spp_lists = []
        for num, data in enumerate(self.grids):
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

    def grid_to_dense(self, spacer='-', spp_sep='\n'):
        '''
        This function converts a the list of gridded data sets into dense 
        data sets and stores them in self.dense_data.  In addition, it
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
        self.dense_data = []
        for i, data in enumerate(self.grids):
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
                        start = data[col][row].find(spp_name)
                        if start != -1:
                            raw = data[col][row][start:].split('-')[1]
                            if raw.find(spp_sep) != -1:
                                tot_spp = raw.split(spp_sep)[0].strip()
                                matrix[spp_name][count] = int(tot_spp)
                            else:
                                tot_spp = raw.split()[0].strip()
                                matrix[spp_name][count] = int(tot_spp)
                        else:
                            matrix[spp_name][count] = 0
                    count += 1
            self.dense_data.append(matrix)
        self.Dense_Object = Dense_Data(self.dense_data)

    def output_data_matricies(self, filenames):
        '''
        This function prints the data within self.dense_data with the given
        filenames.  If self.dense_data has not been filled, error is thrown.


        Parameters
        ----------
        filenames : list
            A list of filenames to which the data will be saved
    
        '''
        if len(self.dense_data) == 0:
            raise Exception("No data in self.dense_data")

        assert len(filenames) == len(self.dense_data), "Number of filenames\
                                 must be the same as the number of datasets"
        for i, data in self.dense_data:
            output_form(data, filenames[i]) 

    
class Dense_Data:
    '''This class handles data that are in the dense format

    MORE DOC STRING

    Note: Need to consider how I might break this class
    '''

    def __init__(self, datalist, delim=','):
        '''
        Going to duck-type datalist to find out whether it is list of data or
        a list of strings!

        '''
        #TODO: What kind of files could break this
        #NOTE: What about missing values?
        if type(datalist) == str:
            datalist = [datalist]

        self.columnar_data = []
        if np.all(np.array([type(x) == str for x in datalist])):
            self.dense_data = []
            for name in datalist:
                self.dense_data.append(csv2rec(name, delimiter=delim))
        elif np.all(np.array([type(x) == np.ndarray for x in datalist])):
            self.dense_data = datalist

    def dense_to_columnar(self, spp_col_num, num_spp):
        '''
        This function uses a function in form_func to convert dense data into
        columnar data. Stores the columnar data as a Columnar Object.

        Parameters
        ----------
        spp_col_num : int
            The column number in the dense array where the spp_names begin

        num_spp : int
            Number of species in plot

        '''

        self.columnar_data = format_dense(self.dense_data, spp_col_num,\
                                                                      num_spp)
        self.Columnar_Object = Columnar_Data(self.columnar_data)

    def output_columnar_data(self, filenames):
        '''
        This function prints the data within self.columnar_data with the given
        filenames.  If self.columnar_data has not been filled, an error is 
        thrown.

        Parameters
        ----------
        filenames : list
            A list of filenames to which the data will be saved
        '''
        if len(self.columnar_data) == 0:
            raise Exception("No data in self.columnar_data")

        assert len(filenames) == len(self.columnar_data), "Number of filenames\
                                 must be the same as the number of datasets"
        for i, data in self.columnar_data:
            output_form(data, filenames[i]) 


class Transect_Data:
    '''
    This class handles data that are similar to the Breeding Bird survey data.
    One column has the species ID, one column has stop and all the other
    columns have transects.  This class can handle data with "n" nestings, not
    just two.  For example, the data could have location, transect and stop.

    The "stop" data should all be in consecutive columns

    '''
    
    def __init__(self, filenames, delimiter=','):
        '''
        '''
        self.columnar_data = []
        self.transect_data = []
        if type(filenames) == str:
            filenames = [filenames]

        for name in filenames:
            data = csv2rec(name, delimiter=delimiter)
            self.transect_data.append(data)

    def transect_to_columnar(self, stop_col_num, stop_name, tot_stops,\
                                                    count_name='count'):
        '''
        This function takes transect data and convertes it into columnar data.
        In addition it saves the columnar data as a Columnat_Data object. 
        

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
        This function assumes that all data in self.transect_data is formatted
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

    def output_columnar_data(self, filenames):
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
        if len(self.columnar_data) == 0:
            raise Exception("No data in self.columnar_data")

        assert len(filenames) == len(self.columnar_data), "Number of filenames\
                                 must be the same as the number of datasets"
        for i, data in self.columnar_data:
            output_form(data, filenames[i]) 



                     

            



            



















       









            





