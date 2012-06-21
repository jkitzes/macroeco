#!/usr/bin/python

'''This module contains 5 separate classes, each built to handle a
canonical data type'''


import numpy as np
from matplotlib.mlab import csv2rec
from form_func import *

__author__ = "Mark Wilber"
__copyright__ = "Copyright 2012, Regents of University of California"
__credits__ = "John Harte"
__license__ = None
__version__ = "0.1"
__maintainer__ = "Mark Wilber"
__email__ = "mqw@berkeley.edu"
__status__ = "Development"

class Census_Data:
    '''
    This class handles data types like the Smithsonian Research Plots.
    Examples include BCIS, LUQU, COCO, SHER.

    Multiple data files must have same format and missing value codes
    '''

    def __init__(self, datanames, delimiter=',', missingd=None,\
                delete_missing=True):
        '''
        This __init__ method takes in data and stores it in rec_arrays.
        IF specified,  it will located missing data points and remove them
        from the data set.

        Parameters
        ----------
        datanames : string or list of strings
            Data filenames

        missingd : dict
            Dictionary mapping munged column names to field values which 
            signify that the field does not contain actual data and should be
            masked, e.g. '0000-00-00' or 'unused'
        
        delete_missing : bool
            If True, deletes all of the missing values. If False, only deletes
            the NaNs from the data.

        '''
        self.data_list = []
        self.data_names = []
        if type(datanames) == str:
            self.data_names.append(datanames)
            self.data_list.append(csv2rec(datanames, delimiter=delimiter,\
                                  missingd=missingd))
        elif type(datanames) == list:
            for file_name in datanames:
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

    def split_up_data_by_field(self):
        pass

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

    def merge_data(self):
        '''
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










            





