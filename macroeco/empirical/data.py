#!/usr/bin/python

'''
Routines for loading census data and metadata.

Classes
-------
- `DataTable` -- data and metadata for a single censused area

DataTable Methods
-------------
- `dataload` -- loads census data from csv and metadata from EML xml file
- `get_metadata` -- load metadata
- `get_sub_table` -- return subset of census table
'''

from __future__ import division
import warnings
import numpy as np
from matplotlib.mlab import csv2rec
from macroeco.utility import metadata


class DataTable():
    '''
    Class to hold data table and metadata.

    Parameters
    ----------
    datapath : str
        Path to data. Location of metadata determined from this path.
    params : dict
        Params dictionary.

    Attributes
    ----------
    table : recarray
        Census data loaded from csv
    meta_dict : dict
        Dictionary of metadata associated with table.
    '''

    def __init__(self, datapath, params):
        '''Initialize XYTable object. See class docstring.'''
        self.table, self.meta = self.dataload(datapath, params)

    def dataload(self, datapath, params):
        '''
        Load table and metadata from files.
        
        Parameters
        ----------
        datapath : str
            Path to data table file.
        params : dict
            Params dictionary.
            
        Returns
        -------
        table : recarray
            Census data loaded from csv.
        meta : dict
            Dictionary of metadata associated with table.
        '''
        
        # Check that file is .csv.
        assert datapath[-3:] == 'csv', 'File must be .csv'

        # Load main table
        table = csv2rec(datapath)

        # Load metadata from file, printing warning and returning dict of zeros 
        # if file does not exist.
        metadatapath = datapath[:-3] + 'xml'

        asklist = []
        for name in table.dtype.names:
            asklist.append((name, 'minimum'))
            asklist.append((name, 'maximum'))
            asklist.append((name, 'precision'))

        try:
            meta_dict = self.get_metadata(asklist, metadatapath)
        except IOError:
            meta_dict = {}
            for id in asklist:
                meta_dict[id] = None
            warnings.warn('No metadata file found.')

        return table, meta_dict

    def get_metadata(self, asklist, metadatapath):
        '''
        Gets the relavent metadata from metadata file.

        Parameters
        ----------
        asklist : list
            List of tuples of column name and attribute, e.g., [('x', 
            'precision'), ('y', 'maximum')].

        metadatapath : str
            Path to the metadata file.
        '''

        meta_raw = metadata.Metadata(metadatapath)
        meta_raw.get_dataTable_values(asklist)
        meta_raw_dict = meta_raw.TableDescr

        meta_dict = {}
        for key, value in meta_raw_dict.iteritems():
            if value is None:
                new_value = value
            else:
                new_value = float(value)

            meta_dict[key] = new_value

        return meta_dict

    def get_sub_table(self, x_st, x_en, y_st, y_en):
        '''
        Returns subtable including only species within rectangle defined by 
        x_st, x_en, y_st, and y_en.
        
        Rectangle is inclusive only on the lower and left sides.
        '''

        x_col = self.head.index('x')
        y_col = self.head.index('y')

        in_sub = np.all([self.table[:,x_col] >= x_st, self.table[:,x_col] < 
                         x_en, self.table[:,y_col] >= y_st, self.table[:,y_col] 
                         < y_en], axis = 0)
        sub_table = self.table[in_sub]
        return sub_table
