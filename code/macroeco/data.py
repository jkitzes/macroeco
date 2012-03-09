#!/usr/bin/python

'''
Routines for loading and converting data.

This module provides functions that load raw macroecological data and convert 
it between forms. Currently, the only data option is class xytable.

Functions
---------
- `xytable_load` -- load xytable data and metadata from file
- `xytable_add_count` -- add count column to xytable if not already present
'''

from __future__ import division
import numpy as np


class XYTable():
    '''
    Class to hold xytable data and parameters.

    Parameters
    ----------
    datapath : str
        Path to data. Location of metadata determined from this path.
    compress : bool
        Flag to compress data so that there is only one row in xytable for each 
        combination of spp and point coordinate.

    Attributes
    ----------
    table : ndarray
        xytable data with cols for spp_code, x, y, etc. No header.
    head : list
        Headers for cols of table.
    meta : dict
        Dictionary of metadata associated with table.
    spp_codes : ndarray
        1D array of all specie codes in patch.
    max_spp_code : int
        Maximum value of species code found in xy.
    '''

    def __init__(self, datapath, compress = False):
        '''
        Initialize XYTable object. See class docstring.
        '''

        self.table, self.head, self.meta = self.xytable_load(datapath, 
                                                             compress)

        self.col_spp_code = self.head.index('spp_code')
        self.col_count = self.head.index('count')

        self.spp_codes = np.unique(self.table[:, self.col_spp_code])
        self.max_spp_code = max(self.spp_codes)


    def xytable_load(self, datapath, compress):
        '''
        Load xytable table and metadata from file.
        
        Parameters
        ----------
        datapath : str
            Path to xytable table file. File must have header row and include 
            cols, at a minimum, named int_code, x, and y.
        compress : bool
            Flag to compress data so that there is only one row in xytable for 
            each combination of spp and point coordinate.
            
        Returns
        -------
        xy_table : 2D ndarray
            Array of xytable data with cols specified in xy_head.
        xy_head : list
            List of headers for cols in xy_data.
        meta : dict
            Dictionary of metadata needed for analysis of xy_data. If no 
            metadata file found, returns None.
        '''
        
        # Check that file is .csv. Not strictly necessary (only required that
        # delimiter be a comma.
        if datapath[-3:] != 'csv':
            raise NotImplementedError('File must be .csv')

        # Attempt to load main table
        xy_table = np.loadtxt(datapath, delimiter = ',', skiprows = 1)

        # Attempt to load header of data file
        with open(datapath) as file:
            xy_head = file.next()
        xy_head = xy_head.replace(',','').split()  # Remove commas, make list

        # Check that there are cols in data named spp_code, x, y
        req_cols = set(['spp_code', 'x', 'y'])
        cols = set(xy_head)
        if cols.intersection(req_cols) != req_cols:
            raise NameError('Data missing required column.')

        # Add count col if not already present
        xy_table, xy_head = self.xytable_add_count(xy_table, xy_head, 
                                                   compress)

        # Load metadata from file, printing warning and returning None if file 
        # does not exist.
        metapath = datapath[:-3] + '.xml'
        try:
            np.loadtxt('not_implemented')  # TODO: Add reader implementation
        except IOError:
            print('No metadata file found, meta filled with 0\'s')
            meta = {'precision': 0, 'xrange': (0, 0), \
                    'yrange': (0, 0)}
            
        # Check meta has required keys
        req_keys = set(['precision', 'xrange', 'yrange'])
        keys = set(meta.keys())
        if keys.intersection(req_keys) != req_keys:
            raise NameError('Metadata missing required parameter.')

        return xy_table, xy_head, meta
        

    def xytable_add_count(self, xy_in_data, xy_in_head, compress):
        '''
        Add count column to xytable if not already present.
        
        Count gives the number of individuals of spp found at the point x, y. 
        Note that the option to compress the data is not yet implemented.

        Parameters
        ----------
        xy_in_data : 2D ndarray
            Array of xytable data with cols specified in xy_in_head
        xy_in_head : list
            List of headers for cols in xy_data
        compress : bool
            Flag to compress data so that there is only one row in xytable for 
            each combination of spp and point coordinate.

        Returns
        -------
        xy_out_data : 2D ndarray
            Array of updated xytable data with cols specified in xy_out_head
        xy_out_head : list
            List of updated headers (with 'count' appended) for cols in xy_data
        '''
        
        # Check if header named count, if so do nothing, else add 1's
        if 'count' not in xy_in_head:
            rows = xy_in_data.shape[0]
            xy_out_data = np.append(xy_in_data, np.ones((rows,1)), axis = 1)
            xy_out_head = xy_in_head + ['count']
        else:
            xy_out_data = xy_in_data
            xy_out_head = xy_in_head

        # If not compress, return xy_out_data, else throw not implemented error
        if not compress:
            return xy_out_data, xy_out_head
        else:
            raise NotImplementedError('Compression of xytable not yet ' + 
                                      'implemented')

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
