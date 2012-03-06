#!/usr/bin/python

'''
Routines for loading and converting data.

This module provides functions that load raw macroecological data and convert 
it between forms.

Functions
---------
- `xytable_load` -- load xytable data and metadata from file
- `xytable_add_count` -- add count column to xytable if not already present
'''


import numpy as np


def xytable_load(filepath, meta = {}):
    '''
    Load xytable data and metadata from file.
    
    Parameters
    ----------
    filepath : str
        Path to xytable data file. File must have header row and include cols, 
        at a minimum, named int_code, x, and y.
    meta : dict
        If exists, returned as dict of metadata needed for analysis of data. If 
        not, meta is loaded from metadata xml file.
        
    Returns
    -------
    xy_data : 2D ndarray
        Array of xytable data with cols specified in xy_head
    xy_head : list
        List of headers for cols in xy_data
    meta : dict
        Dictionary of metadata needed for analysis of xy_data
    '''
    
    # Check that file is .csv. Not strictly necessary (only required that
    # delimiter be a comma.
    if filepath[-3:] != 'csv':
        raise NotImplementedError('File must be .csv')

    # Attempt to load data
    xy_data = np.loadtxt(filepath, delimiter = ',', skiprows = 1)

    # Attempt to load header of data file
    with open(filepath) as file:
        xy_head = file.next()
    xy_head = xy_head.replace(',','').split()  # Remove commas, make into list

    # Check that there are cols in data named spp_code, x, y
    req_cols = set(['spp_code', 'x', 'y'])
    cols = set(xy_head)
    if cols.intersection(req_cols) != req_cols:
        raise NameError('Data missing required column.')

    # Add count col
    xy_data, xy_head = xytable_add_count(xy_data, xy_head, compress = False)

    # If meta not given as argument, load from file
    if not meta:
        raise NotImplementedError('Reading from metadata file not yet '+
                                  'implemented.')

    # Check meta has required keys
    req_keys = set(['precision', 'xrange', 'yrange'])
    keys = set(meta.keys())
    if keys.intersection(req_keys) != req_keys:
        raise NameError('Metadata missing required parameter.')

    return xy_data, xy_head, meta
    

def xytable_add_count(xy_in_data, xy_in_head, compress):
    '''
    Add count column to xytable if not already present.
    
    Count gives the number of individuals of spp found at the point x, y. Note 
    that the option to compress the data is not yet implemented.

    Parameters
    ----------
    xy_in_data : 2D ndarray
        Array of xytable data with cols specified in xy_in_head
    xy_in_head : list
        List of headers for cols in xy_data
    compress : bool
        Flag to compress data so that there is only one row in xytable for each 
        combination of spp and point coordinate.

    Returns
    -------
    xy_out_data : 2D ndarray
        Array of updated xytable data with cols specified in xy_out_head
    xy_out_head : list
        List of updated headers (with 'count' appended) for cols in xy_data
    '''
    
    # Check if already header named count, if so do not add col, else add 1's
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
        raise NotImplementedError('Compression of xytable not yet implemented')
