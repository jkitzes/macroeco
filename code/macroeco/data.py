#!/usr/bin/python

'''
Routines for loading and converting data.

This module provides functions that load raw macroecological data and convert 
it between forms.

Functions
---------
- `load_xytable` -- load xytable data and metadata from file
- `add_count_xytable` -- add count column to xytable if not already present
'''


import numpy as np


def load_xytable(filepath, param_dict):
    '''
    Load xytable data and metadata from file.
    
    Parameters
    ----------
    filepath : str
        Path to xytable data file. File must have header row and include cols, 
        at a minimum, named int_code, x, and y.
    param_dict : dict
        If exists, returned as dict of parameters for analysis of data. If not, 
        param_dict is loaded from metadata xml file.
        
    Returns
    -------
    xy_data : 2D ndarray
        Array of xytable data with cols specified in xy_head
    xy_head : list
        List of headers for cols in xy_data
    '''

    # TODO: Add loading and parsing metadata
    # TODO: Currently only takes in comma delimited, allow other options
    try:
        data = np.loadtxt(filename, delimiter = ',')
    except IOError as detail:
        print detail
    
    return data


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
        List of headers for cols in xy_in_data
    compress : bool
        Flag to compress data so that there is only one row in xytable for each 
        combination of spp and point coordinate.

    Returns
    -------
    xy_out_data : 2D ndarray
        Array of updated xytable data with cols specified in xy_out_head
    xy_out_head : list
        List of headers for cols in xy_out_data
    '''
    
    # Check if already header named count, if so do not add col, else add 1's

    # Check if compress - if not, do nothing, if so, 
    #

    if not compress:
        return xy_out_data, xy_out_head
    else:
        raise NotImplementedError('Compression of xytable not yet implemented')
