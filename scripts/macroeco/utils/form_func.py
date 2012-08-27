#!/usr/bin/python
'''This module contains the functions for formatting data files'''

import os
import numpy as np
import csv
import matplotlib.mlab as plt
import glob
import sys

#Hacking this..Oh well
import format_data
loc = format_data.__file__
gcwd = os.getcwd #get current directory
pd = os.path.dirname #get parent directory
chdir = os.chdir #change directories
jp = os.path.join #Join paths
sys.path.append(pd(pd(loc)))
from data import Metadata

__author__ = "Mark Wilber"
__copyright__ = "Copyright 2012, Regents of University of California"
__credits__ = "John Harte"
__license__ = None
__version__ = "0.1"
__maintainer__ = "Mark Wilber"
__email__ = "mqw@berkeley.edu"
__status__ = "Development"

#Formatting functions
def get_metadata(asklist, folder_name, dataname):
    '''
    This function takes in a list of tuples and returns the appropriate
    metadata in a dictionary

    Parameters
    ----------
    asklist : list
        A list of tuples e.g. [('x', 'precision'), ('y', 'maximum')]

    folder_name : string
        Name of the archival folder where data is located e.g. BCIS

    dataname : string
        Name of the metadata e.g. BCIS_1984.xml (string)

    Returns
    -------
    : dict
        A dictionary containing requested metadata values
    
    '''
    cwd = gcwd()
    chdir(jp(pd(pd(gcwd())), 'archival', folder_name))
    meta = Metadata(dataname, asklist)
    chdir(cwd)
    return meta.get_meta_dict(asklist)

def get_files(filetype, num, direct, globber='_????'):
    '''
    This function gets the filetype files from the data directory
    /archival/direct and returns the names of the filetype files in the 
    directory.

    Parameters
    ----------
    filetype : string
        A string specifying the type of the file, i.e. 'csv' or 'txt'

    num : int
        Expected number of files of type 'direct_????.filetype'

    direct : string 
        The directory within /data/archival/ where the files are. 
        Example 'BCIS' or 'COCO'

    globber : string
        String of what pattern is to be globbed
    
    Returns
    -------
    : list
        A list of strings
    
    '''

    assert direct.find('/') == -1, "%s should not contain a '/'" % (direct)
    cwd = gcwd();
    filedir = jp(pd(pd(gcwd())), 'archival', direct)
    chdir(filedir)
    datafiles = glob.glob(direct + globber + '.' + filetype)
    chdir(cwd)
    if not(len(datafiles) == num):
        raise Exception("Must be exactly {0} {1}_*.{2} file in /archival/{1}"\
                        .format(num, direct, filetype))     
    return datafiles


def open_data(filename, delim, names=None):
    '''
    This functions takes in the filename and returns a rec array.
    
    Parameters
    ----------
    filename : string
        Name of the data file

    delim : string
        File delimiter

    names : list
        A list of columns names. See csv2rec?

    Returns
    -------
    : recarray
        A recarray containing the data from the specified file name

    '''

    data = plt.csv2rec(filename, delimiter=delim, names=names)
    return data

def create_intcodes(speclist, unq_specs, unq_ints):
    '''This function converts each species code into 
    its corresponding integer value.
    
    speclist -- a 1D np.array which contains the occurrences 
    of the species within the plot
        
    unq_specs -- a 1D np.array of the unique species codes 
    within the plot
        
    unq_int -- a 1D np.array of unique integers referring to
    the unique species codes found within the plot
        
    returns: 
        A 1D np.array of integers that is equivalent to
        speclist
        
    '''
    assert len(speclist) > 0, "Species array cannot be empty"
    speclist = speclist.astype(unq_specs.dtype)
    tot_int = np.empty(len(speclist))
    for s in xrange(len(unq_specs)):
        check = (unq_specs[s] == speclist)
        for i in xrange(len(check)):
            if check[i]:
                tot_int[i] = unq_ints[s]
    return tot_int

def output_form(data, filename):
    '''This function writes data as a .csv into the current working directory

    Parameters
    ----------
    data : structures array
        An structured array containing the data to be output

    filename : string
        A string representing the name of the file to be output.

    '''
    savedir = jp(gcwd(), filename.split('.')[0] + '.csv')
    fout = csv.writer(open(savedir, 'w'), delimiter=',')
    fout.writerow(data.dtype.names)
    for i in xrange(len(data)):
        fout.writerow(data[i])

def open_dense_data(filenames, direct, delim=','):
    '''
    This function takes in a list of dense data file names, opens
    them and returns them as list of rec arrays.

    Parameters
    ----------

    filenames : list 
        A list of filenames

    direct : string
        The directory within data/archival/ where the files are.
        Example 'ANBO_2010' or 'LBRI'

    delim : string
        The default file delimiter is ','

    Returns
    -------
    : list
        A list of rec arrays

    '''
    assert direct.find('/') == -1, "%s should not contain a '/'" % (direct)
    filedir = jp(pd(pd(gcwd())), 'archival', direct)
    datayears = []
    for name in filenames:
        data = plt.csv2rec(jp(filedir, name), delimiter=delim)
        datayears.append(data)
    return datayears

def format_dense(datayears, spp_col, num_spp):
    '''
    This function takes a list of data.  This functions interates 
    through the list and formats each year of data and stores the 
    formatted data into a list containing all years of formatted data.

    Parameters
    ----------
    datayears : list
        A list of rec arrays containing all years of data


    spp_col : int
        The column in the dense array where the spp_names begin. 0 is the first
        column.

    num_spp : tuple
        Total number of species in plot. Each element in the tuple is the
        number of species in the corresponding rec array in data year.
        Therefore, len(num_spp) should equal len(datayears).

    Returns
    -------
    : list
        A list of formatted structured arrays.

    '''
    data_formatted = []
    for k, data in enumerate(datayears):
        ls = len(data.dtype.names[spp_col:spp_col + num_spp[k]])
        if len(data.dtype.names[:spp_col + num_spp[k]]) == \
                                                        len(data.dtype.names):
            dtype = data.dtype.descr[:spp_col] + [('spp', 'S22'), ('count',\
                                                                np.float)]
        else:
            dtype = data.dtype.descr[:spp_col] + data.dtype.descr[spp_col + \
                    num_spp[k]:] + [('spp', 'S22'), ('count', np.float)]

        data_out = np.empty(ls * len(data), dtype=dtype)

        for s, name in enumerate(data_out.dtype.names[:-2]):
            cnt = 0
            for i in xrange(len(data)):
                if s == 0:
                    data_out[name][cnt:(ls*(i+1))] = data[name][i]
                    data_out['spp'][cnt:(ls*(i+1))] = np.array\
                                                (data.dtype.names[spp_col:\
                                                spp_col + num_spp[k]])
                    data_out['count'][cnt:(ls*(i+1))] =\
                                    np.array(list(data[i]))[spp_col:spp_col +\
                                    num_spp[k]]
                    cnt = cnt + ls
                else:
                    data_out[name][cnt:(ls*(i+1))] = data[name][i]
                    cnt = cnt + ls
        #Remove all zeros, they are not needed
        data_out = data_out[data_out['count'] != 0]
        data_formatted.append(data_out)
    return data_formatted

def open_nan_data(filenames, missing_value, site, delim, col_labels):
    '''
    This function takes in the filenames with nans data file, removes any
    NaN values for the x and y coordinates and returns a rec array.

    Parameters
    ----------
    
    filename : list
        A list of filenames which point to data with missing values

    missing_value : string
        How a missing value is labeled in the data

    site : string 
        Site name. Ex. 'COCO' or 'BCIS'

    delim : string
        Delimiter for the files

    xylabels : tuple 
        Tuple with x and y column labels, i.e. ('gx', 'gy') or ('x', 'y')

    Returns
    -------
    : list
        list of recarrays

    '''
    #NOTE: Might need to get rid of some more NA fields
    datadir = jp(pd(pd(gcwd())), 'archival', site)
    datayears = []
    for name in filenames:
        data = plt.csv2rec(jp(datadir, name), delimiter=delim,\
                                                    missing=missing_value)
        for label in col_labels:
            notNaN = (False == np.isnan(data[label]))
            data = data[notNaN]
        datayears.append(data)

    return datayears

def fractionate(datayears, wid_len, step, col_names):
    '''
    This function takes in a list of formatted data yearsand converts the grid
    numbers into meter measurements. Forexample, LBRI is a 16x16 grid and each
    cell is labeled withintegers.  However, the length (and width) of a cell is
    0.5m. This function converts each integer cell number to the appropriate
    integer (i.e. for LBRI cell (2,2) (counting from 1) becomes cell
    (0.5, 0.5)).
    
    Parameters
    ----------
    datayears : list 
        A list of formatted structured arrays

    wid_len : tuple
        A tuple containing the width (x) in meters and length (y)
        in meters of the entire plot.

    step : tuple
        The step (or stride length) of the cell width and length 
        (tuple: (x_step, y_step)). It shouldbe given in terms of meters. Also,
        called precision.

    col_names : list 
        The col_names of the structured array that are to be fractionated.

    Returns
    -------
    : list
        A list of converted structured arrays

    Notes
    -----
    This function should be used on formatted data

    '''
    
    frct_array = []
    for data in datayears:
        for i, name in enumerate(col_names):
            nums = np.unique(data[name])
            frac = np.arange(0, wid_len[i], step=step[i])
            #Have to make sure I have the data right type
            ind = list(data.dtype.names).index(name)
            dt = data.dtype.descr
            dt[ind] = (name, 'f8')
            data = data.astype(dt)
            data[name] = create_intcodes(data[name], nums, frac)
        frct_array.append(data)
    return frct_array

def add_data_fields(data_list, fields, values):
    '''
    Add fields to data based on given names and values

    Parameters
    ----------
    data_list : list 
        List of data to which a field will be appended

    fields : list
        List of field names to be added

    values : list
        List of tuples corresponding to fields.  Must be same length as fields.
        These values are added to the new field.

    Returns 
    -------
    : list
        A list containing the structured arrays with the new fields appended

    '''
    #TODO: Check that data in data list share dtypes
    alt_data = []
    for i, data in enumerate(data_list):
        for j, name in enumerate(fields):
            data = add_field(data, [(name, 'S20')])
            data[name] = values[j][i]
        alt_data.append(data)
    return alt_data

def merge_formatted(data_form):
    '''
    Take in a list of formatted data an merge all data in
    the list.  The dtypes of the data in the list must
    be the same

    Parameters
    ----------
    data_form : list 
        List of formatted structured arrays (or rec_arrays)

    Returns
    -------
    : list
        A list containing one merged structured array

    '''
    if len(data_form) == 1:
        return np.array(data_form[0])
    else:
        merged = np.array(data_form[0])
        for i in xrange(1, len(data_form)):
            if merged.dtype != data_form[i].dtype:
                raise TypeError("dtypes of data do not match")
            merged = np.concatenate((merged, np.array(data_form[i])))
        return merged

def add_field(a, descr):
    '''
    Add field to structured array and return new array with empty field

    Parameters
    ----------
    a : structured array
        Orginial structured array
    descr : list 
        dtype of new field i.e. [('name', 'type')]
    
    Returns
    -------
    : structured array
        Structured array with field added
    
    '''

    if a.dtype.fields is None:
        raise ValueError, "'A' must be a structured numpy array"
    b = np.empty(a.shape, dtype=descr + a.dtype.descr)
    for name in a.dtype.names:
        b[name] = a[name]
    return b

    

    
    






