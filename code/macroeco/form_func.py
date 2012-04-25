#!/usr/bin/python
'''This module contains the functions for formatting data files'''
'''TESTING git'''


import os
import numpy as np
import csv
import matplotlib.mlab as plt
import glob
from metadata import *

gcwd = os.getcwd #get current directory
pd = os.path.dirname #get parent directory
chdir = os.chdir #change directories
jp = os.path.join #Join paths

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
    '''This function takes in a list of tuples and returns the appropriate
    metadata in a dictionary

    asklist -- A list of tuples e.g. [('x', 'precision'), ('y', 'maximum')]

    folder_name -- name of the archival folder where data is e.g. BCIS

    dataname -- name of the metadata e.g. BCIS_1984.xml (string)

    
    '''

    cwd = gcwd()
    chdir(jp(pd(pd(gcwd())), 'archival', folder_name))
    meta = Metadata(dataname)
    meta.get_dataTable_values(asklist)
    chdir(cwd)
    return meta.TableDescr


def get_files(filetype, num, direct, globber='_????.'):
    '''This function gets the filetype files from the 
    data directory /archival/direct and returns the 
    names of the filetype files in the directory.  
    filetype -- a string specifying the type of the file, i.e. 'csv' or 'txt'

    num -- expected number of files of type 'direct_????.filetype'

    direct -- the directory within /data/archival/ where the files are.
              example 'BCIS' or 'COCO'

    globber -- String of what pattern is to be globbed (string)
    
    returns:
        A list of strings
    
    '''

    #NOTE:This function is vulnerable to break! More checking needed
    assert direct.find('/') == -1, "%s should not contain a '/'" % (direct)
    cwd = gcwd();
    filedir = jp(pd(pd(gcwd())), 'archival', direct)
    chdir(filedir)
    datafiles = glob.glob(direct + globber + filetype)
    chdir(cwd)
    if not(len(datafiles) == num):
        raise Exception("Must be exactly {0} {1}_*.{2} file in /archival/{1}"\
                        .format(num, direct, filetype))     
    return datafiles


def open_data(filename, delim, names=None):
    '''This functions takes in the data and returns a rec array.

    filename -- name of the file, a string

    delim -- file delimiter, a string

    '''

    data = plt.csv2rec(filename, delimiter=delim, names=names)
    return data

def make_spec_dict(spp_array):
    '''This function takes in an array of all the 
    species occurences in a plot and returns a species 
    dictionary with spp_codes coressponding to species
    names.

    spp_array -- 1D array of all species occurences

    returns:
        2D structured np.array

    '''
    assert len(spp_array) > 0, "Species array cannot be empty"
    unq_specs = np.unique(spp_array)
    unq_ints = np.linspace(0, len(unq_specs) - 1, num=len(unq_specs)) 
    spec_dict = np.empty(len(unq_ints), dtype=[('spp_code', np.int), \
                        ('spp', 'S40')])
    spec_dict['spp'] = unq_specs
    spec_dict['spp_code'] = unq_ints
    return spec_dict

def format_by_year(tot_int, data, col_name='recr', years=3):
    '''This function breaks data up by year, 
        removes dead trees, and flips y values.
        
        tot_int -- 1D np.array of integer codes for each species location
        
        data -- 2D np.array formatted like COCO or SHER

        col_name -- string, a column name in the data to use for formatting
        
        years -- an int specifying how many years of data are held in the
                 parameter data. col_name allows you to extract these years.

        returns:
            A list of formatted 2D np.arrays
               
        '''
    
    datalist = []
    #Eliminating all trees that aren't marked as alive
    #xrange(3) is 3 because this function only looks at 3 years
    for i in xrange(years):
        datalist.append(rm_trees_reduce(tot_int, data, i + 1, col_name))
      
    return datalist

def rm_trees_reduce(tot_int, data, year, col_name):
    '''This function removes all trees that are not
    marked as alive in the SHER_t data set and returns
    the reduced data for a given year.
        
    tot_int -- 1D np.array of integer codes for each species location
        
    data -- 2D np.array of COCO or SHER data
        
    year -- int

    returns:
        A rec array containing the data for a given year
        
    '''
    if col_name == 'recr':
        include = np.bitwise_or((data['recr' + str(year)] == 'A'),\
                                (data['recr' + str(year)] == 'P'))
    if col_name == 'dbh':
        include = (data['dbh' + str(year)] >= 1)
                            
    alive = data[include]
    tot_int_alive = tot_int[include]
    sample = len(tot_int_alive)
    data_out = np.empty(sample, dtype=[('spp_code', np.int), ('x', np.float), \
                                       ('y', np.float)])
    data_out['spp_code'] = tot_int_alive
    data_out['x'] = alive['x']
    data_out['y'] = alive['y']
    
    return data_out



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
    '''This function writes the formatted data into
    the appropriate formatted folder as a .csv file.

    data -- An np structured array containing the the data
            to be output

    filename -- A string representing the name of the file to 
                be output.

    Notes: Must be called within the appropriate formatted directory
    '''
    savedir = jp(gcwd(), filename.split('.')[0] + '.csv')
    fout = csv.writer(open(savedir, 'w'), delimiter=',')
    fout.writerow(data.dtype.names)
    for i in xrange(len(data)):
        fout.writerow(data[i])



def open_dense_data(filenames, direct, delim=','):
    '''This function takes in a list of dense data file names, opens
    them and returns them as list of rec arrays.

    filenames -- a list of filenames

    direct -- The directory within data/archival/ where the files are.
              example 'ANBO_2010' or 'LBRI'

    delim -- Default file delimiter is ','

    returns:
        A list of rec arrays

    '''
    assert direct.find('/') == -1, "%s should not contain a '/'" % (direct)
    filedir = jp(pd(pd(gcwd())), 'archival', direct)
    datayears = []
    for name in filenames:
        data = plt.csv2rec(jp(filedir, name), delimiter=delim)
        datayears.append(data)
    return datayears

def make_dense_spec_dict(datayears):
    '''This function makes and returns a global species dictionary
    for a list of dense data.  It assigns each species in the 
    data a unique integer code. The dense file should be formatted
    like LBRI of ANBO for this to work.

    datayears -- list of rec arrays containing data

    returns:
        A 1D array of tuples AKA a structured array

    '''
    
    spp_array = np.empty(0)
    for data in datayears:
        spp = np.array(data.dtype.names[3:]) #species names
        spp_array = np.concatenate((spp_array, spp))
    spec_dict = make_spec_dict(spp_array)
    return spec_dict

def format_dense(datayears, spec_dict):
    '''This function takes a list of data from 
    as well as a global species dictionary.  This functions interates 
    through the list and formats each year of data and stores the 
    formatted data into a list containing all years of formatted data.

    datayears -- a list of rec arrays containing all years of data

    spec_dict -- a 1D structured array of tuples containing the name of 
                 each species and the corresponding integer.
    
    returns:
        A list of formatted structured arrays.

    NOTE: Do not use this function unless you are sure that your data 
    match the formatting 

    '''

    for data in datayears:
        if set(data.dtype.names[0:3]) != set(['cell', 'column', 'row']):
            raise Exception("Improper labeling of first 3 columns")

    data_formatted = []
    for data in datayears:
        ls = len(data.dtype.names[3:])
        data_out = np.empty(ls * len(data), \
                        dtype=[('spp_code', np.int), ('x', np.float),
                               ('y', np.float), ('count', np.int)])

        tot_int = create_intcodes(np.array(data.dtype.names[3:]), \
                                   spec_dict['spp'], spec_dict['spp_code'])
        cnt = 0 
        for i in xrange(len(data)):
            data_out['x'][cnt:(ls*(i+1))] = data['row'][i]
            data_out['y'][cnt:(ls*(i+1))] = data['column'][i]
            data_out['spp_code'][cnt:(ls*(i+1))] = tot_int
            data_out['count'][cnt:(ls*(i+1))] = np.array(list(data[i]))[3:]
            cnt = cnt + ls

        data_formatted.append(data_out)

    return data_formatted

def open_nan_data(filenames, missing_value, site, delim, xy_labels=('gx', 'gy')):
    '''This function takes in the filenames with nans data file, removes any
    NaN values for the x and y coordinates and returns a rec array.
    
    filename -- list of filenames on data with nans (list)

    missing_value -- How a missing value is labeled in the data (string)

    site -- Site name (string)

    delim -- delimiter for the files (string)

    xylabels -- tuple with x and y column labels, i.e. ('gx', 'gy') or ('x', 'y')
    
    returns:
        a rec array

    '''
    #NOTE: Might need to get rid of some more NA fields
    datadir = jp(pd(pd(gcwd())), 'archival', site)
    datayears = []
    for name in filenames:
        data = plt.csv2rec(jp(datadir, name), delimiter=delim, missing=missing_value)
        xnotNaN = (False == np.isnan(data[xy_labels[0]]))
        data_out = data[xnotNaN]
        ynotNaN = (False == np.isnan(data_out[xy_labels[1]]))
        data_out = data_out[ynotNaN]
        datayears.append(data_out)

    return datayears

def make_multiyear_spec_dict(datayears, sp_field):
    '''This functions makes and returns a global species dictionary across 
    all years of a multi year data set.  It assigns each species that occurs
    over the years of data a unique integer code.

    datayears -- list of recarrays containing data (list)

    sp_field -- field name in rec array for species column (string)

    returns:
        A structured array with field names 'spp' and 'spp_code'
    '''

    spp_array = np.empty(0)
    for data in datayears:
        spp_array = np.concatenate((spp_array, data[sp_field]))
    spec_dict = make_spec_dict(spp_array)
    return spec_dict

def fractionate(datayears, wid_len, step):
    '''This function takes in a list of formatted data years
    and converts the grid numbers into meter measurements. For
    example, LBRI is a 16x16 grid and each cell is labeled with
    integers.  However, the length (and width) of a cell is 0.5m.
    This function converts each integer cell number to the appropriate
    integer (i.e. for LBRI cell (2,2) becomes cell (0.5, 0.5)). 
    
    NOTE: This function should be used on formatted data.

    datayears -- a list of formatted structured arrays

    wid_len -- a tuple containing the width (x) in meters and length (y)
               in meters of the entire plot.

    step -- the step (or stride length) of the cell width and
                    length (tuple: (x_step, y_step)). It should
                    be given in terms of meters.

    returns:
        a list of converted structured arrays

    '''
    
    for data in datayears:
        assert set(data.dtype.names[1:3]) == set(['x', 'y']), "Data must be\
                                            properly formatted"

    for data in datayears:
        x_nums = np.unique(data['x'])
        y_nums = np.unique(data['y'])
        x_frac = np.arange(0, wid_len[0], step=step[0])
        y_frac = np.arange(0, wid_len[1], step=step[1])
        data['x'] = create_intcodes(data['x'], x_nums, x_frac)
        data['y'] = create_intcodes(data['y'], y_nums, y_frac)
    return datayears




