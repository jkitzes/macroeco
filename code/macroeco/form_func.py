#!/usr/bin/python
'''This module contains the functions for formatting data files'''


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


def get_files(filetype, num, direct):
    '''This functoin gets the .txt files from the 
    data directory /archival/BCIS and returns the 
    names of the .txt files in the directory.  It is
    assumed that the .txt files are BCIS data.

    filetype -- a string specifying the type of the file, i.e. 'csv' or 'txt'

    num -- expected number of files of type 'direct_????.filetype'

    direct -- the directory within /data/archival/ where the files are.
              example 'BCIS' of 'COCO'
    
    returns:
        A list of strings
    
    '''

    #NOTE:This function is vulnerable to break! More checking needed
    assert direct.find('/') == -1, "%s should not contain a '/'" % (direct)
    cwd = gcwd();
    filedir = jp(pd(pd(gcwd())), 'archival', direct)
    chdir(filedir)
    datafiles = glob.glob(direct + '_????.' + filetype)
    chdir(cwd)
    if not(len(datafiles) == num):
        raise Exception("Must be exactly {0} {1}_*.{2} file in /archival/{1}"\
                        .format(num, direct, filetype))     
    return datafiles


def open_data(filename, delim):
    '''This functions takes in the data and returns a rec array.

    filename -- name of the file, a string

    delim -- file delimiter, a string

    '''

    data = plt.csv2rec(filename, delimiter=delim)
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
                        ('spp', 'S6')])
    spec_dict['spp'] = unq_specs
    spec_dict['spp_code'] = unq_ints
    return spec_dict

def format_by_year(tot_int, data):
    '''This function breaks data up by year, 
        removes dead trees, and flips y values.
        
        tot_int -- 1D np.array of integer codes for each species location
        
        data -- 2D np.array formatted like COCO or SHER

        returns:
            A list of formatted 2D np.arrays
               
        '''
    
    datalist = []
    #Eliminating all trees that aren't marked as alive
    #xrange(3) is 3 because this function only looks at 3 years
    for i in xrange(3):
        datalist.append(rm_trees_reduce(tot_int, data, i + 1))
      
    return datalist

def rm_trees_reduce(tot_int, data, year):
    '''This function removes all trees that are not
    marked as alive in the SHER_t data set and returns
    the reduced data for a given year.
        
    tot_int -- 1D np.array of integer codes for each species location
        
    data -- 2D np.array of COCO or SHER data
        
    year -- 1 = 1996, 2 = 1998, 3 = 1999 (int)
        
    returns:
        A rec array containing the data for a given year
        
        '''
    include = np.bitwise_or((data['recr' + str(year)] == 'A'),\
                          (data['recr' + str(year)] == 'P'))
                            
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
                        dtype=[('spp_code', np.int), ('x', np.int),
                               ('y', np.int), ('count', np.int)])

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


