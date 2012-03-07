#!/usr/bin/python
'''This module contains the functions for formatting data files'''


import os
import numpy as np
import csv
import matplotlib.mlab as plt
import glob

__author__ = "Mark Wilber"
__copyright__ = "Copyright 2012, Regents of University of California"
__credits__ = "John Harte"
__license__ = None
__version__ = "0.1"
__maintainer__ = "Mark Wilber"
__email__ = "mqw@berkeley.edu"
__status__ = "Development"

#Formatting functions

def get_files(filetype, num, direct):
    '''This functoin gets the .txt files from the 
    data directory /archival/BCIS and returns the 
    names of the .txt files in the directory.  It is
    assumed that the .txt files are BCIS data.

    filetype -- a string specifying the type of the file, i.e. 'csv' or 'txt'

    num -- expected number of files of type 'direct_????.filetype'

    direct -- the directory within /data/archival/ where the files are.
              example 'BCIS'
    
    returns:
        A list of strings
    
    '''

    #NOTE:This function is vulnerable to break! More checking needed
    datadir = os.path.dirname((os.path.dirname(os.getcwd())))
    globdir = datadir + '/archival/' + direct
    datafiles = glob.glob(globdir + '/' + direct + '_????.' + filetype)
    datafiles = [x.split('/')[len(x.split('/')) - 1] for x in datafiles]
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

    '''
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
    
    alive = data[(data['recr' + str(year)] == 'A')]
    tot_int_alive = tot_int[(data['recr' + str(year)] == 'A')]
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
    
    speclist -- a 1D np.array which contains the occurences 
    of the species within the plot
        
    unq_specs -- a 1D np.array of the unique species codes 
    within the plot
        
    unq_int -- a 1D np.array of unique integers referring to
    the unique species codes found within the plot
        
    returns: 
        A 1D np.array of integers that is equivalent to
        speclist
        
    '''
    
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
    '''
    
    savedir = os.getcwd() + '/' + (filename.split('.')[0] + '.csv')
    fout = csv.writer(open(savedir, 'w'), delimiter=',')
    fout.writerow(data.dtype.names)
    for i in xrange(len(data)):
        fout.writerow(data[i])



