#!/usr/bin/python

'''
This script interacts with the user to allow them to format their data into a
form useable by the macroeco software.  There data must be in one of the 4
canonical forms.

'''


from macroeco.utility import format_data
import sys

__author__ = "Mark Wilber"
__copyright__ = "Copyright 2012, Regents of University of California"
__credits__ = "John Harte"
__license__ = None
__version__ = "0.1"
__maintainer__ = "Mark Wilber"
__email__ = "mqw@berkeley.edu"
__status__ = "Development"

#The user may want to format multiple datasets, but for now we are only going
#let them format one dataset.  Though this can include multiple data files.

def format_columnar(filenames):
    '''
    '''
    pass

def format_grid(filenames):
    pass

def format_dense(filenames):
    pass

def format_transect(filenames):
    pass
    

if len(sys.argv) == 1:
    print "No dataset(s) found"
else:
    #If there are multiple data filenames should we assume they are all
    #related
    filenames = [names for names in sys.arv[1:]]
    dataset_type = raw_input("What is the format of your dataset? Your " + \
                            " options are 'columnar', 'grid', 'dense',  or" + \
                            " 'transect' (see documentaion for a" + \
                            " description of each format): ")
    while not(dataset_type is 'columnar' or dataset_type is 'grid' or dataset_type\
        is 'dense' or dataset_type is 'transect'):
        print "Incorrect dataset type"
        dataset_type = raw_input("What is the format of your dataset? Your " + \
                            " options are 'columnar', 'grid', 'dense',  or" + \
                            " 'transect' (see documentaion for a" + \
                            " description of each format): ")
        func_dict = {'columnar': format_columnar, 'grid' : format_grid,\
                     'dense' : format_dense, 'transect' : format_transect}





    




