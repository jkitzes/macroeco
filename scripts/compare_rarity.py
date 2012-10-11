#!/usr/bin/python

'''
Script to compare observed and predicted rarity
'''

__author__ = "Mark Wilber"
__copyright__ = "Copyright 2012, Regents of the University of California"
__credits__ = ["John Harte"]
__license__ = None
__version__ = "0.1"
__maintainer__ = "Mark Wilber"
__email__ = "mqw@berkeley.edu"
__status__ = "Development"

gui_name = '''Rarity Analysis'''

summary = '''Compares a dataset's observed rarity against theoretical rarity'''

explantion = '''This script takes in a dataset(s) and a list of distributions
and examines the observed and predicted rarity for each dataset against each
distribution.  Rarity is defined by the parameter 'rarity_measure'.
The required parameters for this script are the following: 

'subset' : How one would like to initially subset his data (see DataTable 
class docstring). 

'criteria' : How one would like to divide her dataset when caluculating the 
rarity. The theoretical distributions are compared against every dataset that is 
generated from the cutting specified in this parameter.  See Patch class 
docstring for details

'dist_list' : The list of distributions to which one could compare her observed
data.  The full list of distributions can be found in distributions.py

'patch_type' :  Rarity can be compared between either sads or ssads.  This
parameter is a string of either 'sad' or 'ssad'.  Depending on which form
patch_type takes, sads or ssads will be examined.

'rarity_measure' : This parameter specifies what the user would like to
consider rare for a given analysis.  For example, setting this parameter to
[10, 1] would specify that the user would like to consider all counts less than
or equal to 10 as rare and all counts less than or equal to 1 as rare.

For each dataset, a csv file(s) is generated with multiple columns.  The first
column is always the criteria used split the data set and the remaining columns
are observed rarity and predicted rarity.  Each row in the csv file displays
the observed and predicted rarity for the given subset displayed in the
criteria column.  So, if one divides the original data set into two smaller
data sets, the resulting csv file will have two rows (not including the column
names).  The file name specifies what measure of rarity was used. The number of
outputed csv files is equal to the number of values in the rarity_measure
parameter.

'''

required_params = {'subset' : 'Dictionary of initial subsets', 'criteria' :
        'Dictionary of how to split the data', 'dist_list' : 'List of' +\
        'distributions to compare', 'patch_type' : 'Either sad or ssad',
        'rarity_measure' : 'A list of values to consider rare'}

if __name__ == '__main__':

    import logging
    from macroeco.utils.workflow import Workflow
    from macroeco.empirical import Patch
    import macroeco.compare as comp
    from macroeco.utils.form_func import output_form
    import numpy as np
    import os

    wf = Workflow(required_params=required_params, clog=True, 
                                                            svers=__version__)
    
    for data_path, output_ID, params in wf.single_datasets():

        patch = Patch(data_path, subset=params['subset'])

        if params['patch_type'] == 'sad':
            patch_out = patch.sad(params['criteria'], clean=True)
        elif params['patch_type'] == 'ssad':
            patch_out = patch.ssad(params['criteria'])
        else:
            raise ValueError('%s not a recognized patch type' %
                                                        params['patch_type'])

        cmpr = comp.CompareDistribution(patch_out, params['dist_list'], 
                                                   patch=params['patch_type'])
        rarity = cmpr.compare_rarity(cmpr.compare_rads(),
                                            mins_list=params['rarity_measure'])
        
        # Make dtype
        keys = list(rarity.viewkeys())
        dtype = [(kw, np.int) for kw in keys]
        max_len = np.max([len(str(crit)) for crit in cmpr.criteria])
        dtype.insert(0, ('criteria', 'S90')) # arbitrary length
        dtype.insert(0, ('data_name', 'S90')) # arbitrary length

        # Get a list of my minimums
        rare_list = []
        mins = list(rarity['obs'].viewkeys())
        for mn in mins:
            rarity_array = np.empty(len(cmpr.data_list), dtype=dtype)
            rarity_array['criteria'] = cmpr.criteria
            nm = os.path.split(data_path)[1].split('.')[0]
            rarity_array['data_name'] = np.repeat(nm, len(rarity_array))
            for kw in keys:
                rarity_array[kw] = rarity[kw][mn]
            rare_list.append(rarity_array)
        
        # Output results
        for i, rare in enumerate(rare_list):
            output_form(rare, output_ID + '_rarity_<=_' + str(mins[i]))
        
        logging.info('Completed analysis %s\n' % output_ID)
    logging.info("Completed 'compare_rarity.py' script")




        








