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

class global_str:
    subset = '''You should examine the columns in your data set and decide if you
	would like to subset your data in some particular way before the analysis
	begins. It is important to note that only the subsetted data will be analyzed.
	For example,  if you have a column named 'year' in your data set with values
	1998, 1999, and 2000 and you only want to look at the year 2000 for a
	particular analysis, you should select the == operator from the drop down list
	and type 2000 in the value field.  Similarly, you could use <, >, <=, >=, or
	!='''

    criteria = '''You should examine the columns in your dataset and decide if you
	would like to divide the data in a particular way for this analysis. For
	example, if the you have a spatial dataset with x,y coordinates and you are
	interested in examining macroecological metrics for two separate halves of your
	plot along the x coordinate, you could cut the x coordinate in two halves by
	giving the 'x' column a value of 2.  If the column that you would like to
	divide contains discrete values (e.g. year), you could enter the keyword
	'split' and each unique value will be analyzed separately. Conversely, the
	value 'whole' could be given to specify the entire column.  The value 'whole'
	is equivalent to 1 or leaving the value blank.\n\n

	There are four special words that can be used on a given column: 'species',
	'energy', 'count', and 'mass'.  When assigned to a column in your data set, the
	special word 'species' indicates the column that contains your species IDs, the
	special word 'energy' indicates the column that contains some type of energy
	measure, the special word 'mass' indicates a column that contains some type of
	mass measure, and the special word 'count' indicates the column that contains
	your species counts.  In the GUI, these special words can be chosen from the
	dropdown menu next to each column header. The special word 'species' MUST be
	assigned for every analysis.  If the special word 'count' is not assigned, the
	species counts are all assumed to be one.\n\n'''

    rarity_measure = '''This parameter allows you to specify the counts that
	you will consider rare.  If, for example, you want to know how many species in
	your plot have an abundance of 2 or less you would set this parameter to 2. If
	you enter more then one value, each value will be examined. Example input: [2]
	or [2, 5]. The brackets MUST be included.'''

    SAD_distributions = ''' 'logser','logser_ut', 'logser_ut_appx', 'plognorm_lt',
'nbd_lt', 'geo_ser', 'broken_stick', 'lognorm' '''

    SSAD_distributions = ''' 'nbd', 'binm', 'tgeo', 'fgeo', 'fnbd', 'pois' '''

gui_name = '''Rarity Analysis'''

summary = '''Compares a dataset's observed rarity against predicted rarity'''

explanation = '''This script allows you to compare the rarity in your observed
species abundance distributions (SAD) to the rarity predicted by a predicted
SAD. An SAD is a distribution of the number of individuals within each species
for an entire community.  If you were to fully census a community with S
species and N individuals, the N individuals would be distributed among the
S species in a certain way. 

This script outputs csv files containing the headings 'data_name', 'criteria',
and 'obs'.  The remainder of the headings are the names of the distributions to
which you are comparing your observed rarity.  The column 'data_name' contains
the name of the data that you are examining, the column 'criteria' specifies
the criteria that you imposed upon the data (i.e. how you split up the data in
the 'criteria' parameter), and the rest of the columns are the observed and
predicted rarity measurements. The file name of the csv specifies which level
of rarity is being examined in the csv file.  For example, if the file name
contained the phrase 'rarity_<=_10' the csv file is looking at rarity less then
or equal to 10.

For more information on SADs and rarity please see the following reference and
the references therein:

Harte, J. 2011. Maximum Entropy and Ecology: A Theory of Abundance,
Distribution, and Energetics. Oxford University Press.


'''

subset = global_str.subset

criteria = global_str.criteria

predicted_SAD_distributions = '''This parameter is a list of SAD
distributions that you can test against your observed rarity.

You may use any number of the following SAD distributions : %s

Example input: ['logser', 'plognorm_lt'] or ['nbd_lt']. The brackets MUST be
included.''' % (global_str.SAD_distributions)

rarity_measure = global_str.rarity_measure 

required_params = { 'criteria' : criteria,
        'predicted_SAD_distributions' : predicted_SAD_distributions,
        'rarity_measure' : rarity_measure}

optional_params = {'subset' : ('''Dictionary of initial subsets. Optional.
                                Default: ''', {})}


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
	for optpar in optional_params: # TODO: move into Workflow
	    if not optpar in params:
                logging.info('''Default value for {!s}: {!s}'''.format(			        optpar, str(optional_params[optpar][1])))
		params[optpar] = optional_params[optpar][1]

        patch = Patch(data_path, subset=params['subset'])
        
        # NOTE: Only looks at rarity for SADs 
        patch_out = patch.sad(params['criteria'], clean=True)

        cmpr = comp.CompareDistribution(patch_out, 
                 params['predicted_SAD_distributions'], patch='sad')
        rarity = cmpr.compare_rarity(cmpr.compare_rads(),
                                            mins_list=params['rarity_measure'])
       
        # NOTE: We could wrap all this in an output function if we want to
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




        








