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

import macroeco.utils.global_strings as gb 

gui_name = '''Rarity Analysis'''

summary = '''Compares a dataset's observed rarity against predicted rarity'''


predicted_SAD_distributions = '''This parameter is a list of SAD
distributions that you can test against your observed rarity.

You may use any number of the following SAD distributions : %s

Example input: ['logser', 'plognorm_lt'] or ['nbd_lt']. The brackets MUST be
included.''' % (gb.SAD_distributions)

rarity_measure = gb.rarity_measure + ''' In this analysis, the rarity
counts refer to individuals per species.'''

explanation = '''
ANALYSIS EXPLANATION\n
This script allows you to compare the rarity in your observed species abundance
distributions (SAD) to the rarity predicted by a predicted SAD. An SAD is a
distribution of the number of individuals within each species for an entire
community.  If you were to fully census a community with S species and N
individuals, the N individuals would be distributed among the S species in a
certain way. This analysis looks at both observed and predicted SADs and
determines how many species have the number of individuals less than or equal
to the rarity value(s) that you specify in the rarity parameter. For more
information on SADs and rarity please see the reference and the references
therein.

OUTPUT

This script outputs csv files containing the headings 'data_name', 'criteria',
and 'observed'.  The remainder of the headings are the names of the distributions to
which you are comparing your observed rarity.  The column 'data_name' contains
the name of the data that you are examining, the column 'criteria' specifies
the criteria that you imposed upon the data (i.e. how you split up the data in
the 'criteria' parameter), and the rest of the columns are the observed and
predicted rarity measurements. The file name of the csv specifies which level
of rarity is being examined in the csv file.  For example, if the file name
contained the phrase 'rarity_<=_10' the csv file is looking at rarity less then
or equal to 10.

PARAMETER EXPLANATIONS

*** subset ***:

{0}

*** criteria ***:

{1}

*** rarity_measure ***:

{2}

*** predicted_SAD_distributions ***:

{3}


REFERENCES

Harte, J. 2011. Maximum Entropy and Ecology: A Theory of Abundance,
Distribution, and Energetics. Oxford University Press.
'''.format(gb.subset, gb.criteria, gb.rarity_measure,
predicted_SAD_distributions)

required_params = { 'criteria' : gb.short_criteria + gb.req,
        'predicted_SAD_distributions' : predicted_SAD_distributions + gb.req,
        'rarity_measure' : rarity_measure + gb.req}

optional_params = {'subset' : (gb.short_subset + gb.optional, {})}


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
        for optpar in optional_params: #TODO: Move to workflow
            if not optpar in params:
                logging.info("Default value for {!s} in {!s}: {!s}".format(optpar,
                              output_ID,  str(optional_params[optpar][1])))
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
        mins = list(rarity['observed'].viewkeys())
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




        








