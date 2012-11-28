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
The analysis compare_rarity allows you to compare the rarity in your observed
species abundance distributions (SAD) to the rarity predicted by a predicted
SAD. An SAD is a distribution of the number of individuals within each species
for an entire community.  If you were to fully census a community with S
species and N individuals, the N individuals would be distributed among the S
species in a certain way. This analysis looks at both observed and predicted
SADs and determines how many species have the number of individuals less than
or equal to the rarity value(s) that you specify in the rarity parameter. For
more information on SADs and rarity please see the reference and the references
therein.

OUTPUT

The analysis compare_rarity outputs one folder per dataset, a logfile.txt, and,
if possible, a map of the location(s) of the dataset(s). The folder has the
name rarity_values_compare_rarity_* and contains csv files containing the
headings 'data_name', 'criteria', and 'observed'.  The remainder of the
headings are the names of the distributions to which you are comparing your
observed rarity.  The column 'data_name' contains the name of the data that you
are examining, the column 'criteria' specifies the criteria that you imposed
upon the data (i.e. how you split up the data in the 'criteria' parameter), and
the rest of the columns are the observed and predicted rarity measurements. The
file name of the csv specifies which level of rarity is being examined in the
csv file.  For example, if the file name contained the phrase 'rarity_<=_10'
the csv file is looking at rarity less then or equal to 10.

The logfile.txt contains the analysis process information. Please see the
logfile if the analysis fails.

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

required_params = { 'criteria' : gb.req + gb.short_criteria,
        'predicted_SAD_distributions' : gb.req + predicted_SAD_distributions,
        'rarity_measure' : gb.req + rarity_measure}

optional_params = {'subset' : (gb.short_subset + gb.optional, {})}


if __name__ == '__main__':

    import logging
    from macroeco.utils.workflow import Workflow
    from macroeco.empirical import Patch
    import macroeco.compare as comp
    from macroeco.output import OutputRarity, make_directory
    import os

    wf = Workflow(required_params=required_params,
                optional_params=optional_params, clog=True, svers=__version__)

    folder_name = 'Rarity_analysis'
    make_directory(folder_name)
    cwd = os.getcwd()
    
    for data_path, output_ID, params in wf.single_datasets():

        os.chdir(os.path.join(cwd,folder_name))

        patch = Patch(data_path, subset=params['subset'])
        
        # NOTE: Only looks at rarity for SADs 
        patch_out = patch.sad(params['criteria'], clean=True)

        cmpr = comp.CompareSAD(patch_out, 
                 params['predicted_SAD_distributions'], patch=True)
        rarity = cmpr.compare_rarity(mins_list=params['rarity_measure'])
       
        OutputRarity(output_ID).output_rarity(rarity, data_path,
                                    cmpr.observed_data, criteria=cmpr.criteria)
    
        os.chdir(cwd)
        logging.info('Completed analysis %s\n' % output_ID)


    os.chdir(os.path.join(cwd,folder_name))
    fout = open('README_compare_rarity', 'w')
    with fout:
        fout.write(explanation)
    os.chdir(cwd)

    logging.info("Completed 'compare_rarity.py' script")




        








