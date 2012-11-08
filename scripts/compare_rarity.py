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

summary = '''Compares a dataset's observed rarity against predicted rarity'''

class global_str:
    subset = '''You should examine the columns in your data set and decide if
    you would like to subset your data in some particular way before the
    analysis begins. It is important to note that only the subsetted data will
    be analyzed.  For example,  if you have a column named 'year' in your data
    set with values 1998, 1999, and 2000 and you only want to look at the year
    2000 for a particular analysis, you should select the column year from
    left-hand most dropdown list, select the == operator from the operator
    dropdown list and type 2000 in the value field.  Similarly, you could use
    <, >, <=, >=, or != with any column and value in your data.'''

    criteria = '''You should examine the columns in your dataset and decide if
    you would like to divide the data in a particular way for this analysis.
    For example, if you have a spatial dataset with x,y coordinates and you
    are interested in examining macroecological metrics for two separate halves
    of your plot along the x coordinate, you could cut the x coordinate in two
    halves by giving the 'x' column a value of 2.  If the column that you would
    like to divide contains discrete values (e.g. year), you could enter the
    keyword 'split' and each unique value will be analyzed separately.
    Conversely, the value 'whole' could be given to specify the entire column.
    The value 'whole' is equivalent to 1 or leaving the value blank. If you
    would like to divide a given column, please select the word 'division' from
    the dropdown menu and input a value as discussed above.\n\n

    There are four other special words that can be used on a given column:
    'species', 'energy', 'count', and 'mass'.  When assigned to a column in
    your data set, the special word 'species' indicates the column that
    contains your species IDs, the special word 'energy' indicates the column
    that contains some type of energy measure, the special word 'mass'
    indicates a column that contains some type of mass measure, and the special
    word 'count' indicates the column that contains your species counts.  These
    special words can be chosen from the dropdown menu next to each column
    header. The special word 'species' MUST be assigned for every analysis.  If
    the special word 'count' is not assigned, the species counts are all
    assumed to be one.\n\n
    
    If there are columns in your data that are not relevant for this analysis
    leave the value in the dropdown box as 'NA'.  Columns designated 'NA'
    will not influence the analysis.\n\n'''


    rarity_measure = '''This parameter allows you to specify the counts that
    you will consider rare.  If, for example, you want to know how many species
    in your plot have an abundance of 2 or less you would set this parameter to
    2. If you enter more then one value, each value will be examined. Example
    input: [2] or [2, 5]. The brackets MUST be included.'''

    SAD_distributions = ''' 
    'logser' : Fisher's logseries distribution;
    'logser_ut' : Upper-truncated logseries derived from MaxEnt;
    'logser_ut_appx' : Approximation for the upper-truncated logseries;
    'lognorm' : Lognormal distribution;
    'plognorm_lt' : Poisson lognormal distribution with 0 truncated;
    'nbd_lt' : Negative binomial distribution with 0 truncated;
    'geo_ser' : Geometric series distribution;
    'broken_stick' : McArthur's broken stick distribution '''

    SSAD_distributions = ''' 
    'binm' : Binomial distribution;
    'pois' : Poisson distribution;
    'nbd' : Negative binomial distribution;
    'fnbd' : Finite-negative binomial distribution;
    'geo' : Geometric distribution;
    'fgeo' : Finite-geometric distribution;
    'tgeo' : Truncated geometric distrbituion derived from MaxEnt'''


subset = '''Specifications for how you want to subset your data before the
analysis. Note that only the subsetted data will be included in the analysis.
The left-hand dropdown box contains all the columns of your dataset and you may
choose one or more to subset. Please see analysis explanation for more detail
and examples.'''

criteria = '''Specifications for how you want to divide your data during the
analysis. The words you see below are the shared columns of your dataset(s).
You must designate your species column with the special word 'species' found in
the dropdown menu. You are not required to fill any additional columns for this
analysis. Please see analysis explanation for more detail and examples.'''

predicted_SAD_distributions = '''This parameter is a list of SAD
distributions that you can test against your observed rarity.

You may use any number of the following SAD distributions : %s

Example input: ['logser', 'plognorm_lt'] or ['nbd_lt']. The brackets MUST be
included.''' % (global_str.SAD_distributions)

rarity_measure = global_str.rarity_measure + ''' In this analysis, the rarity
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
'''.format(global_str.subset, global_str.criteria, global_str.rarity_measure,
predicted_SAD_distributions)

required_params = { 'criteria' : criteria,
        'predicted_SAD_distributions' : predicted_SAD_distributions,
        'rarity_measure' : rarity_measure}

optional_params = {'subset' : (subset + ''' Optional. Default: ''', {})}


if __name__ == '__main__':

    import logging
    from macroeco.utils.workflow import Workflow
    from macroeco.empirical import Patch
    import macroeco.compare as comp
    from macroeco.output import OutputRarity

    wf = Workflow(required_params=required_params,
                optional_params=optional_params, clog=True, svers=__version__)
    
    for data_path, output_ID, params in wf.single_datasets():

        patch = Patch(data_path, subset=params['subset'])
        
        # NOTE: Only looks at rarity for SADs 
        patch_out = patch.sad(params['criteria'], clean=True)

        cmpr = comp.CompareSAD(patch_out, 
                 params['predicted_SAD_distributions'], patch=True)
        rarity = cmpr.compare_rarity(mins_list=params['rarity_measure'])
       
        OutputRarity(output_ID).output_rarity(rarity, data_path,
                                    cmpr.observed_data, criteria=cmpr.criteria)
        
        logging.info('Completed analysis %s\n' % output_ID)
    logging.info("Completed 'compare_rarity.py' script")




        








