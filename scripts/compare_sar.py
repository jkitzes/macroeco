#!/usr/bin/python

'''
Script to compare sars
'''

__author__ = "Mark Wilber"
__copyright__ = "Copyright 2012, Regents of the University of California"
__credits__ = ["John Harte"]
__license__ = None
__version__ = "0.1"
__maintainer__ = "Mark Wilber"
__email__ = "mqw@berkeley.edu"
__status__ = "Development"

gui_name = ''' Analyze Species-Area Relationships'''

summary = '''Compares a dataset's observed species-area relationships against 
theoretical species-area relationships'''

class global_str:
    subset = '''You should examine the columns in your data set and decide if
    you would like to subset your data in some particular way before the
    analysis begins. It is important to note that only the subsetted data will
    be analyzed.  For example,  if you have a column named 'year' in your data
    set with values 1998, 1999, and 2000 and you only want to look at the year
    2000 for a particular analysis, you should select the == operator from the
    drop down list and type 2000 in the value field.  Similarly, you could use
    <, >, <=, >=, or !='''

    criteria = '''You should examine the columns in your dataset and decide if
    you would like to divide the data in a particular way for this analysis.
    For example, if the you have a spatial dataset with x,y coordinates and you
    are interested in examining macroecological metrics for two separate halves
    of your plot along the x coordinate, you could cut the x coordinate in two
    halves by giving the 'x' column a value of 2.  If the column that you would
    like to divide contains discrete values (e.g. year), you could enter the
    keyword 'split' and each unique value will be analyzed separately.
    Conversely, the value 'whole' could be given to specify the entire column.
    The value 'whole' is equivalent to 1 or leaving the value blank. If you
    would like to divide a given column, please select the word 'division' from
    the GUI dropdown menu and input a value as discussed above.\n\n

    There are four special words that can be used on a given column: 'species',
    'energy', 'count', and 'mass'.  When assigned to a column in your data set,
    the special word 'species' indicates the column that contains your species
    IDs, the special word 'energy' indicates the column that contains some type
    of energy measure, the special word 'mass' indicates a column that contains
    some type of mass measure, and the special word 'count' indicates the
    column that contains your species counts.  In the GUI, these special words
    can be chosen from the dropdown menu next to each column header. The
    special word 'species' MUST be assigned for every analysis.  If the special
    word 'count' is not assigned, the species counts are all assumed to be
    one.\n\n'''

    rarity_measure = '''This parameter allows you to specify the counts that
    you will consider rare.  If, for example, you want to know how many species
    in your plot have an abundance of 2 or less you would set this parameter to
    2. If you enter more then one value, each value will be examined. Example
    input: [2] or [2, 5]. The brackets MUST be included.'''

    SAD_distributions = ''' 'logser','logser_ut', 'logser_ut_appx',
    'plognorm_lt', 'nbd_lt', 'geo_ser', 'broken_stick', 'lognorm' '''

    SSAD_distributions = ''' 'nbd', 'binm', 'tgeo', 'fgeo', 'fnbd', 'pois' '''

subset = '''Specifications for how you want to subset your data before the
analysis.  Only the subsetted data will be included in the analysis.  See
explanation link for more detail.'''

criteria = '''Specifications for how you want to divide your data during the
analysis.  See explanation link for more detail.'''

columns_to_divide = '''This parameter specifies which spatial columns you would
like to divide for the SAR analysis.  For example, if your data had spatial
columns 'x' and 'y' you would enter: ('x', 'y')'''

list_of_divisions_on_columns = '''This parameter specifies how you would like
to divide the columns you named in the columns_to_divide parameter.  For
example, if you specified that you wanted to divide on columns 'x' and 'y' in
the parameter columns_to_divide and you wanted to divide your plot into fourths
and eighths you would input: [(2,2), (2,4)].  The first value splits the plot
into fourths and the second splits the plot into eighths.  The values within
each parentheses are divisions on 'x' and 'y', respectively.  '''

predicted_SAR_curves = '''A list of SAR curves to which you can compare your
observed data.

You may use any number of the following SAR distributions: 'powerlaw',
'mete_sar', 'mete_sar_iter', 'logser_ut_binm', 'lognorm_binm', etc.

Example input: ['mete_sar', 'powerlaw', 'mete_sar_iter']
'''

name = '''A name for the plot that will be generated.

Example input: My SAR Plot
'''



explanation = '''ANALYSIS EXPALANTION\n
This script allows you to compare observed species-area relationships (SAR)
with any number of predicted SARs.  An SAR is a commonly used macroecological
metric which examines the number of species found in a given area. All SAR
curves show increasing species with increasing area, but shape of this
increasing curve differs depending on the theory used to derive it. It is
common practice to examine SAR plots on a log-log scale because the curve often
becomes close to a straight line. Using this script, you can generate the
observed SARs for any nested plots within your data.  For example, if you had a
fully censused plot with spatial coordinates x,y and you wanted to examine an
SAR looking at the anchor area (A), 1/2 * A, 1/4 * A, and 1/8 * A, you
would input the appropriate parameters and this analysis would divide the plot
into halves, fourths, and eighths and take the average number of species across
all of the smaller plots. Therefore, you can have a fractional average number
of species per areas. For additional information on SARs, please see the
provided references and the references therein.

This script ouputs log(species) vs. log(area fraction) plots.  Area fraction
means that the anchor area (largest area) for the plot is 1 and all smaller
subplots are fractions less than one. Each plot will have the observed SAR and
any SAR generated by a predicted SAR specified in the predicted_SAR_curves
parameter. In addition to this plot, this script will output csv files with the
same file name as the plot, but with the name of the predicted SAR appended to
the end.  These files contain the data for the given SAR used to make the plot.
With this data, you can re-plot the data in anyway you chose.

PARAMETER EXPLANATIONS

*** subset ***:

{0}

*** criteria ***:

{1}

For this (compare_sar) analysis, if the columns you filled in the parameter
columns_to_divide are repeated in the parameter criteria, the values that you
assigned in criteria will be ignored. Generating multiple SARs is not currently
implemented. 

*** columns_to_divide **

{2}

*** list_of_divisions_on_columns ***

{3}

*** predicted_SAR_curves ***

{4}

*** name ***

{5}

REFERENCES

Harte, J. 2011. Maximum Entropy and Ecology: A Theory of Abundance,
Distribution, and Energetics. Oxford University Press.

Rosenzweig, M. L. 1995. Species Diversity in Space and Time. Cambridge
University Press.

'''.format(global_str.subset, global_str.criteria, columns_to_divide,
list_of_divisions_on_columns, predicted_SAR_curves, name)


required_params = {'criteria' : criteria,
                   'columns_to_divide' : columns_to_divide,
                   'list_of_divisions_on_columns' :
                   list_of_divisions_on_columns,
                   'predicted_SAR_curves' : predicted_SAR_curves,
                   'name' : name}

optional_params = {'subset' : (subset + ''' Optional. Default: ''', {})}

if __name__ == '__main__':

    import logging
    from macroeco.utils.workflow import Workflow
    from macroeco.empirical import Patch
    import macroeco.compare as comp
    from macroeco.output import SAROutput
    from copy import deepcopy

    wf = Workflow(required_params=required_params, clog=True, 
                                                            svers=__version__)
    
    for data_path, output_ID, params in wf.single_datasets():
        for optpar in optional_params: #TODO: Move to workflow
            if not optpar in params:
                logging.info("Default value for {!s}: {!s}".format(optpar,
                              str(optional_params[optpar][1])))
                params[optpar] = optional_params[optpar][1]

        try:
            params['list_of_divisions_on_columns'].index((1,1))
        except:
            logging.info("Adding base area to parameter " +
                                       "'list_of_divisions_on_columns': (1,1)")
            params['list_of_divisions_on_columns'].append((1,1))
        sad_criteria = deepcopy(params['criteria'])

        for nm in params['columns_to_divide']:
            sad_criteria[nm] = 'whole'
        patch = Patch(data_path, subset=params['subset'])
        sad = patch.sad(sad_criteria, clean=True)
        sar = patch.sar(params['columns_to_divide'], 
                    params['list_of_divisions_on_columns'], params['criteria'])
        cmpr = comp.CompareSARCurve([sar], params['predicted_SAR_curves'],
                                                    [sad[0][1]], patch=True)
        srout = SAROutput(output_ID)
        srout.plot_sars(cmpr.compare_curves(), names=[params['name']])
        logging.info('Completed analysis %s\n' % output_ID)
    logging.info("Completed 'compare_sar' script")




        








