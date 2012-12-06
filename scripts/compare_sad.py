#!/usr/bin/python

'''
Script to compare sads
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

gui_name = '''SAD'''#'''Analysis of Species Abundance Distributions'''

summary = '''Compares a dataset's observed species abundance distribution
against predicted species abundance distributions'''


# NOTE: Need to find a different way to specify which distributions they can
# use
predicted_SAD_distributions = '''A list of  SAD
distributions to which you can compare your observed data. 

You may use any number of the following SAD distributions: {!s} 

Example input: ['logser', 'plognorm_lt'] or ['nbd_lt']. The brackets MUST be
included.'''.format(gb.SAD_distributions)

rarity_measure = gb.rarity_measure + ''' In this analysis, the rarity
counts refer to individuals per species.'''

explanation = '''
ANALYSIS EXPLANATION\n
The analysis compare_sad allows you to compare an observed species abundance distribution
(SAD) against any number of predicted SADs. An SAD is a distribution of the
number of individuals within each species for an entire community.  This
distribution can also be thought of as the probability that a given species in
your community will have 'n' individuals. If you were to fully census a
community with S species and N individuals, the N individuals would be
distributed among the S species in a certain way.  Ecologists often predict
that this "certain way" will take the form of a logseries or lognormal
distribution.  These predicted distributions are among some of the
distributions against which you can compare your observed SAD. For more
information on SADs please see the provided references and the references
therein.

OUTPUT

The analysis compare_sad outputs three folders per dataset, a logfile.txt, and,
if possible, a map of the location(s) of the dataset(s). The three folders have
the names rank_abundance_plots_compare_sad_*, cdf_plots_compare_sad_*, and
summary_statistics_compare_sad_*. These folders contain cumulative density
function (cdf) plots and rank abundance distribution (rad) plots in which the
observed SAD distribution is compared to the distributions given in the
required parameter predicted_SAD_distributions.  For each plot, a corresponding
csv file with the same name as the plot except with a .csv extension is output
containing the data used to make the plot.  In the summary statistics folder,
summary .txt and .csv files are output containing summary values for the plot
and fitting statistics for the different predicted distributions you choose. 

If you choose to split up your data, individual SADs are generated for
each subsection.  For example, if you had data from three different years and
you choose to split your data by year, you would get three SAD comparisons, one
for each year.

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

Magurran, A. E. 1988. Ecological Diversity and Its Measurement. Princeton
University Press.

May, R. M. 1975. Patterns of species abundance and diversity. In Ecology and
Evolution of Communities (eds M. L. Cody and J. M. Diamond), Harvard University
Press.
'''.format(gb.subset, gb.criteria, rarity_measure,
predicted_SAD_distributions)

required_params = {'criteria' : gb.short_criteria + gb.req,
        'predicted_SAD_distributions': predicted_SAD_distributions + gb.req}

optional_params = {'subset' : (gb.short_subset + gb.optional, {}),
                       'rarity_measure' : (rarity_measure + gb.optional, [10])}

if __name__ == '__main__':

    import logging
    from macroeco.utils.workflow import Workflow
    from macroeco.empirical import Patch
    import macroeco.compare as comp
    from macroeco.output import SADOutput, make_directory
    import os

    wf = Workflow(required_params=required_params,
                 optional_params=optional_params, clog=True, svers=__version__)

    folder_name = 'SAD_analysis'
    make_directory(folder_name)
    cwd = os.getcwd()

    for data_path, output_ID, params in wf.single_datasets():
        
        os.chdir(os.path.join(cwd,folder_name))
        
        patch = Patch(data_path, params['subset'])
        sad = patch.sad(params['criteria'], clean=True)

        cmpr = comp.CompareSAD(sad,
                            params['predicted_SAD_distributions'], patch=True)
        summary = cmpr.summary(mins_list=params['rarity_measure'])
        
        sout = SADOutput(output_ID)
        sout.write_summary_table(summary, criteria=cmpr.criteria)
        sout.plot_rads(cmpr.compare_rads(), criteria=cmpr.criteria,
                                        species=cmpr.sad_spp_list)
        sout.plot_cdfs(cmpr.compare_cdfs(), cmpr.observed_data,
                        criteria=cmpr.criteria, species=cmpr.sad_spp_list)

        os.chdir(cwd)
        logging.info('Completed analysis %s\n' % output_ID)
    
    os.chdir(os.path.join(cwd,folder_name))    
    fout = open('README_compare_sad', 'w')
    with fout:
        fout.write(explanation)
    os.chdir(cwd)

    logging.info("Completed 'compare_sad.py' script")




        







