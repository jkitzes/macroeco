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

import macroeco.utils.global_strings as global_str

gui_name = '''Analysis of Species Abundance Distributions'''

summary = '''Compares a dataset's observed species abundance distribution
against predicted species abundance distributions'''

explanation = '''This script allows you to compare an observed species
abundance distribution (SAD) against any number of predicted SADs. An SAD is a
distribution of the number of individuals within each species for an entire
community.  This distribution can also be thought of as the probability that a
given species in your community will have 'n' individuals. If you were to fully
census a community with S species and N individuals, the N individuals would be
distributed among the S species in a certain way.  Ecologists often predict
that this "certain way" will take the form of a logseries or lognormal
distribution.  These predicted distributions are among some of the
distributions against which you can compare your observed SAD. 

This script outputs cumulative density function (cdf) plots and rank abundance
distribution (rad) plots in which the observed SAD distribution is compared to
the distributions given in the required parameter predicted_SAD_distributions.
For each plot, a corresponding csv file with the same name as the plot except
with a .csv extension is output containing the data used to make the plot.  In
addition, a summary .txt file is output containing summary values for the plot
and fitting statistics for the different predicted distributions you
choose. If you choose to split up your data, individual SADs are generated for
each subsection.  For example, if you had data from three different years and
you choose to split your data by year, you would get three SAD comparisons, one
for each year.

For more information on SADs please see the following references and the
references therein:

Magurran, A. E. 1988. Ecological Diversity and Its Measuremnt. Princeton
University Press.

May, R. M. 1975. Patterns of species abundance and diversity. In Ecology and
Evolution of Communities (eds M. L. Cody and J. M. Diamond), Harvard University
Press.

'''

subset = global_str.subset

criteria = global_str.criteria

# NOTE: Need to find a different way to specify which distributions they can
# use
predicted_SAD_distributions = '''This parameter is the list of SAD
distributions to which you can compare your observed data. 

You may use any number of the following SAD distributions: %s 

Example input: ['logser', 'plognorm_lt'] or ['nbd_lt']. The brackets MUST be
included.''' % (global_str.SAD_distributions)

rarity_measure = global_str.rarity_measure + ''' In this analysis, the rarity
counts refer to individuals per species.'''

required_params = {'criteria' : criteria, 
                   'predicted_SAD_distributions' : predicted_SAD_distributions,
                   'rarity_measure' : rarity_measure}

optional_params = {'subset' : subset}

if __name__ == '__main__':

    import logging
    from macroeco.utils.workflow import Workflow
    from macroeco.empirical import Patch
    import macroeco.compare as comp
    from macroeco.output import DistOutput

    wf = Workflow(required_params=required_params,
                  clog=True, svers=__version__)
    
    for data_path, output_ID, params in wf.single_datasets():
        try:
            subset = params['subset']
        except:
            logging.info("Not subsetting anything (default value)")
            subset = {}
        
        patch = Patch(data_path, subset)
        sad = patch.sad(params['criteria'], clean=True)

        cmpr = comp.CompareDistribution(sad,
                            params['predicted_SAD_distributions'], patch='sad')
        rads = cmpr.compare_rads()
        summary = cmpr.summary(rads, mins_list=params['rarity_measure'])
        
        sout = DistOutput(output_ID, 'sad')
        sout.write_summary_table(summary, criteria=cmpr.criteria)
        sout.plot_rads(rads, criteria=cmpr.criteria)
        sout.plot_cdfs(cmpr.compare_cdfs(), cmpr.data_list,
                        criteria=cmpr.criteria)
        logging.info('Completed analysis %s\n' % output_ID)
    logging.info("Completed 'compare_sad.py' script")




        







