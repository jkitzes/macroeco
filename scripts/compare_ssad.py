#!/usr/bin/python

'''
Script to compare ssads
'''

__author__ = "Mark Wilber"
__copyright__ = "Copyright 2012, Regents of the University of California"
__credits__ = ["John Harte"]
__license__ = None
__version__ = "0.1"
__maintainer__ = "Mark Wilber"
__email__ = "mqw@berkeley.edu"
__status__ = "Development"

from macroeco.utils import global_strings

gui_name = '''Analysis of Species-level Spatial Abundance Distributions'''

summary = '''Compares a dataset's observed species-level spatial abundance
distribution against predicted species-level spatial abundance distributions'''

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

# NOTE: Need to find a different way to specify which distributions they can
# use
predicted_SSAD_distributions = '''A list of SSAD
distributions to which you can compare your observed data. 

You may use any number of the following SSAD distributions: {!s} 

Example input: ['binm', 'pois'] or ['fnbd']. The brackets MUST be
included.'''.format(global_strings.SSAD_distributions)

rarity_measure = global_strings.rarity_measure + ''' In this analysis, the rarity
counts refer to the number of individuals of a single species per cell.'''

explanation = '''
ANALYSIS EXPLANATION\n This script allows you to compare an observed
species-level spatial abundance distribution (SSAD) against any number of
predicted SSADs. An SSAD can be thought of as the probability that a given
species with n_o individuals will have n individuals in a cell of size A <=
A_o, where A_o is the anchor scale for a given plot.  Some common predictions
of SSADs include binomial distributions and negative binomial distributions.
These distributions are are among some of the distributions against which you
can compare your observed SSADs. Because SSADs are species specific, a given
plot with 30 species will have 30 associated SSADs regardless of how you divide
it.  For example, if you consider the whole plot, each species has an SSAD with
one sample that is equal to its abundance.  If you were to divide the plot into
sixteenths, each species would have an SSAD with 16 samples which summed to its
total abundance in the plot. For more information on SSADs please see the
provided references and the references therein.

OUTPUT

This script outputs cumulative density function (cdf) plots and rank abundance
distribution (rad) plots in which the observed SSAD distribution for each
species is compared to the distributions given in the required parameter
predicted_SSAD_distributions (png files).  For each plot, a corresponding csv
file with the same name as the plot except with a .csv extension is output
containing the data used to make the plot.  In addition, a summary .txt file is
output containing summary values for the plot and fitting statistics for the
different predicted distributions you choose. If you perform an SSAD analysis
on a plot with 30 species, you will get 30 cdf plots, 30 rad plots, 30 + 30 csv
files, and 30 summary txt files, one for each species.   

PARAMETER EXPLANATIONS

*** subset ***:

{0}

*** criteria ***:

{1}

*** rarity_measure ***:

{2}

*** predicted_SSAD_distributions ***:

{3}

REFERENCES

Harte, J. 2011. Maximum Entropy and Ecology: A Theory of Abundance,
Distribution, and Energetics. Oxford University Press.
'''.format(global_strings.subset, global_strings.criteria, rarity_measure,
predicted_SSAD_distributions)

required_params = {'criteria' : criteria, 'rarity_measure' : rarity_measure,
		           'predicted_SSAD_distributions': predicted_SSAD_distributions}

optional_params = {'subset' : (subset + ''' Optional. Default: ''', {})}

if __name__ == '__main__':

    import logging
    from macroeco.utils.workflow import Workflow
    from macroeco.empirical import Patch
    import macroeco.compare as comp
    from macroeco.output import DistOutput

    wf = Workflow(required_params=required_params, clog=True, 
                                                            svers=__version__)
    
    for data_path, output_ID, params in wf.single_datasets():
        for optpar in optional_params: #TODO: Move to workflow
            if not optpar in params:
                logging.info("Default value for {!s} in {!s}: {!s}".format(optpar,
                              output_ID,  str(optional_params[optpar][1])))
                params[optpar] = optional_params[optpar][1]

        patch = Patch(data_path, subset=params['subset'])
        ssad = patch.ssad(params['criteria'])

        cmpr = comp.CompareSSAD(ssad,
                           params['predicted_SSAD_distributions'], patch=True)
        rads = cmpr.compare_rads()

        sout = DistOutput(output_ID, 'ssad')
        sout.write_summary_table(cmpr.summary(rads,
                   mins_list=params['rarity_measure']), criteria=cmpr.spp_list)
        sout.plot_rads(rads, criteria=cmpr.spp_list)
        sout.plot_cdfs(cmpr.compare_cdfs(), cmpr.data_list,
                        criteria=cmpr.spp_list)
        logging.info('Completed analysis %s\n' % output_ID)
    logging.info("Completed 'compare_ssad.py' script")




        








