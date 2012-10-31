#!/usr/bin/python

''' Script to compare energy distributions '''

__author__ = "Mark Wilber"
__copyright__ = "Copyright 2012, Regents of the University of California"
__credits__ = ["John Harte"]
__license__ = None
__version__ = "0.1"
__maintainer__ = "Mark Wilber"
__email__ = "mqw@berkeley.edu"
__status__ = "Development"


gui_name = '''Analysis of Macroecological Energy Metrics'''

summary = '''Compares a dataset's observed energy metrics against predicted
energy metrics'''

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
    For example, if the you have a spatial dataset with x,y coordinates and you
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
    assumed to be one.
    
    If there are columns in your data that are not relevant for this analysis
    leave the value in the dropdown box as 'NA'.  Columns designated 'NA' will
    not influence the analysis.\n\n'''

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
For this compare_energy anlysis, you must designate your species column with
the special word 'species' and your mass or energy column with the special word
'mass' or 'energy'.  These special words can be found in the dropdown menu. You
are not required to fill any additional columns for this analysis. Please see
analysis explanation for more detail and examples.'''

predicted_SED_distributions = '''This parameter is the list of SED
distributions to which you can compare your observed data. 

You may use any number of the following SED distributions: 'theta'

Example input: ['theta'].  The brackets MUST be included.  '''

predicted_IED_distributions = '''This parameter is the list of IED
distributions to which you can compare your observed data. 

You may use any number of the following IED distributions: 'psi'

Example input: ['psi']. The brackets MUST be included.  '''

energy_metrics = ''' This parameter allows you to specify which energy
metric(s) you would like to examine. Chose one of the following and copy
exactly: ['SED'], ['IED'], or ['SED', 'IED'].'''


explanation = ''' 
ANALYSIS EXPLANATION\n
This script allows you to compare observed energy metrics
against predicted energy metrics. There are two energy metrics that this script
can compare: individual energy distributions (IED) and species-level energy
distributions (SED).  The IED is a distribution of the energy of each
individual within a community.  The energy of each individual can be calculated
from the the biomass using the 3/4 allometric scaling law.  Other proxies for
energy, such as leaf surface area, can be used as well.  The IED is species
blind; it does not consider what species an individual belongs to. For
normalization purposes, each value of the empirical IED is divided by the
smallest empirical energy value. 

The SED is a distribution of the energy of each individual within a species.
So, for example, a community with 30 speices would have 30 SED's, one for each
species.  The entire community's energy distribution is normalized in the exact
same way as the IED and then individuals are placed into their respective
species and the SED's are generated. For more information on energy distributions please see the reference and
references therein.

OUTPUT

This script outputs rank energy distribution (red) plots in which the observed
energy distribution is compared to the distributions given in the required
parameter predicted_SED_distributions and/or predicted_IED_distributions. These
files are output as .png files. For each plot, a corresponding csv file with
the same name as the plot except with a .csv extension is output containing the
data used to make the plot. For all SED plots, the species name is printed on
the top of the plot.  For all of the IED plots, the criteria used to make plot
is printed on the top of the plot. 

PARAMETER EXPLANATIONS

*** subset ***:

{0}

*** criteria ***:

{1}

A column must be designated as 'energy' and/or 'mass' for this analysis. If
not, your analysis will not run.

*** predicted_SED_distributions ***:

{2}

*** predicted_IED_distributions ***:

{3}

*** energy_metrics ***:

{4}

REFERENCES

Harte, J. 2011. Maximum Entropy and Ecology: A Theory of Abundance,
Distribution, and Energetics. Oxford University Press.

'''.format(global_str.subset, global_str.criteria, predicted_SED_distributions,
predicted_IED_distributions, energy_metrics)


required_params = {'criteria' : criteria,
                   'predicted_SED_distributions' : predicted_SED_distributions,
                   'predicted_IED_distributions' : predicted_IED_distributions,
                   'energy_metrics' : energy_metrics}

optional_params = {'subset' : (subset + ''' Subsetting is optional. Default: 
                    ''', {})}

if __name__ == '__main__':

    import logging
    from macroeco.utils.workflow import Workflow
    from macroeco.empirical import Patch
    import macroeco.compare as comp
    from macroeco.output import SEDOutput, IEDOutput

    wf = Workflow(required_params=required_params, clog=True, 
                                                            svers=__version__)
    
    for data_path, output_ID, params in wf.single_datasets():
        for optpar in optional_params: #TODO: Move to workflow
            if not optpar in params:
                logging.info("Default value for {!s} in {!s}: {!s}".format(optpar,
                              output_ID,  str(optional_params[optpar][1])))
                params[optpar] = optional_params[optpar][1]
        # Put data in patch object
        patch = Patch(data_path, subset=params['subset'])

        # Calculate empirical metrics
        sad = patch.sad(params['criteria'])
        if set(params['energy_metrics']) == set(['SED', 'IED']) or\
           set(params['energy_metrics']) == set(['sed', 'ied']) :
            cmengy = patch.ied(params['criteria'])
            spengy = patch.sed(params['criteria'])

            # Make comparison objects 
            cmprt = comp.CompareSED((spengy, cmengy, sad),
                             params['predicted_SED_distributions'], patch=True)
            cmprp = comp.CompareIED((cmengy, sad), 
                             params['predicted_IED_distributions'], patch=True)
            
            # Make output objects and output plots
            sout = SEDOutput(output_ID)
            soup = IEDOutput(output_ID)
            sout.plot_reds(cmprt.compare_reds(), criteria=cmprt.sad_criteria)
            soup.plot_reds(cmprp.compare_reds(), criteria=cmprp.sad_criteria)

        elif set(params['energy_metrics']) == set(['IED']) or\
             set(params['energy_metrics']) == set(['ied']):
            cmengy = patch.ied(params['criteria'])
            cmprp = comp.CompareIED((cmengy, sad), 
                              params['predicted_IED_distributions'],patch=True)
            soup = IEDOutput(output_ID)
            soup.plot_reds(cmprp.compare_reds(), criteria=cmprp.sad_criteria)

        elif set(params['energy_metrics']) == set(['SED']) or\
             set(params['energy_metrics']) == set(['sed']):
            cmengy = patch.ied(params['criteria'])
            spengy = patch.sed(params['criteria'])

            cmprt = comp.CompareSED((spengy, cmengy, sad),
                             params['predicted_SED_distributions'], patch=True)
            sout = SEDOutput(output_ID)
            sout.plot_reds(cmprt.compare_reds(), criteria=cmprt.sad_criteria)

        logging.info('Completed analysis %s\n' % output_ID)
    logging.info("Completed 'compare_energy.py' script")




        








