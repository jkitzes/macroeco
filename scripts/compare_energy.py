#!/usr/bin/python

'''
Script to compare energy distributions
'''

__author__ = "Mark Wilber"
__copyright__ = "Copyright 2012, Regents of the University of California"
__credits__ = ["John Harte"]
__license__ = None
__version__ = "0.1"
__maintainer__ = "Mark Wilber"
__email__ = "mqw@berkeley.edu"
__status__ = "Development"

gui_name = '''Energetics Analysis'''

summary = '''Compares a dataset's observed energy distributions against 
             theoretical energy distributions'''

explantion = ''' This script allows the user to compare observed energy metrics
against predicted energy metrics. There are two energy metrics that this script
can compare: individual energy distributions (ied) and species-level energy
distributions (sed).  The ied is a distribution of the energy of each
individual within a community.  The energy of each individual can be calculated
from the the biomass using the 3/4 allometric scaling law.  Other proxies for
energy, such as leaf surface area, can be used as well.  The ied is species
blind; it does not consider what species an individual belongs to. For
normalization purposes, each value of the empirical ied is divided by the
smallest empirical energy value. 

The sed is a distribution of the energy of each individual within a species.
So, for example, a community with 30 speices would have 30 sed's, one for each
species.  The entire communities energy distribution is normalized in the exact
same way as the ied then individuals are placed into their respective species
and the sed's are generated. 

This script outputs cumulative density function (cdf) plots and rank abundance
distribution (rad) plots in which the observed energy distribution is compared
to the distributions given in the required parameter
predicted_distributions_sed and/or predicted_distributions_ied.  For each plot,
a corresponding csv file with the same name as the plot except with a .csv
extension is output containing the data used to make the plot. For all sed
plots, the species name is printed on the top of the plot
'''

subset = '''You should examine the columns in your data set and decide if you
would like to subset his data in some particular way BEFORE the analysis
begins. It is important to note that only the subsetted data will be analyzed.
For example,  if the user has column names 'year' in his data set with values
1998, 1999, and 2000 and he only wants to look at the year 2000 for a
particular analysis, he should enter ==2000 next as the value in the 'year'
field.  Similarly <, >, <=, >=, !=, etc. operators could be used.'''

criteria = '''You should examine the columns in your dataset and decide if you
would like to divide the data in a particular way FOR THE analysis. For
example, if the you have a spatial dataset with x,y coordinates and you are
interested in examining the ied for two separate halves of your plot along the
x coordinate she could cut the x coordinate in two halves by giving the 'x'
column a value of 2.  If the column that she would like to divide contains
discrete values (e.g. year), she could enter the value 'split' and each unique
discrete value will be analyzed separately. Conversely, the the value 'whole'
could be given could be given to specify the entire column.  The value 'whole'
is equivalent to 1 or leaving the value blank. There are four special words
that can be used on a given column: 'species', 'energy', 'count', and 'mass'.
The value 'species' MUST be given to the column that contains species ids. The
other three special words do not necessarily have to be given.'''


predicted_distributions_sed = '''The list of sed distributions to which you can
compare your observed data. These MUST be entered in the following format:
['example_1', 'example_2'] or ['example_1']. Possible distributions are:
'theta' '''

predicted_distributions_ied = '''The list of ied distributions to which you
could compare your observed data. These MUST be entered in the following
format: ['example_1', 'example_2'] or ['example_1']. P Possible distributions
to use are: 'psi' '''

energy_metrics = '''Type exactly ['sed'], ['ied'], ['sed', 'ied']. This
specifies which energy metric(s) you would like to examine'''

required_params = {'subset' : subset, 'criteria' : criteria,
                   'predicted_distributions_sed' : predicted_distributions_sed,
                   'predicted_distributions_ied' : predicted_distributions_ied,
                   'energy_metrics' : energy_metrics}

if __name__ == '__main__':

    import logging
    from macroeco.utils.workflow import Workflow
    from macroeco.empirical import Patch
    import macroeco.compare as comp
    from macroeco.output import SEDOutput, IEDOutput

    wf = Workflow(required_params=required_params, clog=True, 
                                                            svers=__version__)
    
    for data_path, output_ID, params in wf.single_datasets():

        # Put data in patch object
        patch = Patch(data_path, subset=params['subset'])

        # Calculate empirical metrics
        sad = patch.sad(params['criteria'])
        if set(params['energy_metrics']) == set(['sed', 'ied']):
            cmengy = patch.ied(params['criteria'])
            spengy = patch.sed(params['criteria'])

            # Make comparison objects 
            cmprt = comp.CompareSED((spengy, cmengy, sad),
                                    params['dist_list_sed'], patch=True)
            cmprp = comp.CompareIED((cmengy, sad), params['dist_list_ied'],
                                        patch=True)
            
            # Make output objects and output plots
            sout = SEDOutput(output_ID)
            soup = IEDOutput(output_ID)
            sout.plot_reds(cmprt.compare_reds(), criteria=cmprt.sad_criteria)
            soup.plot_reds(cmprp.compare_reds(), criteria=cmprp.sad_criteria)

        elif set(params['energy_metrics']) == set(['ied']):
            cmengy = patch.ied(params['criteria'])
            cmprp = comp.CompareIED((cmengy, sad), params['dist_list_ied'],
                                        patch=True)
            soup = IEDOutput(output_ID)
            soup.plot_reds(cmprp.compare_reds(), criteria=cmprp.sad_criteria)

        elif set(params['energy_metrics']) == set(['sed']):
            cmengy = patch.ied(params['criteria'])
            spengy = patch.sed(params['criteria'])

            cmprt = comp.CompareSED((spengy, cmengy, sad),
                                    params['dist_list_sed'], patch=True)
            sout = SEDOutput(output_ID)
            sout.plot_reds(cmprt.compare_reds(), criteria=cmprt.sad_criteria)

        logging.info('Completed analysis %s\n' % output_ID)
    logging.info("Completed 'compare_energy.py' script")




        








