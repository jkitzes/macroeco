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

explantion = ''' This script allows the user to compare observed energy
distributions against predicted energy distributions.  This script can compare
the individual community distribution (psi) and/or the species-level energy
distribution (theta). The choice of which energy metrics to compare can be
specified in the 'energy_metrics' parameter. The required parameters are below.

'subset' : How one would like to initially subset his data (see DataTable 
class docstring). 

'criteria' : How one would like to divide her dataset when caluculating the 
energy distributions. The theoretical distributions are compared against every
dataset that is generated from the cutting specified in this parameter.  See 
Patch class docstring for details

'dist_list_sed' : The list of sed distributions to which one could compare her
observed data.  Possible distributions to use are: 'theta'

'data_list_ied' : The list of ied distributions to which one could compare her
observed dara. Possible distributions to use are: 'psi'

'energy_metrics' : A list of either ['sed'], ['ied'], ['sed', 'ied'], or []
specifiying which energy metrics to examine.

This script outputs cdf and rad plots in which the observed energy distribution
is compared the distributions given in dist_list_sed and/or dist_list_ied.  For
each plot, a corresponding csv file with the same name as the plot except with
a .csv extension is output containing the data used to make the plot. For all
sed plots, the species name is printed on the top of the plot

'''

required_params = {'subset' : 'Dictionary of initial subsets', 'criteria' :
        'Dictionary of how to split the data', 'dist_list_sed' : 
    'List of sed distributions', 'dist_list_ied' : 'List of ied distributions',
    'energy_metrics' : 'List of energy metrics to examine'}

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




        








