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

import macroeco.utils.global_strings as gb 

gui_name = '''Energy Metrics'''#'''Analysis of Macroecological Energy Metrics'''

summary = '''Compares a dataset's observed energy metrics against predicted
energy metrics'''


predicted_SED_distributions = '''This parameter is the list of SED
distributions to which you can compare your observed data. 

You may use any number of the following SED distributions: 'theta'

Example input: ['theta'].  The brackets MUST be included.  '''

predicted_IED_distributions = '''This parameter is the list of IED
distributions to which you can compare your observed data. 

You may use any number of the following IED distributions: 'psi'

Example input: ['psi']. The brackets MUST be included.  '''

predicted_ASED_distributions = '''This parameter is the list of ASED
distributions to which you can compare your observed data. 

You may use any number of the following ASED distributions: 'nu'

Example input: ['nu']. The brackets MUST be included.  '''

energy_metrics = ''' This parameter allows you to specify which energy
metric(s) you would like to examine. Example: ['SED'], ['IED'], ['SED', 'IED']
['ASED'], ['ied', 'ASED', 'sed'].'''


explanation = ''' 
ANALYSIS EXPLANATION\n
The compare_energy analysis allows you to compare observed energy metrics
against predicted energy metrics. There are three energy metrics that this
script can compare: individual energy distributions (IED), species-level energy
distributions (SED), and average species energy distribution (ASED).  The IED
is a distribution of the energy of each individual within a community.  The
energy of each individual can be calculated from the biomass using the 3/4
allometric scaling law.  Other proxies for energy, such as leaf surface area,
can be used as well.  The IED is species blind; it does not consider what
species an individual belongs to. For normalization purposes, each value of the
empirical IED is divided by the smallest empirical energy value. An example of
and IED distribution is the 'psi' distribution given in Harte (2011).

The SED is a distribution of the energy of each individual within a species.
For example, a community with 30 speices would have 30 SED's, one for each
species.  The entire community's energy distribution is normalized in the exact
same way as the IED and then individuals are placed into their respective
species and the SED's are generated. For more information on energy
distributions please see the reference and references therein.  An example of
an SED is the 'theta' distribution given in Harte (2011).

The ASED is a distribution of the average energy of each species within a
community. For example, a community with 30 species would have one ASED with
thirty data points; one for each species.  The average energy of a species is
the average energy of all individuals within that species. Similar to the SED
mentioned above, the IED is normalized before the ASED is calculated.  An
example of an ASED is the 'nu' distribution given in Harte (2011).

OUTPUT

This analysis outputs up to three folders per dataset, a logfile.txt, and, if
possible, a .png file with a map of the location of the datasets(s). The three
possible folders begin with the names ied_rank_energy_plots_compare_energy_*,
sed_rank_energy_plots_compare_energy_*, and
ased_rank_energy_plots_compare_energy_* .  Within each folder, there are rank
energy distribution (red) plots in which the observed energy distribution is
compared to the distributions given in the required parameter
predicted_SED_distributions, predicted_IED_distributions, and/or
predicted_ASED_distributions. These files are output as .png files. For each
plot, a corresponding csv file with the same name as the plot except with a
.csv extension is output containing the data used to make the plot. For all SED
plots, the species name specified in the plot title.  For all of the IED and
ASED plots, the criteria used to make plot is printed on the right hand side of
the plot. 

The logfile.txt contains the analysis process information. Please see the
logfile if the analysis fails.

PARAMETER EXPLANATIONS

*** subset ***:

{0}

*** criteria ***:

{1}
A column must be designated as 'energy' and/or 'mass' for this analysis.  If a
column is designated mass, all values will be raised to the 3/4 power in
accordance with metabolic scaling theory and then normalized as described in 
the analysis explanation.

*** predicted_SED_distributions ***:

{2}

*** predicted_IED_distributions ***:

{3}

*** predicted_ASED_distributions ***:

{5}

*** energy_metrics ***:

{4}

REFERENCES

Harte, J. 2011. Maximum Entropy and Ecology: A Theory of Abundance,
Distribution, and Energetics. Oxford University Press.

'''.format(gb.subset, gb.criteria, predicted_SED_distributions,
predicted_IED_distributions, energy_metrics, predicted_ASED_distributions)


required_params = {'criteria' : gb.req + gb.short_criteria,
        'predicted_SED_distributions' :gb.req + predicted_SED_distributions,
        'predicted_IED_distributions' : gb.req + predicted_IED_distributions,
        'energy_metrics' : gb.req + energy_metrics}

optional_params = {'subset' : (gb.short_subset + gb.optional, {})}

if __name__ == '__main__':

    import logging
    from macroeco.utils.workflow import Workflow
    from macroeco.empirical import Patch
    import macroeco.compare as comp
    from macroeco.output import SEDOutput, IEDOutput, ASEDOutput, make_directory
    import os

    wf = Workflow(required_params=required_params,
                optional_params=optional_params, clog=True, svers=__version__)

    folder_name = 'Energy_analysis'
    make_directory(folder_name)
    cwd = os.getcwd()
    
    for data_path, output_ID, params in wf.single_datasets():

        os.chdir(os.path.join(cwd,folder_name))

        # Put data in patch object
        patch = Patch(data_path, subset=params['subset'])

        # Calculate empirical metrics
        sad = patch.sad(params['criteria'], clean=True)
        ied = patch.ied(params['criteria'])

        # Check which distributions where specified
        dist_str = str(params['energy_metrics'])
        sed_there, ied_there, ased_there = (False, False, False)
        if dist_str.find('sed') != -1 or dist_str.find('SED') != -1:
            sed_there = True
        if dist_str.find('ied') != -1 or dist_str.find('IED') != -1:
            ied_there = True
        if dist_str.find('ased') != -1 or dist_str.find('ASED') != -1:
            ased_there = True

        if sed_there:
            sed = patch.sed(params['criteria'])
            cmprt = comp.CompareSED((sed, ied, sad),
                             params['predicted_SED_distributions'], patch=True)
            sout = SEDOutput(output_ID)
            sout.plot_reds(cmprt.compare_rads(), criteria=cmprt.criteria)

        if ied_there:
             cmprp = comp.CompareIED((ied, sad), 
                             params['predicted_IED_distributions'], patch=True)
             soup = IEDOutput(output_ID)
             soup.plot_reds(cmprp.compare_rads(), criteria=cmprp.criteria)

        if ased_there:
            ased = patch.ased(params['criteria'])
            cmpra = comp.CompareASED((ased, ied, sad), ['nu'], patch=True)
            soua = ASEDOutput(output_ID)
            soua.plot_reds(cmpra.compare_rads(), criteria=cmpra.criteria,
                                                   species=cmpra.sad_spp_lists)


        os.chdir(cwd)

        logging.info('Completed analysis %s\n' % output_ID)

    os.chdir(os.path.join(cwd,folder_name))
    fout = open('README_compare_energy', 'w')
    with fout:
        fout.write(explanation)
    os.chdir(cwd)

    logging.info("Completed 'compare_energy.py' script")

