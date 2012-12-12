#!/usr/bin/python

'''
Script to get betas for given S and N
'''

__author__ = "Mark Wilber"
__copyright__ = "Copyright 2012, Regents of the University of California"
__credits__ = ["John Harte"]
__license__ = None
__version__ = "0.1"
__maintainer__ = "Mark Wilber"
__email__ = "mqw@berkeley.edu"
__status__ = "Development"


gui_name = '''Get METE Betas'''#'''Calculate Beta from Given N and S'''

summary = '''Calculate the Lagrange multiplier 'beta' for a given set of S and
N.'''

approximation = '''This parameter can be either 'True' or 'False'.  If it is
'False', the beta parameter will be calculated without an approximation and if
it is 'True', the beta parameter will be calculated with an approximation. When
N and S are large, both methods will return similar results.'''

explanation = '''
ANALYSIS EXPLANATION: \n This script takes in a csv file with columns 'S' and
'N' and calculates the beta value associated with each S and N pair. S
represents the total number of species in a given census and N represents the
total number of individuals.  The beta value is present in
the METE logseries and can be calculated with or without an approximation. The
beta value is a combination of lagrange multipliers that are derived from
maximizing entropy given some constraints. See references for more information
about beta and METE.

This script saves a csv file with columns 'N', 'S', and 'beta', where the
'beta' column contains the beta for the given N and S.  

PARAMETER EXPLANATION

*** approximation ***:

{0}

REFERENCES

Harte, J. 2011. Maximum Entropy and Ecology: A Theory of Abundance,
Distribution, and Energetics. Oxford University Press.

'''.format(approximation)

required_params = {'approximation' : approximation}

optional_params = {'subset' : ('''Not applicable for this analysis ''', {}),
                   'criteria' : ('Not applicable for this analysis', {})}

if __name__ == '__main__':

    import logging
    import numpy as np
    from macroeco.utils.workflow import Workflow
    from macroeco.distributions import logser_ut, logser_ut_appx
    from macroeco.utils.form_func import output_form
    from matplotlib.mlab import csv2rec

    wf = Workflow(required_params=required_params,
                optional_params=optional_params, clog=True, svers=__version__)

    for data_path, output_ID, params in wf.single_datasets():
        n_s = csv2rec(data_path)
        if set(n_s.dtype.names).intersection({'n', 's'}) != {'n', 's'}:
            raise KeyError('File %s does not have the required columns N ' %
                            data_path + 'and S')
        if params['approximation'] == False:
            lgsr = logser_ut(tot_obs=n_s['n'], n_samp=n_s['s'])
        else:
            lgsr = logser_ut_appx(tot_obs=n_s['n'], n_samp=n_s['s'])
        lgsr.pmf(1)
        betas = -np.log(np.array(lgsr.var['x']))
        n_s_new = np.array(zip(n_s['s'], n_s['n'], betas), dtype=[('S',
                               np.float), ('N', np.float), ('beta', np.float)])
        output_form(n_s_new, output_ID + '_betas')

        logging.info('Completed analysis %s\n' % output_ID)
    logging.info("Completed 'get_betas.py' script")
