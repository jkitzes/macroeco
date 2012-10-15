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

gui_name = '''Get Betas'''

summary = '''Calculate the Lagrange multiplier 'beta' for a given set of S and N.'''

explanation = '''This script takes in a csv file with columns 'S' and 'N' and
calculates the beta parameter associated with each 'S' and 'N' pair.  The beta
parameter is present in the METE logseries (see logser_ut and logser_ut_appx
classes) and can be calculated with or without an approximation. 
The parameter needed for this analysis is:

'approximation' : This parameter can be either 'True' or 'False'.  If it is
'False', the beta parameter will be calculated without an approximation and if
it is 'True', the beta parameter will be calculated with an approximation. 

This script saves a csv file with columns 'N', 'S', and 'beta', where the
'beta' column contains the 'beta' for the given 'N' and 'S'.  

'''

required_params = {'approximation' : 'A boolean that determines if beta is ' +\
                    ' approximated'}

if __name__ == '__main__':

    import logging
    import numpy as np
    from macroeco.utils.workflow import Workflow
    from macroeco.distributions import logser_ut, logser_ut_appx
    from macroeco.utils.form_func import output_form
    from matplotlib.mlab import csv2rec

    wf = Workflow(required_params=required_params, clog=True, 
                                                            svers=__version__)

    for data_path, output_ID, params in wf.single_datasets():
        n_s = csv2rec(data_path)
        if set(n_s.dtype.names).intersection({'n', 's'}) != {'n', 's'}:
            raise KeyError('File %s does not have the required columns N ' %
                            data_path + 'and S')
        if params['approximation'] == False:
            lgsr = logser_ut(tot_obs=n_s['n'], n_samp=n_s['s'])
        else:
            lgsr = logser_ut_appx(tot_obs=n_s['n'], n_samp=n_s['s'])
        betas = -np.log(np.array(lgsr.pmf(1)[1]['x']))
        n_s_new = np.array(zip(n_s['s'], n_s['n'], betas), dtype=[('S',
                               np.float), ('N', np.float), ('beta', np.float)])
        output_form(n_s_new, output_ID + '_betas')

        logging.info('Completed analysis %s\n' % output_ID)
    logging.info("Completed 'get_betas.py' script")
