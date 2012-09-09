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

#TODO: Need to fill out docstrings and all that good stuff
gui_name = '''SAR Analysis'''

summary = '''Compares a dataset's observed sar against theoretical sars'''

explantion = '''This script takes in a dataset(s) and list of curves to which
compare the observed datasets sad will be compared.  The required parameters
for the script are the following:

'subset' : How one would like to initially subset his data (see DataTable 
class docstring). 

'div_cols' : A tuple of specifying the columns that will be divided during the
sar analysis. e.g. ('x', 'y')

'div_list' : A list of tuples, in the same order as the columns in div_list, that
specify the divisions one would like to make in order to generate the sar.  A
division of (1,1) indicates the whole plot.

'sar_criteria' : This parameter is a dictionary that specifies which column is
the data set is the column is species counts and which column is the column
with species names. e.g. {'spp' : 'species', 'spp_count' : 'count'.  In this case,
column name 'spp' is the species column and column 'spp_count' is the count
column. See Patch method sar() for more info.

'curve_list' : A list with the name of the SARCurve objects that will be
compared to the empirical SAR.  If none are given, the empirical SAR is just
plotted on its own.

'names' : List with the desired name of the plot.

For each dataset, this script generates a log-log sar plot with the empirical
and theoretical SARs plotted and csv files containing all of the data used to
make the plot.

'''
if __name__ == '__main__':

    import logging
    from macroeco.utils.workflow import Workflow
    from macroeco.empirical import Patch
    import macroeco.compare as comp
    from macroeco.output import SAROutput

    wf = Workflow(clog=True, svers=__version__)
    
    for data_path, output_ID, params in wf.single_datasets():
        try:
            params['div_list'].index((1,1))
        except:
            logging.info("Adding base area to parameter 'div_list': (1,1)")
            params['div_list'].append((1,1))
        sad_criteria = params['sar_criteria']
        for nm in params['div_cols']:
            sad_criteria[nm] = 'whole'
            sad_criteria[nm] = 'whole'

        patch = Patch(data_path, subset=params['subset'])
        sad = patch.sad(sad_criteria)
        sar = patch.sar(params['div_cols'], params['div_list'],
                                                    params['sar_criteria'])
        cmpr = comp.CompareSARCurve([sar], params['curve_list'],
                                                        [sad[1][0][1]])
        srout = SAROutput(output_ID)
        srout.plot_sars(cmpr.compare_curves(), names=params['names'])
        logging.info('Completed analysis %s\n' % output_ID)
    logging.info("Completed 'compare_sar' script")




        








