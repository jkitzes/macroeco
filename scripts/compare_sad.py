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

#TODO: Need to fill out docstrings and all that good stuff
gui_name = '''SAD Analysis'''

summary = '''Compares a dataset's observed sad(s) against theoretical sads'''

explantion = '''This script takes in a dataset(s) and a list of distributions
to which the observed dataset(s)'s sads will be compared.  The required
parameters for this script are the following: 

'subset' : How one would like toinitially subset his data (see DataTable 
class docstring). 

'criteria' : How one would like to divide her dataset when caluculating the 
sad(s). The theoretical distributions are compared against every sad that is 
generated from the cutting specified in this parameter.  See Patch class 
docstring for details

'dist_list' : The list of distributions to which one could compare her observed
sad(s).  The full list of distributions can be found in distributions.py

For each sad generated by the divisions specified in criteria, this script
outputs a summary file, a rank-abundance plot, and data file with the
rank-abundance distributions used to build the rank-abundance plot.'''

if __name__ == '__main__':

    import logging
    from macroeco.utils.workflow import Workflow
    from macroeco.empirical import Patch
    import macroeco.compare as comp
    from macroeco.output import SADOutput

    wf = Workflow(clog=True)
    
    for data_path, output_ID, params in wf.single_datasets():
        patch = Patch(data_path, subset=params['subset'])
        sad = patch.sad(params['criteria'])
        cmpr = comp.CompareDistribution(sad, params['dist_list'], patch=True,
                                                                clean=True)
        sout = SADOutput(output_ID)
        sout.write_summary_table(cmpr.summary(), criteria=cmpr.criteria)
        sout.plot_rads(cmpr.compare_rads(), criteria=cmpr.criteria)
        sout.plot_cdfs(cmpr.compare_cdfs(), cmpr.data_list,
                        criteria=cmpr.criteria)
        logging.info('Completed analysis %s' % output_ID)




        







