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

summary = '''Compares a datasets observed sad(s) against theoretical sads'''

explantion = '''Blank for now'''

if __name__ == '__main__':

    import logging
    from macroeco.utils.workflow import Workflow
    from macroeco.empirical import Patch
    import macroeco.compare as comp
    from macroeco.output import SADOutput

    wf = Workflow()
    
    for data_path, output_ID, params in wf.single_datasets():
        patch = Patch(data_path, subset=params['subset'])
        sad = patch.sad(params['criteria'])
        cmpr = comp.CompareDistribution(sad, params['dist_list'], patch=True,
                                                                clean=True)
        sout = SADOutput(output_ID)
        sout.write_summary_table(cmpr.summary(), criteria=cmpr.criteria)
        sout.plot_rads(cmpr.compare_rads(), criteria=cmpr.criteria)
        logging.info('Completed analysis %s' % output_ID)




        







