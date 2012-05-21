#!/usr/bin/python

'''Script uses reproducible workflow to get summary values from given datasets
and print them to a .txt file.

Parameters in parameters.xls
-----------------------------
Required:

name='distr' value=list of distribution names

    the list specifies which distributions to compare
    to empirical data.
    examples:
    name='distr' value="['mete', 'neg_binom']"

name='grids' value=list of tuples
    
    the list of tuples specifies the divisions of the plot
    examples:
    name='grids' value='[(16,16), (4,4), (2,2)]'

'''

from macroeco.workflow import Workflow
from macroeco import empirical
from macroeco import graph
from macroeco import predict_sad
from macroeco import sad_analysis


__author__ = "Mark Wilber"
__copyright__ = "Copyright 2012, Regents of the University of California"
__license__ = None
__version__ = "0.5"
__maintainer__ = "Mark Wilber"
__email__ = "mqw@berkeley.edu"
__status__ = "Development"

wf = Workflow()
wf.logger.debug('Getting values from given dataset(s)')

for datafile, outputID, params in wf.single_datasets():
    sad = sad_analysis.get_gridded_sad_list(datafile, eval(params['grid']), clean=True)
    for i in xrange(len(sad)):
        for m in xrange(len(sad[i])):
            ID = outputID + '_grid_[' + str(i) + '][' + str(m) + ']'
            graph.write_summary_table(sad[i][m], ID, eval(params['distr']), params)


    


wf.logger.info("Values obtained")


