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

name='grid' value=list of tuples
    
    the list of tuples specifies the divisions of the plot
    examples:
    name='grid' value='[(16,16), (4,4), (2,2)]'

Output
------
Per data set, per run, this analysis outputs as many .txt files as there are cells
made by gridding the given data plot by the values specified in 'grid'.  For example,
if 'grid' is "[(1,1)]", this analysis outputs one .txt files that contains the empirical
values for the entire plot (see sad_analysis.get_values_for_sad) and predicted values
for all the distributions specified in 'distr'.  These are all contained in one .txt file.
If 'grid' is something like "[(2,2),(4,4),(6,6)]", seperate .txt files are made for each tuple 
in the list and each cell within each tuple.  Within each tuple (x,y) in 'grid', 
there are x*y .txt files made. In the output name, one will see something that looks like *_grid_[i][j].  
i represents the corresponding tuple in 'grid' where j represents the cell number in the 
gridded plot, counting down columns then across rows (where 0 is the first entry).  
So for example, if 'grids' is "[(2,2), (4,4)]" the file *_grid_[1][13] corresponds to 
the tuple (4,4) and the 14th cell in the grid counting down columns and across rows.  
The 14th cell corresponds to cell (3,0). Again, each .txt files contains information for
all distributions specified in 'distr' for the given tuple and given cell.

'''

from macroeco.workflow import Workflow
from macroeco import empirical
from macroeco import output
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
            output.write_summary_table(sad[i][m], ID, eval(params['distr']), params)


    


wf.logger.info("Values obtained")


