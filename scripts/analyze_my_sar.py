#!/usr/bin/python

'''Script the uses reproducible worlflow to compare the observed and predicted
sar for given data sets.

NOTES


'''

from macroeco.workflow import Workflow
from macroeco import output
from macroeco import sar_analysis

__author__ = "Mark Wilber"
__copyright__ = "Copyright 2012, Regents of the University of California"
__license__ = None
__version__ = "0.5"
__maintainer__ = "Mark Wilber"
__email__ = "mqw@berkeley.edu"
__status__ = "Development"

wf = Workflow()
wf.logger.debug("Analyzing SAR for given datasets")

for datafile, outputID, params in wf.single_datasets():
    if params['gridded'] == 'T' or params['gridded'] == 'True':
        obs_pred_sad = sar_analysis.get_obs_pred_sar(datafile, eval(params['grid']),\
                                                                        gridded=True)
        output.plot_sar(obs_pred_sad, outputID, params={}, interactive=wf.runs.interactive)

wf.logger.debug("SAR analysis complete")

