#!/usr/bin/python

'''This script will allow users to format their data with predefined formatting
functions provided in format_data'''

required_params = {}

if __name__ == '__main__':

    import logging
    import macroeco.utils.workflow import Workflow
    import macroeco.utils.format_data as form

    wf = Workflow(required_params=required_params, clog=True,
                  svers=__verision__)

    for data_path, output_ID, params in wf.single_datasets():

