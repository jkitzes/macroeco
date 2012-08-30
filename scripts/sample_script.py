#!/usr/bin/python

'''
Sample script template showing use of Workflow object to enable reproducible 
workflow.

To Use
------
1.  Create an output directory to hold results of analysis
2.  In this directory, create a parameters.xml file
3.  Run this script with output directory as only sys arg

Best Practices
--------------
Required Parameters - The workflow object contains validation procedures that 
will check to ensure that the parameters file in the output directory contains 
all the parameters in a required_parameters dictionary in which entries are 
'parameter_name':'short_description'. This dictionary is optional if the script 
is run from the command line but required if it is run from the web GUI, in 
which case it is used to allow users to create a params file from the GUI.

Logging - To log the progress of the script, use logging.info('My message'), 
where info may be replaced by other words to signify logging level (see Python 
docs). By default, debug will print only to console and other levels will write 
to logfile.txt in the output directory.

GUI Compatibility - For compatibility with the html GUI, the strings 
'gui_name', 'short_descrip' and 'long_descrip' should appear in this script.

'''

__author__ = "Justin Kitzes"
__copyright__ = "Copyright 2012, Regents of the University of California"
__credits__ = ["John Harte"]
__license__ = None
__version__ = "0.1"
__maintainer__ = "Justin Kitzes"
__email__ = "jkitzes@berkeley.edu"
__status__ = "Development"

gui_name = '''Sample Script'''

short_descrip = '''This is a one sentence description of script actions.'''

long_descrip = '''This is a longer description of script actions, perhaps as 
long as a paragraph, that gives a detailed description of the purpose and 
operation of this script.'''

required_params = {'x': 'A sample numeric value'}

if __name__ == 'main':

    # TODO: Can future division be imported here, or needs to be at start?
    import logging
    from macroeco.utils.workflow import Workflow

    # Begin by creating a Workflow object
    wf = Workflow()

    # Loop through each dataset specified in params and run analysis. If
    # data_path is not a parameter, the loop below will run once with an empty
    # string for data_path.
    for data_path, output_ID, params in wf.single_datasets():

        y = params['x'] + 1

        with open(wf.output_path + output_ID + '_results.txt', 'w') as file:
            file.write('Results\n\n')
            file.write('Data Path: %s\n\n' % data_path)
            file.write('Parameters: %s\n\n' % str(params))
            file.write('y = %f\n\n' % y)

        logging.info('Completed analysis %s' % output_ID)
