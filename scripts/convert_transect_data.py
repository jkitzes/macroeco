#!/usr/bin/python

'''This format type will allow users to convert transect data into the columnar
form'''

__author__ = "Mark Wilber"
__copyright__ = "Copyright 2012, Regents of the University of California"
__credits__ = ["John Harte"]
__license__ = None
__version__ = "0.1"
__maintainer__ = "Mark Wilber"
__email__ = "mqw@berkeley.edu"
__status__ = "Development"

import macroeco.utils.global_strings as gb

gui_name = '''Convert Transect Data'''

summary = '''Converts and formats transect data'''

information_about_stops = '''\n

A Transect Data specific parameter. Specifies where the stop columns begin
in the data, how many stop columns the data set contains, and the desired name
of the stop column in the formatted data. See explanation for description of
stops.

Example Input:

1.  (3, 3, 'stop')

Specifies that the stop columns begins in the third column and that there are
three stops total.  In addition, specifies that the stop column should be
called 'stop' in the formatted data

2. ([4, 5], [12, 34], 'pit')

Specifies that there are two data sets and in the first dataset the stop
columns begin at column 4 and in the second dataset the stop columns begin at
column 5.  The first data set has 12 stop columns and the second data set has
34 stop columns.  Finally the stop columns will be labeled 'pit' in the
formatted data.

Brackets must be included if you have more than one dataset'''

delimiter = gb.delimiter

replace = '''\n

A Transect Data specific parameter. Specifies a value that you would like to
replace and the value that you would like to replace it with

Example Input:

1.  ('', 0)

Replaces all '' (blank data) with 0

2. ('0', 'NA')

Replaces all 0 with NA. Quotes are required around the items in the
parentheses. '''

explanation = '''
FORMATTING DESCRIPTION

{4} 

We define transect data with given stops along a transect.  This is analogous
to a nested design in which each each stop is a random effect and each transect
is either a random or fixed effect. In this data format, each stop within a
transect has its own column and the species/item count at that given stop
within that given transect is given under each stop column.  A complete
description of transect data can be found in the Documentation.

PROCESS

{5}

OUTPUT

{6}

PARAMETERS

*** information_about_stops ***

{0}

*** delimiter ***

{1}

*** replace **

{2}

{3}
'''.format(information_about_stops, delimiter, replace, gb.columnar_params_med,
gb.explanation_string.format('transect'), gb.process_string.format('transect'),
gb.output_string)

required_params = {'information_about_stops' : gb.req + information_about_stops}
optional_params = {'delimiter' : (gb.optional + delimiter,
                    [',']), 'replace' : (gb.optional + replace, None),
                    'columns_to_split' : (gb.optional + gb.columns_to_split,
                    None), 'change_column_names' : (gb.optional +
                    gb.change_column_names, (None, None)),
                    'add_column_names_and_values' : (gb.optional +
                    gb.add_column_names_and_values, None),
                    'names_of_columns_to_be_removed' : (gb.optional +
                    gb.names_of_columns_to_be_removed, None), 'merge_data' :
                    (gb.optional + gb.merge_data, 'No'), 'subset' :
                    (gb.optional + gb.subset, {})}


if __name__ == '__main__':

    import logging
    from macroeco.utils.workflow import Workflow
    import macroeco.utils.format_data as form

    wf = Workflow(required_params=required_params,
                 optional_params=optional_params, clog=True, svers=__version__)
    
    for data_paths, output_IDs, params, run_name, script_name in\
                                                             wf.all_datasets():

        gb.check_columnar_params(params, 'transect')
    
        transect_data = form.Transect_Data(data_paths, \
                        delim=params['delimiter'][0],
                        replace=params['replace'], archival=False)

        # Convert transect data into columnar form
        transect_data.transect_to_columnar(
                                       params['information_about_stops'][0] - 1
                                         , params['information_about_stops'][1]
                                        , params['information_about_stops'][2])

        # Format Columnar Data
        columnar_obj = transect_data.Columnar_Object

        columnar_obj.split_up_data_by_field(params['columns_to_split'])

        columnar_obj.change_column_names(params['change_column_names'][0],
                                         params['change_column_names'][1])

        columnar_obj.add_fields_to_data_list(params['add_column_names_and_values'])

        columnar_obj.subset_data(params['subset'])

        columnar_obj.remove_columns(params['names_of_columns_to_be_removed'])

        for data_path in data_paths:
            logging.info('Converted and formatted the transect data %s to' % 
                         data_path + ' columnar data')
        
        # Merge data into one data file
        if params['merge_data'] == 'Yes' or params['merge_data'] == 'yes':
            columnar_obj.output_merged_data('{0}_{1}_merged_data'.
                                                format(script_name, run_name))
            logging.info("Merged and saved all data in run '%s' to file" % run_name
                   + ' {0}_{1}_merged_data.csv'.format(script_name, run_name)) 
        else:
            for out in [ot + '_transect_to_columnar' for ot in output_IDs]:
                logging.info('Saving columnar data as %s' % out)

            columnar_obj.output_columnar_data([output_ID +
                        '_transect_to_columnar' for output_ID in output_IDs])

            

    logging.info("Completed 'convert_transect_data.py' format type")


        
