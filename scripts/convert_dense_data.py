#!/usr/bin/python

'''This script will allow users to convert gridded data into the columnar
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

gui_name = '''Convert Dense Data'''

summary = '''Converts and formats dense data'''

number_of_first_species_column = '''\n

A Dense Data specific parameter. The column number in the dataset where the
first species occurs. 

For example, if a dataset has the columns (x1, x1, spp1, spp2) and 'spp1' is
the first species, this parameter would have the value of 3.'''

number_of_species_in_census = '''\n

A Dense Data specific parameter. The total number of species in the dataset(s).

Example Input:

1. 34

Indicates there are 34 species in the single dataset or all datasets

2. [34, 56, 102]

Indicated that there are 34 species in the first dataset, 56 species in the
second dataset, and 102 species in the third dataset'''

replace_missing_with_value = '''\n

A Dense Data specific parameter. Specify a missing value in the dataset and a
value with which to replace it.

Example Input:

1. ('', 0)

Replace '' (blank data value) with 0

2. ('NA', 0)

Replace 'NA' with 0.'''

explanation = '''
FORMATTING DESCRIPTION

{3} 

We define dense data as data with column headers specifying a given cell (row,
column) and a column for each species in the data set. See Documentation for a
complete description of grid data.

PROCESS

{5}

OUTPUT

{4}

PARAMETERS

*** number_of_first_species_column ***

{0}

*** number_of_species_in_census ***

{1}

{2}
'''.format(number_of_first_species_column, number_of_species_in_census, 
           gb.columnar_params_med, gb.explanation_string.format('dense'),
           gb.output_string, gb.process_string.format('dense'))


required_params = {'number_of_first_species_column' :
                    gb.req + number_of_first_species_column, 
                    'number_of_species_in_census' :
                    gb.req + number_of_species_in_census}

optional_params = {'delimiter' : (gb.optional + gb.delimiter, [',']), 
                    'replace_missing_with_value': (gb.optional +
                    replace_missing_with_value , None), 'columns_to_split' :
                    (gb.optional + gb.columns_to_split, None),
                    'change_column_names' : (gb.optional +
                    gb.change_column_names, (None, None)),
                    'add_column_names_and_values' : (gb.optional +
                    gb.add_column_names_and_values, None),
                    'names_of_columns_to_be_removed' : (gb.optional +
                    gb.names_of_columns_to_be_removed, None),
                    'merge_data' : (gb.optional + gb.merge_data, 'No'),
                    'subset' : (gb.optional + gb.subset, {})}

if __name__ == '__main__':

    import logging
    from macroeco.utils.workflow import Workflow
    import macroeco.utils.format_data as form

    wf = Workflow(required_params=required_params,
                 optional_params=optional_params,
                 clog=True, svers=__version__, short_output_name=True)
    
    for data_paths, output_IDs, params, run_name, script_name in\
                                                             wf.all_datasets():

        gb.check_columnar_params(params, 'dense')

        dense_data = form.Dense_Data(data_paths, params['delimiter'][0],
                          params['replace_missing_with_value'], archival=False)
        
        if type(params['number_of_species_in_census']) == int:
            params['number_of_species_in_census'] =\
                                       (params['number_of_species_in_census'],)

        # Convert to Columnar Data
        dense_data.dense_to_columnar(params['number_of_first_species_column'] - 1,
                      params['number_of_species_in_census'], archival=False)

        # Format Columnar Data
        columnar_obj = dense_data.Columnar_Object

        columnar_obj.split_up_data_by_field(params['columns_to_split'])

        columnar_obj.change_column_names(params['change_column_names'][0],
                                         params['change_column_names'][1])

        columnar_obj.subset_data(params['subset'])

        columnar_obj.add_fields_to_data_list(params['add_column_names_and_values'])

        columnar_obj.remove_columns(params['names_of_columns_to_be_removed'])

        for data_path in data_paths:
            logging.info('Converted and formatted the dense data %s to' % 
                         data_path + ' columnar data')
        
        # Merge data into one data file
        if params['merge_data'] == 'Yes' or params['merge_data'] == 'yes':
            columnar_obj.output_merged_data('{0}_{1}_merged_data'.
                                                format(script_name, run_name))
            logging.info("Merged and saved all data in run '%s' to file" % run_name
                   + ' {0}_{1}_merged_data.csv'.format(script_name, run_name)) 
        else:
            columnar_obj.output_columnar_data([output_ID +
                        '_dense_to_columnar' for output_ID in output_IDs])

            for out in [ot + '_dense_to_columnar' for ot in output_IDs]:
                logging.info('Saving columnar data as %s' % out)

    logging.info("Completed 'convert_dense_data.py' script")



        
