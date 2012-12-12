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

gui_name = '''Convert Grid Data'''

summary = '''Converts and formats grid data'''

# Grid parameter descriptions
truncation_symbols = '''\n

Truncate the grid at and after the first occurence of a given symbol.

Example input:

1. '%'

The contents of all cells in the data set are truncated at and after '%'. 

'''

remove_replace_values = '''\n

A Grid Data specific parameter. Remove a given value from every cell and
replace it with another value.

Example Input:

1. [('species1', 'GYSPHY')]

Remove the string 'species1' and replace it with the string 'GYSPHY'

2. [('missing', ''), ('f', 't')]

Remove the string 'missing' and replace with nothing. Remove all 'f' and
replace with 't'.

The parantheses, brackets, and quotation marks are required.'''

char_btwn_species_and_count = '''\n

A Grid Data specific parameter. The character separating a species name from 
its count. 

Example Input:

1. ['-'] 

The character '-' separates a species from its count. For example,
PACNEO - 1; the count of species PACNEO is 1.

2. [',']

The character ',' separates a species from its count. For example,
PACNEO , 1; the count of species PACNEO is 1'''

char_btwn_species = '''\n

A Grid Data specific parameter. The character that separates two species.

Example Input:

1. ['/\/n']

The character /\/n (newline character) separates two species. For example,
PACNEO - 1/\/nATRTRI - 5

2. ['+']

The character '+' separates two species. For example,
PACNEO - 1+ARTTRI - 5
'''

explanation = '''FORMATTING DESCRIPTION

{5} 

We define grid data as data that has the same physical layout as a census grid,
but in digital format. So, for example, a 16 x 16 census grid would be
represented by a 16 x 16 csv file where each cell contains properly formatted
counts. See Documentation for a complete description of grid data.

PROCESS

{7}

OUTPUT

{6}

PARAMETERS

*** truncation_symbols ***

{0}

*** remove_replace_values ***

{1}

*** char_btwn_species_and_count **

{2}

*** char_btwn_species ***

{3}

{4}
'''.format(truncation_symbols, remove_replace_values,
char_btwn_species_and_count, char_btwn_species, gb.columnar_params_small,
gb.explanation_string.format('grid'), gb.output_string,
gb.process_string.format('grid'))

required_params = {}

optional_params = {'truncation_symbols' : (gb.optional + truncation_symbols, None), 
                    'remove_replace_values': (gb.optional +
                    remove_replace_values, [(None, None)]),
                    'char_btwn_species_and_count' : (gb.optional +
                    char_btwn_species_and_count, ['-']), 'char_btwn_species':
                    (gb.optional + char_btwn_species, ['\n']),
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
    
    # What about formatting for multiple data sets simultaneously?  In this
    # case it would be nice to have Workflow yield all the datasets in a given
    # run.
    for data_paths, output_IDs, params, run_name, script_name in\
                                                             wf.all_datasets():

        # Check the parameter formats. Raises error if improper formatting
        gb.check_columnar_params(params, 'grid') 

        # Each script handles one specific type of data.  This script deals
        # with gridded data.

        grid_data = form.Grid_Data(data_paths, archival=False)
        
        # Allowing user to truncated grid cells.  They can do it multple times.
        grid_data.truncate_grid_cells(params['truncation_symbols'])

        #  User can remove and replace multiple things.  This should be a list
        #  of tuples with each tuple having two elements. 
        for rm_rp in params['remove_replace_values']:
            grid_data.remove_and_replace(rm_rp[0], rm_rp[1])

        # Convert gridded data to columnar data
        grid_data.grid_to_dense(spacer=params['char_btwn_species_and_count'][0],
                                spp_sep=params['char_btwn_species'][0],
                                archival=False)
        num_spp = [len(spp_list) for spp_list in grid_data.unq_spp_lists]
        dense_obj = grid_data.Dense_Object
        dense_obj.dense_to_columnar(3, tuple(num_spp), archival=False)

        # Format Columnar Data
        columnar_obj = dense_obj.Columnar_Object

        columnar_obj.split_up_data_by_field(params['columns_to_split'])


        columnar_obj.change_column_names(params['change_column_names'][0],
                                         params['change_column_names'][1])

        columnar_obj.subset_data(params['subset'])

        columnar_obj.add_fields_to_data_list(params['add_column_names_and_values'])

        columnar_obj.remove_columns(params['names_of_columns_to_be_removed'])

        for data_path in data_paths:
            logging.info('Converted and formatted the grid data %s to' % 
                         data_path + ' columnar data')
        
        # Merge data into one data file
        if params['merge_data'] == 'Yes' or params['merge_data'] == 'yes':
            columnar_obj.output_merged_data('{0}_{1}_merged_data'.
                                                format(script_name, run_name))
            logging.info("Merged and saved all data in run '%s' to file" % run_name
                   + ' {0}_{1}_merged_data.csv'.format(script_name, run_name)) 
        else:
            columnar_obj.output_columnar_data([output_ID +
                        '_grid_to_columnar' for output_ID in output_IDs])

            for out in [ot + '_grid_to_columnar' for ot in output_IDs]:
                logging.info('Saving columnar data as %s' % out)

    logging.info("Completed 'convert_grid_data.py' script")
        




        




