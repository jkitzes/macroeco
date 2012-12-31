#!/usr/bin/python

'''This format type will allow users to convert columnar data into a new form of
columnar data'''

__author__ = "Mark Wilber"
__copyright__ = "Copyright 2012, Regents of the University of California"
__credits__ = ["John Harte"]
__license__ = None
__version__ = "0.1"
__maintainer__ = "Mark Wilber"
__email__ = "mqw@berkeley.edu"
__status__ = "Development"

import macroeco.utils.global_strings as gb

gui_name = '''Convert Columnar Data'''

summary = '''Converts and formats columnar data'''

explanation = '''
FORMATTING DESCRIPTION

{0} 

PROCESS

The formatting process is as follows:

1. The specified columnar data is loaded\n
2. Any columnar data-specific formatting parameters are applied to the columnar
data\n
3. The columnar data is output

 
OUTPUT

{1}

PARAMETERS

{2}
'''.format(gb.explanation_string.format('columnar'), gb.output_string,
                                                gb.columnar_params_full)

required_params = {}
optional_params = {'delimiter' : (gb.optional + gb.delimiter,
                    [',']), 'missing_values_from_a_given_column' : (gb.optional
                    + gb.missing_values_from_a_given_column , None),
                    'delete_missing_values' : (gb.optional +
                    gb.delete_missing_values, False), 'columns_to_split' :
                    (gb.optional + gb.columns_to_split, None),
                    'change_column_names' : (gb.optional +
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

        # Check that columnar parameters are generally formatted correctly
        gb.check_columnar_params(params, 'columnar')
    
        columnar_obj = form.Columnar_Data(data_paths, params['delimiter'][0],
                                  params['missing_values_from_a_given_column'],
                                  params['delete_missing_values'],
                                  archival=False)

        # The order in which we do these operations can make a difference
        columnar_obj.split_up_data_by_field(params['columns_to_split'])

        columnar_obj.change_column_names(params['change_column_names'][0],
                                         params['change_column_names'][1])

        columnar_obj.subset_data(params['subset'])

        columnar_obj.add_fields_to_data_list(params['add_column_names_and_values'])

        columnar_obj.remove_columns(params['names_of_columns_to_be_removed'])

        for data_path in data_paths:
            logging.info('Formatted the columnar data %s' % 
                         data_path)
        
        # Merge data into one data file
        if params['merge_data'] == 'Yes' or params['merge_data'] == 'yes':

            try:

                columnar_obj.output_merged_data('{0}_{1}_merged_data'.
                                                format(script_name, run_name))
                logging.info("Merged and saved all data in run '%s' to file" % run_name
                   + ' {0}_{1}_merged_data.csv'.format(script_name, run_name))
            except:
                logging.warning("Failed to merge data. Column names and/or types"
                             + " do not match. Outputting each dataset"
                             + " separately.")
                params['merge_data'] = 'No'

        if params['merge_data'] == 'No' or params['merge_data'] == 'no':
            for out in [ot + '_formatted_columnar' for ot in output_IDs]:
                logging.info('Saving columnar data as %s' % out)

            columnar_obj.output_columnar_data([output_ID +
                        '_formatted_columnar' for output_ID in output_IDs])

    logging.info("Completed 'convert_columnar_data.py' format type")


