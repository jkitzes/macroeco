#!/usr/bin/python

'''This script will allow users to convert transect data into the columnar
form'''

__author__ = "Mark Wilber"
__copyright__ = "Copyright 2012, Regents of the University of California"
__credits__ = ["John Harte"]
__license__ = None
__version__ = "0.1"
__maintainer__ = "Mark Wilber"
__email__ = "mqw@berkeley.edu"
__status__ = "Development"

ds = ''' Optional. Default: '''

information_about_stops = '''temp'''
delimiter = '''temp'''
replace = '''temp'''

# Columnar parameter descriptions
columns_to_split = '''temp'''
change_column_names = '''temp'''
add_column_names_and_values = '''temp'''
names_of_columns_to_be_removed = '''temp'''
how_and_where_to_fractionate = '''temp'''
merge_data = '''temp'''
subset = '''temp'''

required_params = {'information_about_stops' : information_about_stops}
optional_params = {'delimiter' : (delimiter + ds,
                    [',']), 'replace' : (replace + ds, None),
                    'columns_to_split' : (columns_to_split + ds, None),
                    'change_column_names' : (change_column_names + ds, (None,
                    None)), 'add_column_names_and_values' :
                    (add_column_names_and_values + ds, (None, None)),
                    'names_of_columns_to_be_removed' :
                    (names_of_columns_to_be_removed + ds, None),
                    'how_and_where_to_fractionate' :
                    (how_and_where_to_fractionate + ds, (None, None, None)),
                    'merge_data' : (merge_data + ds, 'No'), 'subset' : (subset
                    + ds, {})}


if __name__ == '__main__':

    import logging
    from macroeco.utils.workflow import Workflow
    import macroeco.utils.format_data as form

    wf = Workflow(required_params=required_params,
                 optional_params=optional_params, clog=True, svers=__version__)
    
    for data_paths, output_IDs, params, run_name, script_name in\
                                                             wf.all_datasets():
    
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

        columnar_obj.add_fields_to_data_list(params['add_column_names_and_values'][0],
                                             params['add_column_names_and_values'][1])

        columnar_obj.subset_data(params['subset'])

        columnar_obj.remove_columns(params['names_of_columns_to_be_removed'])

        columnar_obj.fractionate_data(params['how_and_where_to_fractionate'][0]
                                    , params['how_and_where_to_fractionate'][1]
                                   , params['how_and_where_to_fractionate'][2])

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

            

    logging.info("Completed 'convert_transect_data.py' script")


        
