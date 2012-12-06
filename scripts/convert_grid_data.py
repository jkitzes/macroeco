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

ds = ''' Optional. Default: '''

# Grid parameter descriptions
truncation_symbols = '''temp'''
remove_replace_values = '''temp'''
char_btwn_species_and_count = '''temp'''
char_btwn_species = '''temp'''

# Columnar parameter descriptions
columns_to_split = '''temp'''
change_column_names = '''temp'''
add_column_names_and_values = '''temp'''
names_of_columns_to_be_removed = '''temp'''
how_and_where_to_fractionate = '''temp'''
merge_data = '''temp'''
subset = '''temp'''

required_params = {}

optional_params = {'truncation_symbols' : (truncation_symbols + ds, None), 
                    'remove_replace_values': (remove_replace_values + ds,
                    [(None, None)]), 'char_btwn_species_and_count' :
                    (char_btwn_species_and_count + ds, ['-']),
                    'char_btwn_species': (char_btwn_species + ds, ['\n']),
                    'columns_to_split' : (columns_to_split + ds, None),
                    'change_column_names' : (change_column_names + ds, (None,
                    None)), 'add_column_names_and_values' :
                    (add_column_names_and_values + ds, None),
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
    
    # What about formatting for multiple data sets simultaneously?  In this
    # case it would be nice to have Workflow yield all the datasets in a given
    # run.
    for data_paths, output_IDs, params, run_name, script_name in\
                                                             wf.all_datasets():
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

        columnar_obj.fractionate_data(params['how_and_where_to_fractionate'][0]
                                    , params['how_and_where_to_fractionate'][1]
                                   , params['how_and_where_to_fractionate'][2])

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
        




        




