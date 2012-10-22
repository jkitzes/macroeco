#!/usr/bin/python

'''
Script to compare sars
'''

__author__ = "Mark Wilber"
__copyright__ = "Copyright 2012, Regents of the University of California"
__credits__ = ["John Harte"]
__license__ = None
__version__ = "0.1"
__maintainer__ = "Mark Wilber"
__email__ = "mqw@berkeley.edu"
__status__ = "Development"

class global_str:
    subset = '''You should examine the columns in your data set and decide if you
	would like to subset your data in some particular way before the analysis
	begins. It is important to note that only the subsetted data will be analyzed.
	For example,  if you have a column named 'year' in your data set with values
	1998, 1999, and 2000 and you only want to look at the year 2000 for a
	particular analysis, you should select the == operator from the drop down list
	and type 2000 in the value field.  Similarly, you could use <, >, <=, >=, or
	!='''

    criteria = '''You should examine the columns in your dataset and decide if you
	would like to divide the data in a particular way for this analysis. For
	example, if the you have a spatial dataset with x,y coordinates and you are
	interested in examining macroecological metrics for two separate halves of your
	plot along the x coordinate, you could cut the x coordinate in two halves by
	giving the 'x' column a value of 2.  If the column that you would like to
	divide contains discrete values (e.g. year), you could enter the keyword
	'split' and each unique value will be analyzed separately. Conversely, the
	value 'whole' could be given to specify the entire column.  The value 'whole'
	is equivalent to 1 or leaving the value blank.\n\n

	There are four special words that can be used on a given column: 'species',
	'energy', 'count', and 'mass'.  When assigned to a column in your data set, the
	special word 'species' indicates the column that contains your species IDs, the
	special word 'energy' indicates the column that contains some type of energy
	measure, the special word 'mass' indicates a column that contains some type of
	mass measure, and the special word 'count' indicates the column that contains
	your species counts.  In the GUI, these special words can be chosen from the
	dropdown menu next to each column header. The special word 'species' MUST be
	assigned for every analysis.  If the special word 'count' is not assigned, the
	species counts are all assumed to be one.\n\n'''

    rarity_measure = '''This parameter allows you to specify the counts that
	you will consider rare.  If, for example, you want to know how many species in
	your plot have an abundance of 2 or less you would set this parameter to 2. If
	you enter more then one value, each value will be examined. Example input: [2]
	or [2, 5]. The brackets MUST be included.'''

    SAD_distributions = ''' 'logser','logser_ut', 'logser_ut_appx', 'plognorm_lt',
'nbd_lt', 'geo_ser', 'broken_stick', 'lognorm' '''

    SSAD_distributions = ''' 'nbd', 'binm', 'tgeo', 'fgeo', 'fnbd', 'pois' '''

gui_name = ''' Analyze Species-Area Relationships'''

summary = '''Compares a dataset's observed species-area relationships against 
theoretical species-area relationships'''

explanation = '''This script takes in a dataset(s) and list of curves to which
compare the observed datasets SAR will be compared.  The required parameters
for the script are the following:

'subset' : How one would like to initially subset his data (see DataTable 
class docstring). 

'div_cols' : A tuple of specifying the columns that will be divided during the
sar analysis. e.g. ('x', 'y')

'div_list' : A list of tuples, in the same order as the columns in div_list, 
that specify the divisions one would like to make in order to generate the sar.
A division of (1,1) indicates the whole plot.

'sar_criteria' : This parameter is a dictionary that specifies which column in
the data set is the column of species counts and which column is the column
of species names. E.g. {'spp' : 'species', 'spp_count' : 'count'}. In this case,
column name 'spp' is the species column and column 'spp_count' is the count
column. See Patch method sar() for more info.

'curve_list' : A list with the name of the SARCurve objects that will be
compared to the empirical SAR.  If none are given, the empirical SAR is just
plotted on its own.

'names' : List with the desired name of the plot.

For each dataset, this script generates a log-log sar plot with the empirical
and theoretical SARs plotted and csv files containing all of the data used to
make the plot.

'''

subset = global_str.subset

criteria = global_str.criteria

required_params = {'subset' : 'Initial subsetting of the data','div_cols' : 
                   'Tuple of column names to divide', 'criteria' :
                   'Dictionary with sar criteria', 'curve_list' : 'List of' +\
                   ' curves to compare', 'names' : 'List of plot titles',
		   'div_list': 'List of pairs of division values'}

if __name__ == '__main__':

    import logging
    from macroeco.utils.workflow import Workflow
    from macroeco.empirical import Patch
    import macroeco.compare as comp
    from macroeco.output import SAROutput

    wf = Workflow(required_params=required_params, clog=True, 
                                                            svers=__version__)
    
    for data_path, output_ID, params in wf.single_datasets():
        try:
            params['div_list'].index((1,1))
        except:
            logging.info("Adding base area to parameter 'div_list': (1,1)")
            params['div_list'].append((1,1))
        sad_criteria = params['criteria']
        for nm in params['div_cols']:
            sad_criteria[nm] = 'whole'
        
        patch = Patch(data_path, subset=params['subset'])
        sad = patch.sad(sad_criteria)
        sar = patch.sar(params['div_cols'], params['div_list'],
                                                    params['criteria'])
        cmpr = comp.CompareSARCurve([sar], params['curve_list'],
                                                    [sad[0][1]], patch=True)
        srout = SAROutput(output_ID)
        srout.plot_sars(cmpr.compare_curves(), names=params['names'])
        logging.info('Completed analysis %s\n' % output_ID)
    logging.info("Completed 'compare_sar' script")




        








