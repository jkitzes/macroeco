#!/usr/bin/python

'''This python file contains global strings used in the scripts.  Consolidated
in this script for easy maintenance'''

subset = '''\nYou should examine the columns in your data set and decide if you
would like to subset your data in some particular way before the analysis
begins. It is important to note that only the subsetted data will be analyzed.
For example,  if you have a column named 'year' in your data set with values
1998, 1999, and 2000 and you only want to look at the year 2000 for a
particular analysis, you should select the column year from left-hand most
dropdown list, select the == operator from the operator dropdown list and type
2000 in the value field.  Similarly, you could use <, >, <=, >=, or != with any
column and value in your data.'''

criteria = '''\nYou should examine the columns in your dataset and decide if you
would like to divide the data in a particular way for this analysis.  For
example, if you have a spatial dataset with x,y coordinates and you are
interested in examining macroecological metrics for two separate halves of your
plot along the x coordinate, you could cut the x coordinate in two halves by
giving the 'x' column a value of 2.  

If the column that you would like to divide contains discrete values (e.g.
year), you could enter the keyword 'split' and each unique value will be
analyzed separately.  Conversely, the value 'whole' could be given to specify
the entire column.  The value 'whole' is equivalent to 1 or leaving the value
blank. If you would like to divide a given column, please select the word
'division' from the dropdown menu and input a value as discussed above.\n\n

There are four other special words that can be used on a given column:
'species', 'energy', 'count', and 'mass'.  When assigned to a column in your
data set, the special word 'species' indicates the column that contains your
species IDs, the special word 'energy' indicates the column that contains some
type of energy measure, the special word 'mass' indicates a column that
contains some type of mass measure, and the special word 'count' indicates the
column that contains your species counts.  These special words can be chosen
from the dropdown menu next to each column header. The special word 'species'
MUST be assigned for every analysis.  If the special word 'count' is not
assigned, the species counts are all assumed to be one.\n\n
    
If there are columns in your data that are not relevant for this analysis leave
the value in the dropdown box as 'NA'.  Columns designated 'NA' will not
influence the analysis.\n\n'''


rarity_measure = '''\nThis parameter allows you to specify the counts that you
will consider rare.  If, for example, you want to know how many species in your
plot have an abundance of 2 or less you would set this parameter to 2. If you
enter more then one value, each value will be examined. Example input: [2] or
[2, 5]. The brackets MUST be included.'''

SAD_distributions = '''\n 
'logser' : Fisher's logseries distribution;
'logser_ut' : Upper-truncated logseries derived from MaxEnt;
'logser_ut_appx' : Approximation for the upper-truncated logseries;
'lognorm' : Lognormal distribution;
'plognorm_lt' : Poisson lognormal distribution with 0 truncated;
'nbd_lt' : Negative binomial distribution with 0 truncated;
'geo_ser' : Geometric series distribution;
'broken_stick' : McArthur's broken stick distribution '''

SSAD_distributions = '''\n 
'binm' : Binomial distribution;
'pois' : Poisson distribution;
'nbd' : Negative binomial distribution;
'fnbd' : Finite-negative binomial distribution;
'geo' : Geometric distribution;
'fgeo' : Finite-geometric distribution;
'tgeo' : Truncated geometric distrbituion derived from MaxEnt'''

short_subset = '''\nSpecifications for how you want to subset your data before the
analysis. Note that only the subsetted data will be included in the analysis.
The left-hand dropdown box contains all the columns of your dataset and you may
choose one or more to subset. Please see analysis explanation for more detail
and examples.'''

short_criteria = '''\nSpecifications for how you want to divide your data during
the analysis. The words you see below are the shared columns of your
dataset(s).  You must designate your species column with the special word
'species' found in the dropdown menu. You are not required to fill any
additional columns for this analysis. Please see analysis explanation for more
detail and examples.'''

optional = ''' Optional parameter. Default value:  '''

req = '''Required parameter. '''

#### Formatting strings #### 

explanation_string = '''This formatting script loads {0} datasets and
reformats them into columnar data using the parameters that you specify below.
We define columnar data as a dataset that has distinct column headers and has
rows that describe the attributes of a single entity (often a species).  For
example, a row could describe the location spatial location of a species, the
total number of individuals of that species at that spatial location,
attributes about that location, the date the species was censuses, etc. All of
these atttributes are specified by the column headers. Please see the website
http://www.ctfs.si.edu/plots/summary/ for examples of columnar data.

'''

output_string = '''This formatting script outputs a formatted csv data file to
specified folder within ../macroeco/data/formatted.  You can specify the name
of the output formatted file(s).  If you do not, the script will hard code them
with the script name, run name, and some appended string.

'''

process_string = '''
The formatting process is as follows: 

1. The specified {0} data is loaded\n
2. Any {0} data-specific formatting parameters are applied to the {0}
data\n
3. The {0} data is converted into columnar data\n
4. Any columnar data-specific formatting parameters are applied to the columnar
data\n
5. The columnar data is output\n
'''

delimiter = '''\nThe file delimiter used in the data files.

Example input:

1. [','] 

Where ',' is the file delimiter.

2. ['+'] 

Where '+' is the file delimiter.

The brackets and quotes MUST be include'''

missing_values_from_a_given_column = '''\nSpecifies what is a
missing data value in any given column in the data set.  The input must be
formatted as a pythonic dictionary.

Example input:  

1. {'count' : 'NA', 'year' : ''}

This input says that a the data column 'count' has missing values 'NA' and the
data column 'year' has missing values '' (blank). The brackets and semicolons
are required for this parameter'''


delete_missing_values = '''\nEither True or False.  If True, the missing values
specified in the missing_values_from_a_given_column parameter are removed from
the formatted data (your archival data remains unchanged). If False, only NaN
values are removed from the formatted data.

Chose either: True or False.'''

subset = '''\nA permanent subset to the formatted data, {'column_name':
'condition'}, which will limit all analysis to records in which column_name
meets the condition. The condition is a formatted as ('comparison operator',
'value').  Possible comparison operators are '==', '!=', '<, '>', '<=', '>='.  
Please note that your archival data will remain unchanged. 

Subsetting examples:

1. {'year': ('==' ,  2005), 'x': [('>' ,  20), ('<' ,  40)]}

Restricts analysis to year 2005 and x values between 20 and 40. Note that for
multiple conditions for a column square brackets MUST be included 
(i.e. x : [('>', 20), ('<', 40)]).  For a single condition on a column they are
optional (i.e. 'year': ('==', 2005)). 

2. {'name' : ('==', 'John')}

Includes only rows in which column 'name' equals 'John'. When subsetting on a
string, the string should be quoted (i.e. ('==', 'John')) '''

columns_to_split = '''\nUse this if you want to split your single dataset into
multiple datasets based on given column names.  For example, if you have a 
dataset with column names ('x1', 'x2', 'x3','g', 'h') and you want to make
three datasets with column names ('x1', 'g', 'h'), ('x2', 'g', 'h'), and ('x3',
'g', 'h') you could type ['x1', 'x2', 'x3'] and your single data would be made
into three datasets with the columns given above. Notice that ALL columns that
are not specified are included in each new dataset.

Example input:

1. ['x1', 'x2', 'x3'] OR [('x1',) ('x2',), ('x3',)]

Makes three datasets where each one contains only one of the specified columns.
All columns that are not specified are included.  The brackets ([]) MUST be
included.

2. [('x1', 'y1'), ('x2', y2'), ('x3', 'y3'), ('x4', 'y4')]

Makes four datasets where each data set contains only one of the above pairs
x,y pairs.  For example, the first data set would have columns ('x1', 'y1', ...
all unspecified columns) but it would not have columns 'x2', 'x3', 'x4, 'y2',
'y3', or 'y4'. '''

change_column_names = '''\nSpecifies the column names that you wish to change and
the names that you wish to change them to.  This parameter is useful if you
wish to merge data sets.

Example input:

1. (['fred', 'jane'], ['mark', 'mary']) or ['fred', 'jane'], ['mark', 'mary'] 

Changes column 'fred' to 'mark' and column 'jane' to 'mary' in all datasets.
The brackets are required.

2. ([('x1', 'x2', 'x3'), 'h1'], ['x', 'h']) 

Changes columns 'x1', 'x2', 'x3' to 'x' and column 'h1' to 'h'.  All
brackets are required.'''

add_column_names_and_values = '''\nSpecifies additional columns that you want to
add to the data and the values the column will take for each dataset. 

Example input:

1. {'year' : (1998, 1999, 2000), 'name' : ('Fred', 'George', 'Ron')}

Adds the column 'year' and 'name' to all datasets.  In this example, there are
three data sets and the values of 'year' for the first, second, and third
dataset are set to 1998, 1999, and 2000, respectively.  Simiarly, the values of
column 'name' for the first, second, and third dataset are set to 'Fred',
'George', and 'Ron', respectively. The length of values to be assigned (i.e.
(1998, 1999, 2000)) MUST equal the number of datasets. All brackets and
punctuation must be included

2. {'year' : (1998,)}

Adds the columns 'year' with a value of 1998 to the one and only dataset being
considered.

'''
names_of_columns_to_be_removed = '''\nRemove any number of columns from the
dataset by specifying the column names.

Example Input:

1. 'name'

Removes the column 'name' from all data sets

2. ['name', 'species', 'date']

Remove the columns 'name', 'species', and 'date' from all data sets
'''

merge_data = '''\nEither Y/yes or N/no.  If Y/yes, attempts to merge all of the
data into one dataset.  If the merge is successful, only the single merged data
file will be output. If the merge cannot be completed an error will be
displayed.  If N/no, no merge will be attempted and all datasets will be
output.'''

columnar_params_full =\
'''
*** delimiter ***

{0}

*** missing_values_from_a_given_column ***

{1}

*** delete_missing_values *** 

{2}

*** columns_to_split ***

{3}

*** change_column_names ***

{4}

*** add_column_names_and_values ***

{5}

*** names_of_columns_to_be_removed ***

{6}

*** merge_data ***

{7}

*** subset ***

{8}

'''.format(delimiter, missing_values_from_a_given_column,
delete_missing_values, columns_to_split, change_column_names,
add_column_names_and_values, names_of_columns_to_be_removed, merge_data,
subset)

columnar_params_med =\
'''
*** delimiter ***

{0}

*** columns_to_split ***

{1}

*** change_column_names ***

{2}

*** add_column_names_and_values ***

{3}

*** names_of_columns_to_be_removed ***

{4}

*** merge_data ***

{5}

*** subset ***

{6}

'''.format(delimiter, columns_to_split, change_column_names,
add_column_names_and_values, names_of_columns_to_be_removed, merge_data,
subset)

columnar_params_small =\
'''
*** columns_to_split ***

{0}

*** change_column_names ***

{1}

*** add_column_names_and_values ***

{2}

*** names_of_columns_to_be_removed ***

{3}

*** merge_data ***

{4}

*** subset ***

{5}

'''.format(columns_to_split, change_column_names,
add_column_names_and_values, names_of_columns_to_be_removed, merge_data,
subset)




def  check_columnar_params(params, script):
    '''This function checks that all of the parameters required to convert
    columnar data have the correct types. This test does not completely
    validate parameters. Just check the first level type.
    
    Parameters
    ----------
    params : dict
        Parameter dictionary
    script : str
        Either 'grid', 'dense', 'columnar', or 'transect'.
    
    '''
    
    # Can't check names_of_columns_to_be_removed because it can be a string.
    if script == 'grid':
        prms_types = [('columns_to_split', type([])),
                  ('change_column_names', type((2,))),
                  ('add_column_names_and_values', type({})),
                  ('merge_data', str),
                  ('subset', type({}))]

    elif script != 'columnar':
        prms_types = [('delimiter' , type([])),
                  ('columns_to_split', type([])),
                  ('change_column_names', type((2,))),
                  ('add_column_names_and_values', type({})),
                  ('merge_data', str),
                  ('subset', type({}))]

    else:
        prms_types = [('delimiter' , type([])),
                  ('missing_values_from_a_given_column', type({})),
                  ('delete_missing_values', type(True)),
                  ('columns_to_split', type([])),
                  ('change_column_names', type((2,))),
                  ('add_column_names_and_values', type({})),
                  ('merge_data', str),
                  ('subset', type({}))]
    
    for i, pair in enumerate(prms_types):

        if type(params[pair[0]]) != pair[1]:
            if params[pair[0]] != None:
                raise TypeError("Parameter '%s' must be a %s not a %s." % (pair[0],
                        str(pair[1]), str(type(params[pair[0]]))) + 
                        " Please check the formatting of '%s': %s " % (pair[0], 
                        str(params[pair[0]])))






