#!/usr/bin/python

'''This python file contains global strings used in the scripts.  Consolidated
in this script for easy maintenance'''

subset = '''You should examine the columns in your data set and decide if you
would like to subset your data in some particular way before the analysis
begins. It is important to note that only the subsetted data will be analyzed.
For example,  if you have a column named 'year' in your data set with values
1998, 1999, and 2000 and you only want to look at the year 2000 for a
particular analysis, you should select the column year from left-hand most
dropdown list, select the == operator from the operator dropdown list and type
2000 in the value field.  Similarly, you could use <, >, <=, >=, or != with any
column and value in your data.'''

criteria = '''You should examine the columns in your dataset and decide if you
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


rarity_measure = '''This parameter allows you to specify the counts that you
will consider rare.  If, for example, you want to know how many species in your
plot have an abundance of 2 or less you would set this parameter to 2. If you
enter more then one value, each value will be examined. Example input: [2] or
[2, 5]. The brackets MUST be included.'''

SAD_distributions = ''' 
'logser' : Fisher's logseries distribution;
'logser_ut' : Upper-truncated logseries derived from MaxEnt;
'logser_ut_appx' : Approximation for the upper-truncated logseries;
'lognorm' : Lognormal distribution;
'plognorm_lt' : Poisson lognormal distribution with 0 truncated;
'nbd_lt' : Negative binomial distribution with 0 truncated;
'geo_ser' : Geometric series distribution;
'broken_stick' : McArthur's broken stick distribution '''

SSAD_distributions = ''' 
'binm' : Binomial distribution;
'pois' : Poisson distribution;
'nbd' : Negative binomial distribution;
'fnbd' : Finite-negative binomial distribution;
'geo' : Geometric distribution;
'fgeo' : Finite-geometric distribution;
'tgeo' : Truncated geometric distrbituion derived from MaxEnt'''

short_subset = '''Specifications for how you want to subset your data before the
analysis. Note that only the subsetted data will be included in the analysis.
The left-hand dropdown box contains all the columns of your dataset and you may
choose one or more to subset. Please see analysis explanation for more detail
and examples.'''

short_criteria = '''Specifications for how you want to divide your data during
the analysis. The words you see below are the shared columns of your
dataset(s).  You must designate your species column with the special word
'species' found in the dropdown menu. You are not required to fill any
additional columns for this analysis. Please see analysis explanation for more
detail and examples.'''

optional = ''' Optional parameter. Default value:  '''

req = '''Required parameter. '''

def append_value(params, value, index=None):
    '''
    Append the string 'Required' to values in dictionary

    Parameters
    ----------
    params : dict

    value : string

    index : int

    Returns
    -------
    : dict
    '''

    for key in params.iterkeys():
        if index == None:
            params[key] += value
        else:
            params[key][index] += value 
    return params


