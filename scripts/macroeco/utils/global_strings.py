#!/usr/bin/python

'''This python file contains global strings used in the scripts.  Consolidated
in this script for easy maintenance'''

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
