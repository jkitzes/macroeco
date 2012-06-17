#!/usr/bin/env python

'''
Makes minimal metadata for the user
'''

from macroeco.utility import metadata_writer
import sys


__author__ = "Mark Wilber"
__copyright__ = "Copyright 2012, Regents of University of California"
__credits__ = "John Harte"
__license__ = None
__version__ = "0.1"
__maintainer__ = "Mark Wilber"
__email__ = "mqw@berkeley.edu"
__status__ = "Development"

#The user may want to make minimal metadata for multiple files
if len(sys.argv) == 1:
    print "No data files included.  Minimal metadata not made"
else:
    for i in xrange(len(sys.argv)):
        if i > 0:
            metawriter = metadata_writer.MetaWriter(sys.argv[i])
            traitlist = []
            print "Examining file '" + metawriter.filename + "'..."
            for name in metawriter.datafile.dtype.names:
                spatial = raw_input("Is column name '" + name +\
                                    "' spatially explicit? ")
                if spatial == "Yes" or spatial == "yes" or spatial == "y" or\
                   spatial == "Y":
                    while True:
                        minimum = raw_input("Please enter the minimum value" +\
                                            " of column '" + name + "': ")
                        maximum = raw_input("Please enter the maximum value" +\
                                            " of column '" + name + "': ")
                        precision = raw_input("Please enter the precision" +\
                                            " of column '" + name + "': ")
                        try:
                            minimum = float(minimum)
                            maximum = float(maximum)
                            precision = float(precision)
                            break #This might not work
                        except ValueError:
                            print "Maximum, minimum, and precision must all" +\
                                  " be real numbers"
                    traits = (name, {'minimum' : str(minimum),\
                                     'maximum' : str(maximum),\
                                     'precision' : str(precision)})
                    traitlist.append(traits)
            metawriter.add_attribute_traits(traitlist)
            metawriter.write_meta_data()
                    
                        
                        
                
                   
                   
                
                



