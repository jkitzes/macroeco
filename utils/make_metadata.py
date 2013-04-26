#!/usr/bin/env python

'''
Makes minimal metadata for the user
'''

import metadata_writer
import sys


#The user may want to make minimal metadata for multiple files
if len(sys.argv) == 1:
    print "No data files included.  Minimal metadata not made"
else:
    for i in xrange(len(sys.argv)):
        if i > 0:
            metawriter = metadata_writer.MetaWriter(sys.argv[i])
            traitlist = []
            typelist = []
            print "Examining file '" + metawriter.filename + "'..."
            for name in metawriter.column_names:
                cat = raw_input("Is column name '" + name +\
                                    "' categorical? ")
                if cat == "No" or cat == "no" or cat == "n" or\
                   cat == "N":
                    types = (name, {'cat' : False})
                    typelist.append(types)
                    spatial = raw_input("Is column name '" + name +\
                                    "' spatially explicit? ")
                    if spatial == "Yes" or spatial == "yes" or spatial == 'Y'\
                            or spatial == 'y':
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
                        traits = (name, {'precision': str(precision),
                                     'minimum' : str(minimum),\
                                     'maximum' : str(maximum)})
                        traitlist.append(traits)

                else:
                    types = (name, {'cat' : True})
                    typelist.append(types)

            metawriter.add_attribute_types(typelist)
            metawriter.add_attribute_traits(traitlist)
            metawriter.write_meta_data()
                    
                        
                        
                
                   
                   
                
                



