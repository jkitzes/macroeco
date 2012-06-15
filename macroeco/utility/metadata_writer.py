#!/usr/bin/env python

'''
This module contains a minimal metadata writer class for quickly making
metadata

'''


import xml.etree.ElementTree as ET
from matplotlib.mlab import csv2rec

__author__ = "Mark Wilber"
__copyright__ = "Copyright 2012, Regents of University of California"
__credits__ = "John Harte"
__license__ = None
__version__ = "0.1"
__maintainer__ = "Mark Wilber"
__email__ = "mqw@berkeley.edu"
__status__ = "Development"

class MetaWriter:
    '''
    Writes a metadata file based on the given filename and user input
    
    '''

    def __init__(self, datapath):
        '''
        Class takes in a datafile path name and creates an xml tree using the
        column heading of the recarray generated from the csv file.

        Parameters
        ----------
        datapath : string
            Datafile name

        '''
        assert datapath[-4:] == '.csv', "%s must end in .csv" % (datapath)
        self.filename = datapath.split('.')[0]
        self.datafile = csv2rec(datapath)
        self.root = ET.Element('dataset')
        self.dataTable = ET.SubElement(self.root, 'dataTable')
        self.attributeList = ET.SubElement(self.dataTable, 'attributeList')
        self.attributes = []
        for i, name in enumerate(self.datafile.dtype.names):
            attribute = ET.SubElement(self.attributeList, 'attribute')
            attributeName = ET.SubElement(attribute, 'attributeName')
            attributeName.text = name
            self.attributes.append(attribute)

    def add_attribute_traits(self, traitlist):
        '''
        Adds traits to the attributes contained in self.attributes as specified
        by the traitlist.  Traitlist is a list of tuples with each tuple
        containting two elements: the attribute name (string) and a dictionary
        of traits to be added to the attribute.

        Example of traitlist:

            [('x', {'minimum' : '0', 'maximum' : '100'}), ('y', {'precision' :
            '0.1'})]
        
        '''

        for item in traitlist:
            for attribute in self.attributes:
                tree = ET.ElementTree(attribute)
                if (tree.findall('./')[0].text == item[0]):
                    for key in item[1].iterkeys():
                        trait = ET.SubElement(attribute, key)
                        trait.text = item[1][key]

    def write_meta_data(self):
        '''
        Writes out the xml tree that is contained in self.root and saves and
        .xml file in the currect working directory under the given filename 

        
        '''

        tree = ET.ElementTree(self.root)
        tree.write(self.filename + '.xml')






