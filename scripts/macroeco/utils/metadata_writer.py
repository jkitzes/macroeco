#!/usr/bin/env python

'''
This module contains a minimal metadata writer class for quickly making
metadata

'''


import xml.etree.ElementTree as ET
import os


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
        fin = open(datapath, 'r')
        self.column_names = fin.readline().strip().split(',')
        fin.close()
        self.root = ET.Element('eml:eml')
        self.root.attrib = {'packageID' : self.filename, 'system' : 'knb',
            "xmlns:eml" : "eml://ecoinformatics.org/eml-2.1.0", 'xmlns:xsi': 
            "http://www.w3.org/2001/XMLSchema-instance", "xsi:schemaLocation"
            : "eml://ecoinformatics.org/eml-2.1.0 eml.xsd"}
        self.dataset = ET.SubElement(self.root, 'dataset')
        self.dataTable = ET.SubElement(self.dataset, 'dataTable')
        self.entityName = ET.SubElement(self.dataTable, 'entityName')
        self.entityName.text = os.path.split(datapath)[1]
        self.attributeList = ET.SubElement(self.dataTable, 'attributeList')
        self.attributes = []
        self.attributeTypes = []
        for i, name in enumerate(self.column_names):
            attribute = ET.SubElement(self.attributeList, 'attribute')
            attributeName = ET.SubElement(attribute, 'attributeName')
            attributeType = ET.SubElement(attribute, 'unknown')
            attributeName.text = name
            self.attributes.append(attribute)
            self.attributeTypes.append(attributeType)

    def add_attribute_traits(self, traitlist):
        '''
        Adds traits to the attributes contained in self.attributes as specified
        by the traitlist.  Traitlist is a list of tuples with each tuple
        containting two elements: the attribute name (string) and a dictionary
        of traits to be added to the attribute.
        
        Parameters
        ----------
        traitlist : list
            A list of 2 element tuples where the first element contains a
            string and the second element conatins a dict. See example in
            docstring.

        Example of traitlist:

            [('x', {'minimum' : '0', 'maximum' : '100'}), ('y', {'precision' :
            '0.1'})]
        
        '''

        for item in traitlist:
            for attribute in self.attributes:
                tree = ET.ElementTree(attribute)
                for child in tree.findall('attributeName'):
                    if child.text == item[0]:
                        #TODO:Cleaner way to do this than with if?
                        if len(tree.findall('unknown')) == 1:
                            for key in item[1].iterkeys():
                                att_type = tree.findall('unknown')[0]
                                trait = ET.SubElement(att_type, key)
                                trait.text = str(item[1][key])
                        if len(tree.findall('ordinal')) == 1:
                            for key in item[1].iterkeys():
                                att_type = tree.findall('ordinal')[0]
                                trait = ET.SubElement(att_type, key)
                                trait.text = str(item[1][key])
                        if len(tree.findall('interval')) == 1:
                            for key in item[1].iterkeys():
                                att_type = tree.findall('interval')[0]
                                trait = ET.SubElement(att_type, key)
                                trait.text = str(item[1][key])

    def add_attribute_types(self, typelist):
        '''
        Sets the type of the attribute to either ordinal (categorical) or
        interval (categorical). Initialized in constructor as unknown.

        Parameters
        ----------
        typelist : list
            A list of tuples.  Each tuple contains 2 elements: a string and a
            dict. The dict must contain the keywork cat (categorical) or a 
            KeyError will be thrown.

        Example of typelist:

            [('x', {'cat' : True}), ('y' : {'cat' : True}), ('year',
            {'cat' : False}]

        '''

        for item in typelist:
            for attribute in self.attributes:
                tree = ET.ElementTree(attribute)
                att = tree.findall('attributeName')[0]
                if (att.text == item[0]):
                    if item[1]['cat'] == True:
                        if len(tree.findall('unknown')) == 1:
                            att_type = tree.findall('unknown')[0]
                            att_type.tag = 'ordinal'
                        elif len(tree.findall('interval')) == 1:
                            att_type = tree.findall('interval')[0]
                            att_type.tag = 'ordinal'
                    elif item[1]['cat'] == False:
                        if len(tree.findall('unknown')) == 1:
                            att_type = tree.findall('unknown')[0]
                            att_type.tag = 'interval'
                        elif len(tree.findall('ordinal')) == 1:
                            att_type = tree.findall('ordinal')[0]
                            att_type.tag = 'interval'

    def write_meta_data(self):
        '''
        Writes out the xml tree that is contained in self.root and saves and
        .xml file in the currect working directory under the given filename 

        
        '''

        tree = ET.ElementTree(self.root)
        tree.write(self.filename + '.xml')






