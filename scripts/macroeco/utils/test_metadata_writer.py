#!/usr/bin/python
#Testing metadata_writer.py

import unittest
from metadata_writer import *
import numpy as np
import xml.etree.ElementTree as ET

import os
gcwd = os.getcwd #get current directory
pd = os.path.dirname #get parent directory
chdir = os.chdir #change directories
jp = os.path.join

class TestMetadataWriter(unittest.TestCase):
    '''Tests the MetaWriter class'''

    def setUp(self):
        
        self.meta1 = open('meta1.csv', 'w')
        self.meta1.write('''cell,row,column,spp,year
                        0, 1, 2, 3, 4''')
        self.meta1.close()

    def tearDown(self):
        os.remove('meta1.csv')

    def test_metawriter(self):
        mt = MetaWriter('meta1.csv')
        att = [('row', {'maximum' : '2', 'minimum' : '0', 'precision' :\
                        '0.1'}), ('column', {'maximum' : '45', 'minimum' :\
                        '0', 'precision' : '1'})]

        # Check that all attributes are in xml tree
        self.assertTrue(len(mt.attributeList.findall('./')) == 5)
        
        # Check that all types are ordinal by default
        measure = mt.attributeList.findall('./attribute/measurementScale')
        for i, m in enumerate(measure):
            temp = m.findall('./')
            self.assertTrue(temp[0].tag == 'ordinal')

        # Check that it adds correct attribute types
        types = [('cell', {'cat' : True}), ('row', {'cat' : False}), ('column', {'cat'
        : False}), ('spp' , {'cat' : True})]

        mt.add_attribute_types(types)
        order = ['ordinal', 'interval', 'interval', 'ordinal', 'ordinal']
        for i, att in enumerate(mt.attributes):
            temp = att.findall('./measurementScale/' + order[i])
            self.assertTrue(len(temp) == 1)

        # Check that it overwrites types if they are changed
        types = [('cell', {'cat' : False}), ('row', {'cat' : True}), ('column', {'cat'
        : True}), ('spp' , {'cat' : False})]

        mt.add_attribute_types(types)

        mt.add_attribute_types(types)
        order = ['interval', 'ordinal', 'ordinal', 'interval', 'ordinal']
        for i, att in enumerate(mt.attributes):
            temp = att.findall('./measurementScale/' + order[i])
            self.assertTrue(len(temp) == 1)

        # Check that max, min and precision are set correctly

        types = [('cell', {'cat' : True}), ('row', {'cat' : False}), ('column', {'cat'
        : False}), ('spp' , {'cat' : True})]

        mt.add_attribute_types(types)

        att_list = [('row', {'minimum' : 0, 'maximum' : 400, 'precision' : 3,
        'random' : 'harry'}), ('column', {'maximum' : 5}), ('spp', {'precision'
        : 4})]

        mt.add_attribute_traits(att_list)

        # spp should have no precision even though we tried to add it
        have = mt.attributeList.findall(".//attributeName")
        names = [nm.text for nm in have]
        ind = names.index('spp')
        maybe =\
            mt.attributes[ind].findall('./measurementScale/ordinal/precision')
        self.assertTrue(len(maybe) == 0)
        
        # cell should have no precision
        have = mt.attributeList.findall(".//attributeName")
        names = [nm.text for nm in have]
        ind = names.index('cell')
        maybe =\
        mt.attributes[ind].findall('./measurementScale/ordinal/precision')
        self.assertTrue(len(maybe) == 0)

        # Precision of row should be three 
        have = mt.attributeList.findall(".//attributeName")
        names = [nm.text for nm in have]
        ind = names.index('row')
        maybe =\
        mt.attributes[ind].findall('./measurementScale/interval/precision')
        self.assertTrue(maybe[0].text == "3")
        
        # Precision of column should be 0
        have = mt.attributeList.findall(".//attributeName")
        names = [nm.text for nm in have]
        ind = names.index('column')
        maybe =\
        mt.attributes[ind].findall('./measurementScale/interval/precision')
        self.assertTrue(maybe[0].text == "0")
        
        # Maximum is set right
        have = mt.attributeList.findall(".//attributeName")
        names = [nm.text for nm in have]
        ind = names.index('column')
        maybe =\
        mt.attributes[ind].findall('./measurementScale/interval/numericDomain/bounds/maximum')
        self.assertTrue(maybe[0].text == "5")

         

if __name__ == '__main__':
    unittest.main()
