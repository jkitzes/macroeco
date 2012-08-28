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
        self.meta1.write('''cell, row, column, spp, year
                        0, 1, 2, 3, 4''')
        self.meta1.close()

    def tearDown(self):
        os.remove('meta1.csv')

    def test_metawriter(self):
        '''
        '''
        mt = MetaWriter('meta1.csv')
        att = [('row', {'maximum' : '2', 'minimum' : '0', 'precision' :\
                        '0.1'}), ('column', {'maximum' : '45', 'minimum' :\
                        '0', 'precision' : '1'})]
        mt.add_attribute_traits(att)
        self.assertTrue(len(mt.attributeList.findall('./')) == 5)
        for child in mt.attributeList:
           self.assertTrue(len(child.findall('./')) == 2) 
        att = [('year', {'hello' : 5})]
        mt.add_attribute_traits(att)
        ET.dump(mt.root)
        types = [('cell' , {'cat' : False}), ('year', {'cat': True}), \
                 ('spp', {'cat' : False}), ('row', {'cat' : False}), \
                 ('column', {'cat' : False})]
        mt.add_attribute_types(types)
        for child in mt.attributeList:
            pres = False
            for item in child.findall('./'):
                if item.tag == 'ordinal' or item.tag == 'interval':
                    pres = True
            self.assertTrue(pres)

        types = [('year', {'cat' : False})]
        mt.add_attribute_types(types)
        for child in mt.attributeList:
            attnm = child.findall('attributeName')[0]
            if attnm.text == 'year':
                self.assertTrue(len(child.findall('interval')) == 1)
        types = [('year', {'ct' : False})]
        self.assertRaises(KeyError, mt.add_attribute_types, types)

if __name__ == '__main__':
        unittest.main()



