'''
Unit tests for data.py
'''

import unittest
import os
import numpy as np
from matplotlib.mlab import csv2rec
from macroeco.data import DataTable, Metadata

class TestDataTable(unittest.TestCase):

    def setUp(self):
        '''Write test xytable csv file.'''

        self.xyfile1 = open('xyfile1.csv','w')
        self.xyfile1.write('''spp_code, x, y, count
                       0, 0, 0, 1
                       0, 0, 0, 2
                       0, 0, 1, 1
                       1, 0, 0, 1
                       1, 1, 0, 2''')
        self.xyfile1.close()
        self.xyarr1 = csv2rec('xyfile1.csv')

    def tearDown(self):
        os.remove('xyfile1.csv')

    def test_error_if_file_type_not_csv(self):
        self.assertRaises(AssertionError, DataTable, 'file.txt')
        
    def test_meta_None_if_no_meta_file(self):
        xy1 = DataTable('xyfile1.csv')
        self.assertEqual(xy1.meta, None)

    def test_table_is_correct(self):
        xy1 = DataTable('xyfile1.csv')
        np.testing.assert_array_equal(xy1.table, self.xyarr1)

    def test_get_subtable(self):
        xy1 = DataTable('xyfile1.csv')
        xy1.meta = {('x', 'maximum'): 1,
                    ('x', 'minimum'): 0,
                    ('x', 'precision'): 1,
                    ('y', 'maximum'): 1,
                    ('y', 'minimum'): 0,
                    ('y', 'precision'): 1}

        # Whole table
        sub = xy1.get_subtable({})
        np.testing.assert_array_equal(sub, self.xyarr1)

        sub = xy1.get_subtable({'x': ('>=0','<2'), 'y': ('>=0','<2')})
        np.testing.assert_array_equal(sub, self.xyarr1)

        # Subset
        sub = xy1.get_subtable({'spp_code': '==0'})
        np.testing.assert_array_equal(sub, self.xyarr1[0:3])

        sub = xy1.get_subtable({'spp_code': '==0', 'x': '>0'})
        np.testing.assert_array_equal(sub, self.xyarr1[2])

class TestMetadata(unittest.TestCase):
    
    def setUp(self):
        '''Write test data and metadata file.'''

        self.xyfile1 = open('xyfile1.csv','w')
        self.xyfile1.write('''x, y
                       0, 0
                       0, 0
                       0, 0
                       1, 0
                       1, 1''')
        self.xyfile1.close()

        self.xymeta = open('xyfile1.xml','w')
        self.xymeta.write('''<?xml version="1.0" encoding="UTF-8"?>
<eml:eml packageId="macroeco.129.1" system="knb"                          
xmlns:eml="eml://ecoinformatics.org/eml-2.1.0" 
xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" 
xsi:schemaLocation="eml://ecoinformatics.org/eml-2.1.0 
eml.xsd"><dataset>

<title>Unittest XML</title> 

<coverage><geographicCoverage><geographicDescription>NA</geographicDescription>
<boundingCoordinates><westBoundingCoordinate>-79.5915</westBoundingCoordinate>
<eastBoundingCoordinate>-79.5915</eastBoundingCoordinate>
<northBoundingCoordinate>8.975</northBoundingCoordinate>
<southBoundingCoordinate>10</southBoundingCoordinate>
</boundingCoordinates>
</geographicCoverage></coverage>

<dataTable><attributeList><attribute>
<attributeName>y</attributeName><ordinal /></attribute><attribute>
<attributeName>cell</attributeName><interval /></attribute><attribute>
<attributeName>x</attributeName><interval><minimum>0.0</minimum>
<maximum>99.9</maximum><precision>0.1</precision></interval></attribute>
</attributeList></dataTable></dataset></eml:eml>''')
        self.xymeta.close()

    def tearDown(self):
        os.remove('xyfile1.csv')
        os.remove('xyfile1.xml')

    def test_metadata_correct_read(self):
        # Should read values correctly from sample file, including None for
        # attributes that do not exist and elements that do not exist.
        xy1 = DataTable('xyfile1.csv')
        self.assertEqual(len(xy1.meta), 8)
        self.assertEqual(xy1.meta, {('x', 'maximum'): 99.9,
                                    ('x', 'minimum'): 0.0,
                                    ('x', 'precision'): 0.1,
                                    ('x', 'type'): 'interval',
                                    ('y', 'maximum'): None,
                                    ('y', 'minimum'): None,
                                    ('y', 'precision'): None,
                                    ('y', 'type'): 'ordinal'})

    def test_physical_coverage(self):
        meta = Metadata('xyfile1.csv', [])
        edges = meta.get_physical_coverage()
        self.assertEqual(edges, [8.975, -79.5915, 10, -79.5915])

    def test_title(self):
        meta = Metadata('xyfile1.csv', [])
        self.assertEqual(meta.get_title(), 'Unittest XML')
