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
        self.xyfile1.write('''spp_code, x, y
                       0, 0, 0
                       0, 0, 0
                       0, 0, 1
                       1, 0, 0
                       1, 1, 0''')
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

#    def test_accurate_sub_tables(self):
#        xy2 = DataTable('xyfile2.csv')  # More complex table
#
#        sub1 = xy2.get_sub_table(0, 1, 0, 1)  # One cell
#        np.testing.assert_array_equal(sub1, self.xyarr2[[0, 1, 3]])
#        sub2 = xy2.get_sub_table(0, 1, 0, 2)  # Half plot
#        np.testing.assert_array_equal(sub2, self.xyarr2[[0, 1, 2, 3]])
#        sub3 = xy2.get_sub_table(0, 2, 0, 2)  # Whole plot
#        np.testing.assert_array_equal(sub3, self.xyarr2)
        

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

<dataTable><attributeList>

<attribute><attributeName>x</attributeName>
<precision>0.1</precision>
<numericDomain><numberType>real</numberType>
<bounds><minimum exclusive="false">0</minimum>
<maximum exclusive="false">99.9</maximum></bounds>
<type>metric</type>
</numericDomain>
</attribute>

</attributeList></dataTable>
</dataset></eml:eml>''')
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
                                    ('x', 'minimum'): 0,
                                    ('x', 'precision'): 0.1,
                                    ('x', 'type'): 'metric',
                                    ('y', 'maximum'): None,
                                    ('y', 'minimum'): None,
                                    ('y', 'precision'): None,
                                    ('y', 'type'): None})

    def test_physical_coverage(self):
        meta = Metadata('xyfile1.csv', [])
        edges = meta.get_physical_coverage()
        self.assertEqual(edges, [8.975, -79.5915, 10, -79.5915])

    def test_title(self):
        meta = Metadata('xyfile1.csv', [])
        self.assertEqual(meta.get_title(), 'Unittest XML')


if __name__ == '__main__':
    unittest.main()
