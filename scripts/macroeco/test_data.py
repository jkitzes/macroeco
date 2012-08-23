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
        ''' Write test xytable csv files with and without count. '''

        self.xyfile1 = open('xyfile1.csv','w')
        self.xyfile1.write('''spp_code, x, y
                       0, 0, 0
                       0, 0, 0
                       0, 0, 1
                       1, 0, 0
                       1, 1, 0''')
        self.xyfile1.close()
        self.xyarr1 = csv2rec('xyfile1.csv')
        self.xyhead1 = ['spp_code', 'x', 'y']

        self.xyfile2 = open('xyfile2.csv','w')
        self.xyfile2.write('''spp_code, x, y, count
                       0, 0, 0, 1
                       0, 0, 0, 1
                       0, 0, 1, 1
                       1, 0, 0, 1
                       1, 1, 0, 1''')
        self.xyfile2.close()
        self.xyarr2 = csv2rec('xyfile2.csv')
        self.xyhead2 = ['spp_code', 'x', 'y', 'count']

        self.xyfile3 = open('xyfile3.csv','w')
        self.xyfile3.write('''bad, x, y, count
                       0, 0, 0, 1
                       0, 0, 0, 1
                       0, 0, 1, 1
                       1, 0, 0, 1
                       1, 1, 0, 1''')
        self.xyfile3.close()

        self.xyfile4 = open('xyfile4.csv','w')
        self.xyfile4.write('''spp_code, x, y, count
                       0, .1, .1, 2
                       0, .1, .2, 1
                       0, .1, .3, 1
                       2, .1, .2, 1
                       2, .2, .3, 1''')
        self.xyfile4.close()

        self.xymeta = {'precision': 1, 'xrange': (0,1), 'yrange': (0,1)}
        self.xymeta4 = {'precision': .1, 'xrange': (.1,.2), 'yrange': (.1,.3)}

    def tearDown(self):
        os.remove('xyfile1.csv')
        os.remove('xyfile2.csv')
        os.remove('xyfile3.csv')
        os.remove('xyfile4.csv')
        return

    #
    # xytable_load
    #

    def test_error_if_file_type_not_csv(self):
        self.assertRaises(AssertionError, DataTable, 'file.txt')
        
    def test_None_if_no_meta_file(self):
        xy1 = DataTable('xyfile1.csv')
        self.assertEqual(xy1.meta, None)

    # TODO: Test writing and reading metadata file

    # TODO: Test throws error if loaded metadata is missing required key

    #
    # sub_table
    #

    def test_accurate_sub_tables(self):
        xy2 = DataTable('xyfile2.csv')  # More complex table

        sub1 = xy2.get_sub_table(0, 1, 0, 1)  # One cell
        np.testing.assert_array_equal(sub1, self.xyarr2[[0, 1, 3]])
        sub2 = xy2.get_sub_table(0, 1, 0, 2)  # Half plot
        np.testing.assert_array_equal(sub2, self.xyarr2[[0, 1, 2, 3]])
        sub3 = xy2.get_sub_table(0, 2, 0, 2)  # Whole plot
        np.testing.assert_array_equal(sub3, self.xyarr2)
        

class TestMetadata(unittest.TestCase):
    # Uses BCIS archival metadata as sample
    # TODO: Path below may not work on all OS
    
    def setUp(self):
        self.meta = Metadata('../../data/archival/BCIS/BCIS_1995.csv')        

    # TODO: Will need to add check for datatype here also
    def test_BCIS_Archival_metadata(self):
        asklist = [('gx','precision'),
                        ('gy','precision'),
                        ('gx','minimum'),
                        ('gx','maximum'),
                        ('gy','maximum')]
        values = ['0.1', '0.1','0','999.9','499.9']
        expect = dict(zip(asklist, values))

        self.meta.get_dataTable_values(asklist)
        for ask in asklist:
            assert self.meta.TableDescr[ask] == expect[ask]

    def test_region(self):
        edges = self.meta.get_coverage_region()
        assert edges == [9.15,-79.85,9.15,-79.85]


if __name__ == '__main__':
    unittest.main()
