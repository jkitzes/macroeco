'''
Unit tests for data.py
'''

import unittest
import os
from data import *

class TestXYTable(unittest.TestCase):

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
        self.xyarr1 = np.loadtxt('xyfile1.csv', delimiter = ',', skiprows = 1)
        self.xyhead1 = ['spp_code', 'x', 'y']

        self.xyfile2 = open('xyfile2.csv','w')
        self.xyfile2.write('''spp_code, x, y, count
                       0, 0, 0, 1
                       0, 0, 0, 1
                       0, 0, 1, 1
                       1, 0, 0, 1
                       1, 1, 0, 1''')
        self.xyfile2.close()
        self.xyarr2 = np.loadtxt('xyfile2.csv', delimiter = ',', skiprows = 1)
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
        self.xymeta_none = {'precision': 0, 'xrange': (0, 0), 
                            'yrange': (0, 0)}

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
        self.assertRaises(NotImplementedError, XYTable, 'file.txt', False)

    def test_error_if_missing_required_header(self):
        self.assertRaises(NameError, XYTable, 'xyfile3.csv', False)  
        
    def test_return_meta_dict_with_0_if_no_meta_file(self):
        xy = XYTable('xyfile1.csv', False)
        self.assertEqual(self.xymeta_none, xy.meta)

    # TODO: Test writing and reading metadata file

    # TODO: Test throws error if loaded metadata is missing required key

    def test_attribute_header_column_indexes_correct(self):
        xy4 = XYTable('xyfile4.csv', False)  # More complex table
        self.assertEqual(xy4.col_spp_code, 0)
        self.assertEqual(xy4.col_count, 3)

    def test_attribute_spp_codes_correct(self):
        xy4 = XYTable('xyfile4.csv', False)  # More complex table
        np.testing.assert_array_equal(xy4.spp_codes, np.array((0, 2)))
        self.assertEqual(xy4.max_spp_code, 2)

    #
    # xytable_add_count
    #

    # This test effectively tests the entire xytable class as well
    def test_adding_count_works_properly(self):
        xy1 = XYTable('xyfile1.csv', False)  # Simple table without count
        xy2 = XYTable('xyfile2.csv', False)  # Simple table with count

        np.testing.assert_array_equal(xy1.table, self.xyarr2)
        self.assertEqual(xy1.head, self.xyhead2)
        self.assertEqual(xy1.meta, self.xymeta_none)

        np.testing.assert_array_equal(xy2.table, self.xyarr2)
        self.assertEqual(xy2.head, self.xyhead2)
        self.assertEqual(xy2.meta, self.xymeta_none)

    def test_error_if_compress_true(self):
        self.assertRaises(NotImplementedError, XYTable, 'xyfile1.csv', True)

    #
    # sub_table
    #

    def test_accurate_sub_tables(self):
        xy2 = XYTable('xyfile2.csv', False)  # More complex table

        sub1 = xy2.get_sub_table(0, 1, 0, 1)  # One cell
        np.testing.assert_array_equal(sub1, self.xyarr2[[0, 1, 3]])
        sub2 = xy2.get_sub_table(0, 1, 0, 2)  # Half plot
        np.testing.assert_array_equal(sub2, self.xyarr2[[0, 1, 2, 3]])
        sub3 = xy2.get_sub_table(0, 2, 0, 2)  # Whole plot
        np.testing.assert_array_equal(sub3, self.xyarr2)
        

class MorphoMetadata(unittest.TestCase):
    def setUp(self):
        self.meta = metadata.Metadata('../../data/archival/BCIS/BCIS_1995.csv')        

    def test_tablefromBCIS_Archival(self):
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
