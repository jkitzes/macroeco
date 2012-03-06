'''
Unit tests for data.py
'''

import unittest
import os
from data import *

class TestXytable(unittest.TestCase):
    ''' Tests both load and add_count functions. '''

    def setUp(self):
        ''' Write test xytable files with and without count. '''

        self.xyfile1 = open('xyfile1.csv','w')
        self.xyfile1.write('''spp_code, x, y
                       0, 0, 0
                       0, 0, 0
                       0, 1, 1
                       1, 0, 0
                       1, 1, 0''')
        self.xyfile1.close()
        self.xyarr1 = np.loadtxt('xyfile1.csv', delimiter = ',', skiprows = 1)

        self.xyfile2 = open('xyfile2.csv','w')
        self.xyfile2.write('''spp_code, x, y, count
                       0, 0, 0, 1
                       0, 0, 0, 1
                       0, 1, 1, 1
                       1, 0, 0, 1
                       1, 1, 0, 1''')
        self.xyfile2.close()
        self.xyarr2 = np.loadtxt('xyfile2.csv', delimiter = ',', skiprows = 1)

        self.xyfile3 = open('xyfile3.csv','w')
        self.xyfile3.write('''bad, x, y, count
                       0, 0, 0, 1
                       0, 0, 0, 1
                       0, 1, 1, 1
                       1, 0, 0, 1
                       1, 1, 0, 1''')
        self.xyfile3.close()

        self.xyhead1 = ['spp_code', 'x', 'y']
        self.xyhead2 = ['spp_code', 'x', 'y', 'count']
        self.xymeta = {'precision': 1, 'xrange': (0,1), 'yrange': (0,1)}

    def tearDown(self):
        os.remove('xyfile1.csv')
        os.remove('xyfile2.csv')
        os.remove('xyfile3.csv')
        return

    #
    # xytable_load
    #

    def test_error_if_file_type_not_csv(self):
        self.assertRaises(NotImplementedError, xytable_load, 'file.txt', {})

    def test_error_if_missing_required_header(self):
        self.assertRaises(NameError, xytable_load, 'xyfile3.csv', self.xymeta)  
        
    def test_return_input_meta_if_given(self):
        data, head, meta = xytable_load('xyfile1.csv', self.xymeta)
        self.assertEqual(self.xymeta, meta)
    
    def test_error_if_input_meta_empty_or_not_given(self):
        self.assertRaises(NotImplementedError, xytable_load, 'xyfile1.csv', {})
        self.assertRaises(NotImplementedError, xytable_load, 'xyfile1.csv')

    def test_error_if_meta_missing_required_key(self):
        self.assertRaises(NameError, xytable_load, 'xyfile1.csv', {'a': 'a'})

    def test_returns_correct_data_and_meta(self):
        data1, head, meta = xytable_load('xyfile1.csv', self.xymeta)
        data2, head, meta = xytable_load('xyfile2.csv', self.xymeta)

    #
    # xytable_add_count
    #

    def test_return_orig_data_if_has_count(self):
        data, head = xytable_add_count(self.xyarr2, self.xyhead2, compress = 
                                       False)
        np.testing.assert_array_equal(data, self.xyarr2)
        self.assertEquals(head, self.xyhead2)

    def test_return_appended_data_head(self):
        data, head = xytable_add_count(self.xyarr1, self.xyhead1, compress = 
                                       False)
        np.testing.assert_array_equal(data, self.xyarr2)
        self.assertEquals(head, self.xyhead2)

    def test_error_if_compress_true(self):
        self.assertRaises(NotImplementedError, xytable_add_count, self.xyarr1, 
                          self.xyhead1, True)


if __name__ == '__main__':
    unittest.main()
