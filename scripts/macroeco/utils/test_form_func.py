#!/usr/bin/python
#Testing form_func.py

import unittest
from form_func import *
import numpy as np
import os
from matplotlib.mlab import csv2rec
gcwd = os.getcwd #get current directory
pd = os.path.dirname #get parent directory
chdir = os.chdir #change directories
jp = os.path.join 

class TestFormFunc(unittest.TestCase):
    '''Tests the functions with in form_func.py'''

    def setUp(self):
        self.spp_array1 = np.array(['as', 'as', 'as', 'as', 'as'])
        self.spp_array2 = np.array([2,2,3,5,6,3,4,5,7,8])
        self.spp_array3 = np.array(['as','bn', 'as', 'ty', 'bn'])
        self.spp_array4 = np.array([])

        self.arch1 = open('arch1.csv', 'w')
        self.arch1.write('''cell, row, column, AGR, THY, FTW, REW
                        1, 1, 1, 0, 1, 1, 0
                        2, 1, 2, 3, 3, 0, 1
                        3, 2, 1, 0, 0, 0, 0
                        4, 2, 2, 1, 5, 1, 0''')
        self.arch1.close()
        self.arch2 = open('arch2.csv', 'w')
        self.arch2.write('''cell, row, column, AGR, THY, FTW, REW
                        1, 1, 1, 0, 1, 1, 0
                        2, 1, 2, 3, 3, 0, 1
                        3, 2, 1, 0, 0, 0, 0
                        4, 2, 2, 1, 5, 1, 0''')
        self.arch2.close()

    def tearDown(self):
        os.remove('arch1.csv')
        os.remove('arch2.csv')        

    def test_create_intcodes(self):
        unq_specs = np.unique(self.spp_array1)
        unq_ints = np.linspace(0, len(np.unique(self.spp_array1)) - 1,\
                               num=len(np.unique(self.spp_array1)))
        tot_int = create_intcodes(self.spp_array1, unq_specs, unq_ints)
        self.assertTrue(len(tot_int) == 5)
        self.assertTrue(np.unique(tot_int)[0] == .0)
        self.assertTrue(np.all(np.equal(tot_int, np.array([.0,.0,.0,.0,.0]))))
        unq_specs = np.unique(self.spp_array2)
        unq_ints = np.linspace(0, len(np.unique(self.spp_array2)) - 1, \
                            num=len(np.unique(self.spp_array2)))
        tot_int = create_intcodes(self.spp_array2, unq_specs, unq_ints)
        self.assertTrue(len(tot_int) == len(self.spp_array2))
        self.assertTrue(np.all(np.equal(np.unique(tot_int), 
                                                np.linspace(0,6,num=7))))
        self.assertRaises(AssertionError, create_intcodes, self.spp_array4, 
                                unq_specs, unq_ints)

    def test_add_field(self):
        data = csv2rec('arch1.csv')
        data_added = add_field(data, [('test', np.int)])
        names = np.array(data_added.dtype.names)
        self.assertTrue(sum(names == 'test') == 1)

    def test_merge_formatted(self):
        data1 = csv2rec('arch1.csv')
        data2 = csv2rec('arch2.csv')
        dl = [data1, data2]
        merged = merge_formatted(dl)
        self.assertTrue(sum(merged['rew']) == 2)
        self.assertTrue(sum(merged['column']) == 12)

    def test_add_data_fields(self):
        data1 = csv2rec('arch1.csv')
        data2 = csv2rec('arch2.csv')
        dl = [data1, data2]
        alt_data = add_data_fields(dl, ['year'], [(1998, 2002)])
        self.assertTrue(np.all(alt_data[0]['year'] == '1998'))
        self.assertTrue(np.all(alt_data[1]['year'] == '2002'))
        alt_data = add_data_fields(dl, ['year', 'why'], [(1998, 2002),
                                   ('h','a')])
        self.assertTrue(np.all(alt_data[0]['why'] == 'h'))

    def test_fractionate(self):
        data1 = csv2rec('arch1.csv')
        data2 = csv2rec('arch2.csv')
        dl = [data1, data2]
        fr = fractionate(dl, (10, 10), (5, 5), ['row', 'column'])
        self.assertTrue(fr[0]['row'][3] == 5)
        self.assertTrue(fr[1]['column'][2] == 0)

    def test_format_dense(self):
        data1 = csv2rec('arch1.csv')
        data2 = csv2rec('arch2.csv')
        dl = [data1, data2]
        form = format_dense(dl, 3, (4,4))
        self.assertTrue(np.all(form[0]['count'][:4] == np.array([1,1,3,3])))
        self.assertTrue(np.all(form[1]['count'] ==
                                               np.array([1,1,3,3,1,1,5,1])))



        




