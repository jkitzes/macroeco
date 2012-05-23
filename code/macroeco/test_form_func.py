#!/usr/bin/python
#Testing form_func.py

import unittest
from form_func import *
import numpy as np
import os
import matplotlib.mlab as plt
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
        self.arch2.write('''cell, row, colun, AGR, THY, FTW, REW
                        1, 1, 1, 0, 1, 1, 0
                        2, 1, 2, 3, 3, 0, 1
                        3, 2, 1, 0, 0, 0, 0
                        4, 2, 2, 1, 5, 1, 0''')
        self.arch2.close()





    def tearDown(self):
        os.remove('arch1.csv')
        os.remove('arch2.csv')



    def test_make_spec_dict(self):
        self.assertRaises(AssertionError, make_spec_dict, self.spp_array4)
        spec_dict = make_spec_dict(self.spp_array1)
        self.assertTrue(len(spec_dict) == 1)
        self.assertTrue(spec_dict['spp_code'][0] == 0)
        self.assertTrue(spec_dict['spp'][0] == 'as')
        spec_dict = make_spec_dict(self.spp_array2)
        self.assertTrue(len(spec_dict) == 7)
        self.assertTrue(sum(np.equal(spec_dict['spp_code'], \
                            np.array([0,1,2,3,4,5,6]))) == 7)
        spec_dict = make_spec_dict(self.spp_array3)
        self.assertTrue(len(spec_dict) == 3)
        

    def test_create_intcodes(self):
        unq_specs = np.unique(self.spp_array1)
        unq_ints = np.linspace(0, len(np.unique(self.spp_array1)) - 1, num=len(\
                              np.unique(self.spp_array1)))
        tot_int = create_intcodes(self.spp_array1, unq_specs, unq_ints)
        self.assertTrue(len(tot_int) == 5)
        self.assertTrue(np.unique(tot_int)[0] == .0)
        self.assertTrue(np.all(np.equal(tot_int, np.array([.0,.0,.0,.0,.0]))))
        unq_specs = np.unique(self.spp_array2)
        unq_ints = np.linspace(0, len(np.unique(self.spp_array2)) - 1, num=len(\
                              np.unique(self.spp_array2)))
        tot_int = create_intcodes(self.spp_array2, unq_specs, unq_ints)
        self.assertTrue(len(tot_int) == len(self.spp_array2))
        self.assertTrue(np.all(np.equal(np.unique(tot_int), np.linspace(0,6,num=7))))
        self.assertRaises(AssertionError, create_intcodes, self.spp_array4, unq_specs,\
                                                          unq_ints)

    def test_make_dense_spec_dict(self):
        datayears = [plt.csv2rec('arch1.csv')]
        spec_dict = make_dense_spec_dict(datayears)
        self.assertTrue(len(spec_dict) == 4)
        self.assertTrue(set(spec_dict['spp']) == set(['agr', 'thy', 'ftw', 'rew'])) 
        self.assertTrue(set(spec_dict['spp_code']) == set([0,1,2,3]))
    

if __name__ == '__main__':
    unittest.main()

