#!/usr/bin/python

'''
Unit tests for cmn_analysis.py

'''

from __future__ import division
import unittest
import os
gcwd = os.getcwd
pd = os.path.dirname
jp = os.path.join
from empirical import *
from cmn_analysis import *
import numpy as np
import random


class TestPatch(unittest.TestCase):

    def setUp(self): 

        self.xyfile6 = open('xyfile6.csv', 'w')
        self.xyfile6.write('''spp_code, x, y, count
                        0, 0, 0, 1
                        1, 0, 0, 1
                        2, 0, 0, 0
                        3, 0, 0, 3
                        0, 0, 1, 0
                        1, 0, 1, 4
                        2, 0, 1, 0
                        3, 0, 1, 1
                        0, 1, 0, 1
                        1, 1, 0, 0
                        2, 1, 0, 3
                        3, 1, 0, 1
                        0, 1, 1, 0
                        1, 1, 1, 1
                        2, 1, 1, 3
                        3, 1, 1, 1''')
        self.xyfile6.close()
        self.xymeta6 = {'precision': 1, 'xrange':(0,1), 'yrange':(0,1)}
        self.gridtest = Patch('xyfile6.csv')
        self.gridtest.xy.meta = self.xymeta6
        self.gridtest.set_attributes()

    def tearDown(self):
        os.remove('xyfile6.csv')

    def test_get_common_arrays(self):
        common_arrays = get_common_arrays(self.gridtest, [(2,2)])
        common = self.gridtest.QS_grid([(2,2)])
        self.assertTrue(np.array_equal(common_arrays[0]['dist'], np.unique(common[0][:,2])))
        self.assertTrue(np.array_equal(common_arrays[0]['cmn'], np.array([11/15, 8/15])))
        common_arrays = get_common_arrays(self.gridtest, [(1,1)])

if __name__ == '__main__':
    unittest.main()

