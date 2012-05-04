'''
Unit tests for empirical.py
'''

from __future__ import division
import unittest
import os
from empirical import *
import numpy as np


class TestPatch(unittest.TestCase):

    def setUp(self):
        self.xyfile5 = open('xyfile5.csv','w')
        self.xyfile5.write('''spp_code, x, y, count
                       0, .1, .1, 2
                       0, .1, .2, 1
                       0, .1, .3, 1
                       2, .1, .2, 1
                       2, .2, .3, 1''')
        self.xyfile5.close()
        self.xymeta5 = {'precision': .1, 'xrange': (.1,.2), 'yrange': (.1,.3)}

        self.pat = Patch('xyfile5.csv')
        self.pat.xy.meta = self.xymeta5  # TODO: Replace with reader 
        self.pat.set_attributes()

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

        self.xyfile7 = open('xyfile7.csv', 'w')
        self.xyfile7.write('''spp_code, x, y, count
                        0, 1, 1, 1
                        1, 1, 1, 1
                        2, 1, 1, 0
                        3, 1, 1, 3
                        0, 1, 2, 0
                        1, 1, 2, 4
                        2, 1, 2, 0
                        3, 1, 2, 1
                        0, 2, 1, 1
                        1, 2, 1, 0
                        2, 2, 1, 3
                        3, 2, 1, 1
                        0, 2, 2, 0
                        1, 2, 2, 1
                        2, 2, 2, 3
                        3, 2, 2, 1''')
        self.xyfile7.close()
        self.xymeta7 = {'precision': 1, 'xrange':(1,2), 'yrange':(1,2)}
        self.gridtest2 = Patch('xyfile7.csv')
        self.gridtest2.xy.meta = self.xymeta7
        self.gridtest2.set_attributes()






    def tearDown(self):
        os.remove('xyfile5.csv')
        os.remove('xyfile6.csv')
        os.remove('xyfile7.csv')
        pass

    #
    # init and set_attributes
    #

    def test_patch_init_and_attributes(self):
        self.assertEqual(self.pat.width, 0.2)
        self.assertEqual(self.pat.height, 0.3)
        self.assertEqual(self.pat.x_min, .1)
        self.assertEqual(self.pat.x_max, .3)
        self.assertEqual(self.pat.y_min, .1)
        self.assertEqual(self.pat.y_max, .4)
        self.assertEqual(self.pat.S, 2)
        self.assertEqual(self.pat.N, 6)

        np.testing.assert_array_equal(self.pat.n0_vect, np.array((4, 0, 2)))
        self.assertEqual(np.shape(self.pat.n0_vect)[0], 
                         self.pat.xy.max_spp_code + 1)

    #
    # get_sub_sad and get_sub_abund
    #

    def test_get_sub_sad_and_get_sub_abund(self):
        # Correct result for n0_vect also effectively tests these together
        sa1 = self.pat.get_sub_abund(self.pat.xy.table[[0, 4]])
        np.testing.assert_array_equal(sa1, np.array((2, 0, 1)))

        ss1 = self.pat.get_sub_sad(.1, .3, .1, .3, '')
        ss2 = self.pat.get_sub_sad(.1, .3, .1, .3, 'SAR')
        ss3 = self.pat.get_sub_sad(.1, .3, .1, .3, 'EAR')
        np.testing.assert_array_equal(ss1, np.array((3, 0, 1)))
        np.testing.assert_array_equal(ss2, np.array((2)))
        np.testing.assert_array_equal(ss3, np.array((0)))

        ss1 = self.pat.get_sub_sad(.1, .3, .2, .4, '')
        ss2 = self.pat.get_sub_sad(.1, .3, .2, .4, 'SAR')
        ss3 = self.pat.get_sub_sad(.1, .3, .2, .4, 'EAR')
        np.testing.assert_array_equal(ss1, np.array((2, 0, 2)))
        np.testing.assert_array_equal(ss2, np.array((2)))
        np.testing.assert_array_equal(ss3, np.array((1)))

    #
    # sad_grid
    #

    def test_sad_grid_div_errors(self):
        self.assertRaises(IndexError, self.pat.sad_grid, [(1,0)]) # No 0
        self.assertRaises(IndexError, self.pat.sad_grid, [(1,4)]) # Too big
        self.assertRaises(IndexError, self.pat.sad_grid, [(3,3)]) # Not fact

    def test_sad_grid_answer_full_patch(self):
        np.testing.assert_array_equal(self.pat.sad_grid([(1,1)], '')[0][0], 
                                      self.pat.n0_vect)

    def test_sad_grid_answer_half_patch(self):
        out1 = self.pat.sad_grid([(2,1)], '')[0]
        ans1 = np.array(([4,0,1],[0,0,1]))
        np.testing.assert_array_equal(out1, ans1)

        out2 = self.pat.sad_grid([(2,1)], 'SAR')[0]
        ans2 = np.array((2,1))
        np.testing.assert_array_equal(out2, ans2)

        out3 = self.pat.sad_grid([(2,1)], 'EAR')[0]
        ans3 = np.array((1,0))
        np.testing.assert_array_equal(out3, ans3)

    def test_sad_grid_answer_fully_divided(self):
        out1 = self.pat.sad_grid([(2,3)], '')[0]
        ans1 = np.array(([2,0,0],[1,0,1],[1,0,0],[0,0,0],[0,0,0],[0,0,1]))
        np.testing.assert_array_equal(out1, ans1)

        out2 = self.pat.sad_grid([(2,3)], 'SAR')[0]
        ans2 = np.array((1,2,1,0,0,1))
        np.testing.assert_array_equal(out2, ans2)

        out3 = self.pat.sad_grid([(2,3)], 'EAR')[0]
        ans3 = np.array((0,0,0,0,0,0)) # Note correct ans for spp w 0 individs
        np.testing.assert_array_equal(out3, ans3)

    def test_sad_grid_answer_multiple_divs(self):
        out = self.pat.sad_grid([(2,1),(2,3)], '')
        ans1 = np.array(([4,0,1],[0,0,1]))
        ans2 = np.array(([2,0,0],[1,0,1],[1,0,0],[0,0,0],[0,0,0],[0,0,1]))

        np.testing.assert_array_equal(out[0], ans1)
        np.testing.assert_array_equal(out[1], ans2)

    #
    # misc functions
    #

    def test_divisible(self):
        self.assertTrue(divisible(.3, .1, 3))
        self.assertTrue(divisible(1000, .1, 50))
        self.assertTrue(divisible(16, 1, 4))
        self.assertTrue(divisible(16, 1, 16))

        self.assertFalse(divisible(.3, .1, 4))
        self.assertFalse(divisible(.3, .1, 0))
        self.assertFalse(divisible(1000, .1, 32))
        self.assertFalse(divisible(1000, .1, 128))
        self.assertFalse(divisible(1000, .1, 256))
        self.assertFalse(divisible(16, 1, 5))

    def test_distance(self):
        self.assertEquals(distance((0,0),(1,1)), np.sqrt(2))
        self.assertEquals(distance((1,1),(2,2)), np.sqrt(2))
        self.assertEquals(distance((1,1),(1,6)), 5)

    def test_rnd(self):
        self.assertEquals(rnd(.3), .3)
        self.assertEquals(rnd(.2 + .1), .3)  # .2 + .1 != .3 without round

    def test_QS_grid(self):
        common = self.gridtest.QS_grid([(1,1)])
        self.assertTrue(type(common) == type([]))
        self.assertTrue(len(common) == 1)
        self.assertTrue(len(common[0]) == 0)
        common = self.gridtest.QS_grid([(2,2)])
        self.assertTrue(len(common) == 1)
        self.assertTrue(len(common[0]) == 6)
        chi = np.array([4/5, 2/3, 2/3, 2/5, 4/5, 2/3])
        self.assertTrue(np.array_equal(chi, common[0][:,3]))
        self.assertRaises(IndexError, self.gridtest.QS_grid, [(6,6)])
        self.assertRaises(IndexError, self.gridtest.QS_grid, [(6,0)])
        common = self.gridtest.QS_grid([(1,1),(1,2), (2,1), (2,2)])
        self.assertTrue(len(common) == 4)
        self.assertTrue(np.array_equal(chi, common[3][:,3]))
        self.assertTrue((float(common[1][:,3]) == 6/7) and (float(common[2][:,3]) == 6/7))
        common = self.gridtest2.QS_grid([(1,1),(1,2), (2,1), (2,2)])
        self.assertTrue(len(common) == 4)
        self.assertTrue(np.array_equal(chi, common[3][:,3]))
        self.assertTrue((float(common[1][:,3]) == 6/7) and (float(common[2][:,3]) == 6/7))





if __name__ == '__main__':
    unittest.main()

