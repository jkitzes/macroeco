'''
Unit tests for empirical.py
'''

import unittest
import os
from empirical import *


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


    def tearDown(self):
        #os.remove('xyfile5.csv')
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


if __name__ == '__main__':
    unittest.main()

