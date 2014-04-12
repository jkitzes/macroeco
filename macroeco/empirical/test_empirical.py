from __future__ import division
import os
from configparser import ConfigParser

from numpy.testing import (TestCase, assert_equal, assert_array_equal,
                           assert_almost_equal, assert_array_almost_equal,
                           assert_allclose, assert_, assert_raises)
from pandas.util.testing import (assert_frame_equal)

import macroeco.empirical as emp
import macroeco.empirical._empirical as _emp
import numpy as np
import pandas as pd
import scipy.stats as stats


class Patches(TestCase):

    def setUp(self):
        local_path = os.path.dirname(os.path.abspath(__file__))

        self.meta1_path = os.path.join(local_path, 'test_meta1.txt')
        self.table1_path = os.path.join(local_path, 'test_table1.csv')
        self.table1 = pd.DataFrame.from_csv(self.table1_path, index_col=False)
        self.meta1 = ConfigParser()
        self.meta1.read(self.meta1_path)
        self.pat1 = emp.Patch(self.meta1_path)  # No subset
        self.cols1 = 'spp_col:spp; count_col:count; x_col:x; y_col:y'
        self.A1 = 0.2 * 0.3


class TestPatch(Patches):

    def test_load_data_meta(self):
        assert_array_equal(self.pat1.table, self.table1)
        assert_equal(self.pat1.meta, self.meta1)

    def test_subset_numeric(self):
        pat1 = emp.Patch(self.meta1_path, 'x>=0.2')
        assert_array_equal(pat1.table, self.table1[self.table1.x >= 0.2])

        self.meta1['x']['min'] = '0.2'
        assert_equal(pat1.meta, self.meta1)

    def test_subset_categorical(self):
        pat1 = emp.Patch(self.meta1_path, "spp=='b'")
        assert_array_equal(pat1.table, self.table1[self.table1['spp']=='b'])
        assert_equal(pat1.meta, self.meta1)  # Meta should not change

    def test_multiple_subset(self):
        # Only first element in table remains
        pat1 = emp.Patch(self.meta1_path, "spp=='a' ; y < 0.2")
        assert_array_equal(pat1.table.iloc[0], self.table1.iloc[0])
        assert_equal(len(pat1.table), 1)

        self.meta1['y']['max'] = '0.1'
        assert_equal(pat1.meta, self.meta1)


class TestSAD(Patches):

    def test_simple(self):
        # Falling back on spp_col in metadata, so count 1 for each row
        sad = emp.sad(self.pat1, None, None)
        assert_equal(sad[0][1]['y'], [3,2])

    def test_simple_with_cols(self):
        # Specify count and spp_col here
        sad = emp.sad(self.pat1, self.cols1, None)
        assert_equal(sad[0][1]['y'], [4,3])

    def test_two_way_split(self):
        # Complete split generates 6 results
        sad = emp.sad(self.pat1, self.cols1, 'x:2; y:3')
        assert_equal(len(sad), 6)

        # Goes through x then y
        assert_equal(sad[0][1]['spp'].values, 'a')
        assert_equal(sad[0][1]['y'].values, 2)
        assert_equal(sad[1][1]['y'].values, [1,1])
        assert_equal(sad[5][1]['spp'].values, 'b')
        assert_equal(sad[0][1]['y'].values, 2)

    def test_one_way_uneven_split(self):
        # 0.2 should fall in second division of y
        sad = emp.sad(self.pat1, self.cols1, 'y:2')
        print sad
        assert_equal(len(sad), 2)
        assert_equal(sad[0][1]['spp'].values, ['a'])
        assert_equal(sad[0][1]['y'].values, [2])
        assert_equal(sad[1][1]['spp'].values, ['a','b'])
        assert_equal(sad[1][1]['y'].values, [2,3])

    def test_split_categorical(self):
        sad = emp.sad(self.pat1, self.cols1, 'year:split; x:2')
        assert_equal(sad[0][1]['y'].values, 3)
        assert_equal(sad[1][1]['y'].values, [])
        assert_equal(sad[2][1]['y'].values, [1,1])
        assert_equal(sad[3][1]['y'].values, [2])

    def test_clean(self):
        # No a in second split on x
        sad = emp.sad(self.pat1, self.cols1, 'x:2', clean=False)
        assert_equal(len(sad[1][1]), 2)  # Both spp when clean False

        sad = emp.sad(self.pat1, self.cols1, 'x:2', clean=True)
        assert_equal(len(sad[1][1]), 1)  # Only 'b' when clean True


class TestSSAD(Patches):

    def test_no_splits(self):
        # Just total abundance by species
        ssad = emp.ssad(self.pat1, self.cols1, None)
        assert_equal(ssad[0][1]['y'], [4])
        assert_equal(ssad[1][1]['y'], [3])

    def test_with_split(self):
        ssad = emp.ssad(self.pat1, self.cols1, 'x:2')
        assert_equal(ssad[0][1]['y'], [4,0])  # spp a
        assert_equal(ssad[1][1]['y'], [1,2])  # spp b


class TestSAR(Patches):

    def test_no_splits(self):
        sar = emp.sar(self.pat1, self.cols1, None, '1,1; 2,1; 2,3')
        assert_almost_equal(sar[0][1]['x'],
                            [1*self.A1, 0.5*self.A1, 1/6*self.A1])
        assert_equal(sar[0][1]['y'], [2, 1.5, (1+2+1+0+0+1)/6.])

    def test_with_split(self):
        sar = emp.sar(self.pat1, self.cols1, 'year:split', '2,1; 1,3')
        assert_almost_equal(sar[0][1]['x'], [0.5*self.A1, 1/3.*self.A1])
        assert_almost_equal(sar[1][1]['x'], [0.5*self.A1, 1/3.*self.A1])
        assert_equal(sar[0][1]['y'], [0.5, 2/3.])
        assert_equal(sar[1][1]['y'], [3/2., 1])

    def test_single_division(self):
        sar = emp.sar(self.pat1, self.cols1, None, '2,1')
        assert_almost_equal(sar[0][1]['x'], [0.5*self.A1])
        assert_equal(sar[0][1]['y'], [1.5])


class TestEAR(Patches):

    def test_no_splits(self):
        sar = emp.sar(self.pat1, self.cols1, None, '1,1; 2,1; 2,3', ear=True)
        assert_equal(sar[0][1]['y'], [2, 0.5, 0])

    def test_with_split(self):
        sar = emp.sar(self.pat1, self.cols1, 'year:split', '2,1;1,3', ear=True)
        assert_equal(sar[0][1]['y'], [0.5, 0])
        assert_equal(sar[1][1]['y'], [0.5, 1/3.])


class TestCommGrid(Patches):

    def test_no_splits_Sorensen(self):
        comm = emp.comm_grid(self.pat1, self.cols1, None, '2,1')
        assert_almost_equal(comm[0][1]['x'], [0.1])
        assert_equal(comm[0][1]['y'], [2./(2+1)])

    def test_no_splits_Jaccard(self):
        comm = emp.comm_grid(self.pat1, self.cols1, None, '2,1',
                             metric='Jaccard')
        assert_almost_equal(comm[0][1]['x'], [0.1])
        assert_equal(comm[0][1]['y'], [1/2.])

    def test_with_split(self):
        comm = emp.comm_grid(self.pat1, self.cols1, 'year:split', '2,1')
        assert_equal(comm[0][1]['y'], [0])
        assert_equal(comm[1][1]['y'], [2/3.])

    def test_y_division_even(self):
        comm = emp.comm_grid(self.pat1, self.cols1, '', '1,3')
        assert_equal(comm[0][1]['pair'], ['(0.15 0.1) - (0.15 0.2)',
                                          '(0.15 0.1) - (0.15 0.3)',
                                          '(0.15 0.2) - (0.15 0.3)'])
        assert_almost_equal(comm[0][1]['x'], [0.1, 0.2, 0.1])
        assert_equal(comm[0][1]['y'], [2/3., 2/3., 1.])

    def test_x_y_division_uneven_y(self):
        comm = emp.comm_grid(self.pat1, self.cols1, '', '2,2')
        print comm
        assert_equal(comm[0][1]['pair'], ['(0.1 0.125) - (0.1 0.275)',
                                          '(0.1 0.125) - (0.2 0.125)',
                                          '(0.1 0.125) - (0.2 0.275)',
                                          '(0.1 0.275) - (0.2 0.125)',
                                          '(0.1 0.275) - (0.2 0.275)',
                                          '(0.2 0.125) - (0.2 0.275)'])
        assert_almost_equal(comm[0][1]['x'], [0.15, 0.1, 0.180278, 0.180278,
                                              0.1, 0.15], 6)
        assert_equal(comm[0][1]['y'], [2/3., 0, 0, 0, 2/3., 0])

    def test_x_y_division_uneven_y_jaccard(self):
        comm = emp.comm_grid(self.pat1, self.cols1, '', '2,2',metric='Jaccard')
        assert_equal(comm[0][1]['y'], [1/2., 0, 0, 0, 1/2., 0])


class TestProduct():

    def test_product_with_order(self):
        # Several places rely on product to sequentially loop first -> last
        expected = [[1,5], [1,6], [1,7], [2,5], [2,6], [2,7]]
        assert_equal(_emp._product([1,2],[5,6,7]), expected)


class TestDistance():

    def test_cartesian_distance(self):
        assert_equal(_emp._distance((0,0),(2,2)), np.sqrt(8))


class TestDecDegDistance():

    def test_ucberkeley_to_sf(self):
        # Latlong: http://www.findlatitudeandlongitude.com
        # Dist: http://www.movable-type.co.uk/scripts/latlong.html (17.37 km)
        berkeley = (37.87133, -122.259293)
        sf = (37.780213, -122.419968)
        assert_almost_equal(_emp._decdeg_distance(berkeley, sf), 17.37, 1)


class TestEmpiricalCDF():

    def test_sorted_data(self):
        test_data = [1, 1, 1, 1, 2, 3, 4, 5, 6, 6]
        ans = [.4, .4, .4, .4, .5, .6, .7, .8, 1, 1]
        res = emp.empirical_cdf(test_data)
        assert_array_equal(ans, res['ecdf'])

    def test_unsorted_data(self):
        test_data = [6, 6, 1, 1, 5, 1, 1, 2, 3, 4]
        ans = [.4, .4, .4, .4, .5, .6, .7, .8, 1, 1]
        res = emp.empirical_cdf(test_data)
        assert_array_equal(ans, res['ecdf'])  # Result sorted
        assert_array_equal(np.sort(test_data), res['data'])  # Data sorted

    def test_all_data_same(self):
        test_data = [3, 3, 3, 3]
        ans = [1, 1, 1, 1]
        res = emp.empirical_cdf(test_data)
        assert_array_equal(ans, res['ecdf'])




#     def test_comm_sep(self):

#         # Create result recarray
#         comm = self.pat9.comm_sep({'plot1': (0,0), 'plot2': (0,1),
#                                    'plot3': (3,4)},
#                                   {'spp_code': 'species', 'count': 'count'})

#         # Create result recarray with dec degree locs
#         comm_decdeg = self.pat9.comm_sep({'plot1': (9.1,79.0),
#                                    'plot2': (9.2,79.5), 'plot3': (12.7,50)},
#                                   {'spp_code': 'species', 'count': 'count'},
#                                          loc_unit='decdeg')

#         # Check distances
#         dist_sort = np.sort(comm['dist'])
#         np.testing.assert_array_almost_equal(dist_sort, np.array((1,4.242,5)),
#                                              3)

#         # Check distances dec degree
#         # TODO: Find exact third party comparison formula - formulas online use
#         # different radii, etc. and give approx same answer
#         dist_sort = np.sort(comm_decdeg['dist'])
#         #np.testing.assert_array_almost_equal(dist_sort,
#         #                                     np.array((56.058,3193.507,
#         #                                               3245.820)), 3)

#         # Check species in each plot
#         spp_sort = np.sort(np.array(list(comm['spp-a']) + list(comm['spp-b'])))
#         np.testing.assert_array_equal(spp_sort, np.array((0,0,3,3,4,4)))

#         # Check Sorensen - 2 zeros from empty plot1
#         sor_sort = np.sort(comm['sorensen'])
#         np.testing.assert_array_almost_equal(sor_sort,
#                                              np.array((0,0,0.571428571)), 5)

#         # Check Jaccard - 2 zeros from empty plot1
#         jac_sort = np.sort(comm['jaccard'])
#         np.testing.assert_array_almost_equal(jac_sort, np.array((0,0,0.4)), 5)

#     def test_o_ring(self):

#         # Check standard case, no min max, no edge correction, no criteria
#         # Tests that distances and repeats for count col are correct
#         result_list = self.pat1.o_ring(('x','y'), [0,.11,.2],
#                                      {'spp_code': 'species', 'count': 'count'})

#         np.testing.assert_array_equal(result_list[0][2][0], np.array((8,4)))
#         np.testing.assert_array_equal(result_list[0][2][1], np.array((2,4)))

#         # Check standard case, no min max, no edge correction, with division
#         result_list = self.pat1.o_ring(('x','y'), [0,.11,.2],
#                                      {'spp_code': 'species', 'count': 'count',
#                                       'y': 2})

#         # - First half of y, both species
#         np.testing.assert_array_equal(result_list[0][2][0], np.array((6,0)))
#         np.testing.assert_array_equal(result_list[0][2][1], np.array((0,0)))

#         # - Second half of y, both species
#         np.testing.assert_array_equal(result_list[1][2][0], np.array((0,0)))
#         np.testing.assert_array_equal(result_list[1][2][1], np.array((2,0)))

#         # Check edge correction - check only first species
#         # Almost equal required due to float division
#         result_list = self.pat1.o_ring(('x','y'), [0,.05,.1],
#                                      {'spp_code': 'species', 'count': 'count'},
#                                           edge_correct=True)
#         np.testing.assert_array_almost_equal(result_list[0][2][0],
#                                              np.array((8,18)))

#         # Check density - check only second species
#         print 'here '
#         result_list = self.pat1.o_ring(('x','y'), [0,.05,.1],
#                                      {'spp_code': 'species', 'count': 'count'},
#                                           density=True)
#         np.testing.assert_array_almost_equal(result_list[0][2][1],
#                                              np.array((1358.12218105,0)))
