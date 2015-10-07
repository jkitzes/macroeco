from __future__ import division
import os
from configparser import ConfigParser

import unittest
from numpy.testing import (TestCase, assert_equal, assert_array_equal,
                           assert_almost_equal, assert_array_almost_equal,
                           assert_allclose, assert_, assert_raises)
from pandas.util.testing import (assert_frame_equal)

import macroeco.empirical as emp
import macroeco.empirical._empirical as _emp
import numpy as np
import pandas as pd
import scipy.stats as stats

# Check whether shapely is installed
try:
    import shapely.geometry as geo
    shapely_missing = False
except:
    shapely_missing = True

class Patches(TestCase):

    def setUp(self):
        local_path = os.path.dirname(os.path.abspath(__file__))

        self.meta1_path = os.path.join(local_path, 'test_meta1.txt')
        self.meta2_path = os.path.join(local_path, 'test_meta2.txt')
        self.table1_path = os.path.join(local_path, 'test_table1.csv')
        self.table1 = pd.DataFrame.from_csv(self.table1_path, index_col=False)
        self.meta1 = ConfigParser()
        self.meta1.read(self.meta1_path)
        self.pat1 = emp.Patch(self.meta1_path)  # No subset
        self.pat2 = emp.Patch(self.meta2_path)  # No subset
        self.cols1 = 'spp_col:spp; count_col:count; x_col:x; y_col:y'
        self.cols2 = 'spp_col:spp; count_col:count; x_col:mean; y_col:y'
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

    def test_subset_count(self):
        # Subsetting on count should work
        pat1 = emp.Patch(self.meta1_path, subset="count > 2")
        assert_equal(pat1.table['count'].iloc[0], 3)
        assert_equal(len(pat1.table), 1)


class TestSAD(Patches):

    def test_simple(self):
        # Falling back on spp_col in metadata, so count 1 for each row
        sad = emp.sad(self.pat1, None, None)
        assert_array_equal(sad[0][1]['y'], [3,2])

    def test_simple_with_cols(self):
        # Specify count and spp_col here
        sad = emp.sad(self.pat1, self.cols1, None)
        assert_array_equal(sad[0][1]['y'], [4,4])

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
        assert_equal(len(sad), 2)
        assert_equal(sad[0][1]['spp'].values, ['a'])
        assert_equal(sad[0][1]['y'].values, [2])
        assert_equal(sad[1][1]['spp'].values, ['a','b'])
        assert_equal(sad[1][1]['y'].values, [2,4])

    def test_split_categorical(self):
        sad = emp.sad(self.pat1, self.cols1, 'year:split; x:2')
        assert_equal(sad[0][1]['y'].values, 3)
        assert_equal(sad[1][1]['y'].values, [])
        assert_equal(sad[2][1]['y'].values, [1,1])
        assert_equal(sad[3][1]['y'].values, [3])

    def test_clean(self):
        # No a in second split on x
        sad = emp.sad(self.pat1, self.cols1, 'x:2', clean=False)
        assert_equal(len(sad[1][1]), 2)  # Both spp when clean False

        sad = emp.sad(self.pat1, self.cols1, 'x:2', clean=True)
        assert_equal(len(sad[1][1]), 1)  # Only 'b' when clean True

    def test_split_panda_default_column_names(self):
        # Columns can be named as key words in pandas
        sad = emp.sad(self.pat2, self.cols2, splits="mean:2", clean=False)
        assert_equal(len(sad[1][1]), 2)

        sad = emp.sad(self.pat2, self.cols2, splits="mean:2; y:3", clean=True)
        assert_equal(len(sad[1][1]), 2)


class TestSSAD(Patches):

    def test_no_splits(self):
        # Just total abundance by species
        ssad = emp.ssad(self.pat1, self.cols1, None)
        assert_array_equal(ssad[0][1]['y'], [4])
        assert_array_equal(ssad[1][1]['y'], [4])

    def test_with_split(self):
        ssad = emp.ssad(self.pat1, self.cols1, 'x:2')
        assert_array_equal(ssad[0][1]['y'], [4,0])  # spp a
        assert_array_equal(ssad[1][1]['y'], [1,3])  # spp b


class TestSAR(Patches):

    def test_no_splits(self):
        sar = emp.sar(self.pat1, self.cols1, None, '1,1; 2,1; 2,3')
        assert_array_almost_equal(sar[0][1]['x'],
                            [1*self.A1, 0.5*self.A1, 1/6*self.A1])
        assert_array_equal(sar[0][1]['y'], [2, 1.5, (1+2+1+0+0+1)/6.])

    def test_with_split(self):
        sar = emp.sar(self.pat1, self.cols1, 'year:split', '2,1; 1,3')
        assert_array_almost_equal(sar[0][1]['x'], [0.5*self.A1, 1/3.*self.A1])
        assert_array_almost_equal(sar[1][1]['x'], [0.5*self.A1, 1/3.*self.A1])
        assert_array_equal(sar[0][1]['y'], [0.5, 2/3.])
        assert_array_equal(sar[1][1]['y'], [3/2., 1])

    def test_single_division(self):
        sar = emp.sar(self.pat1, self.cols1, None, '2,1')
        assert_array_almost_equal(sar[0][1]['x'], [0.5*self.A1])
        assert_array_equal(sar[0][1]['y'], [1.5])

    def test_empty_equals_split_subset(self):
        sar_empty = emp.sar(self.pat1, self.cols1, "", '1,1')
        sar_split = emp.sar(self.pat1, self.cols1, "x:1; y:1", '1,1')
        print sar_empty
        print sar_split
        assert_frame_equal(sar_empty[0][1].sort(axis=1),
                                            sar_split[0][1].sort(axis=1))


class TestEAR(Patches):

    def test_no_splits(self):
        sar = emp.sar(self.pat1, self.cols1, None, '1,1; 2,1; 2,3', ear=True)
        assert_array_equal(sar[0][1]['y'], [2, 0.5, 0])

    def test_with_split(self):
        sar = emp.sar(self.pat1, self.cols1, 'year:split', '2,1;1,3', ear=True)
        assert_array_equal(sar[0][1]['y'], [0.5, 0])
        assert_array_equal(sar[1][1]['y'], [0.5, 1/3.])


class TestCommGrid(Patches):

    def test_no_splits_Sorensen(self):
        comm = emp.comm_grid(self.pat1, self.cols1, None, '2,1')
        assert_almost_equal(comm[0][1]['x'], [0.1])
        assert_array_equal(comm[0][1]['y'], [2./(2+1)])

    def test_no_splits_Jaccard(self):
        comm = emp.comm_grid(self.pat1, self.cols1, None, '2,1',
                             metric='Jaccard')
        assert_almost_equal(comm[0][1]['x'], [0.1])
        assert_array_equal(comm[0][1]['y'], [1/2.])

    def test_with_split(self):
        comm = emp.comm_grid(self.pat1, self.cols1, 'year:split', '2,1')
        assert_array_equal(comm[0][1]['y'], [0])
        assert_array_equal(comm[1][1]['y'], [2/3.])

    def test_y_division_even(self):
        comm = emp.comm_grid(self.pat1, self.cols1, '', '1,3')
        assert_array_equal(comm[0][1]['pair'], ['(0.15 0.1) - (0.15 0.2)',
                                          '(0.15 0.1) - (0.15 0.3)',
                                          '(0.15 0.2) - (0.15 0.3)'])
        assert_array_almost_equal(comm[0][1]['x'], [0.1, 0.2, 0.1])
        assert_array_equal(comm[0][1]['y'], [2/3., 2/3., 1.])

    def test_x_y_division_uneven_y(self):
        comm = emp.comm_grid(self.pat1, self.cols1, '', '2,2')
        print comm
        assert_array_equal(comm[0][1]['pair'], ['(0.1 0.125) - (0.1 0.275)',
                                          '(0.1 0.125) - (0.2 0.125)',
                                          '(0.1 0.125) - (0.2 0.275)',
                                          '(0.1 0.275) - (0.2 0.125)',
                                          '(0.1 0.275) - (0.2 0.275)',
                                          '(0.2 0.125) - (0.2 0.275)'])
        assert_array_almost_equal(comm[0][1]['x'], [0.15, 0.1, 0.180278, 0.180278,
                                              0.1, 0.15], 6)
        assert_array_equal(comm[0][1]['y'], [2/3., 0, 0, 0, 2/3., 0])

    def test_x_y_division_uneven_y_jaccard(self):
        comm = emp.comm_grid(self.pat1, self.cols1, '', '2,2',metric='Jaccard')
        assert_array_equal(comm[0][1]['y'], [1/2., 0, 0, 0, 1/2., 0])


@unittest.skipIf(shapely_missing, "shapely not present, skipping O-ring test")
class TestORing(Patches):
    # TODO: Main may fail with error if dataframe has no records when trying to
    # fit or make plot.

    def test_spp_no_present_returns_empty_df(self):
        o_ring = emp.o_ring(self.pat1, self.cols1, '', 'nothere', [0,.1,.2])
        assert_frame_equal(o_ring[0][1], pd.DataFrame(columns=['x','y']))

    def test_one_individual_returns_zeros(self):
        self.pat1.table = self.pat1.table[2:4]  # Leave 1 'a' and 1 'b'
        o_ring = emp.o_ring(self.pat1, self.cols1, '', 'a', [0,.1,.2])
        assert_array_equal(o_ring[0][1]['y'], [0, 0])

    def test_no_density_a(self):
        # Points on bin edge may be allocated ambiguously due to floating point
        # issues - testing here with slightly offset edges
        o_ring = emp.o_ring(self.pat1, self.cols1, '', 'a', [0,.101,.201,.301],
                            density=False)
        assert_array_almost_equal(o_ring[0][1]['x'], [0.0505, 0.151, 0.251])
        assert_array_almost_equal(o_ring[0][1]['y'], [8, 4, 0])

    def test_no_density_b(self):
        o_ring = emp.o_ring(self.pat1, self.cols1, '', 'b', [0,.1,.2,.3],
                            density=False)
        assert_array_almost_equal(o_ring[0][1]['x'], [0.05, 0.15,0.25])
        assert_array_almost_equal(o_ring[0][1]['y'], [6, 6, 0])

    def test_with_split_a(self):
        o_ring = emp.o_ring(self.pat1, self.cols1, 'y:2', 'a', [0,.1,.2],
                            density=False)
        assert_array_equal(o_ring[0][1]['y'], [2, 0])  # Bottom
        assert_array_equal(o_ring[1][1]['y'], [2, 0])  # Top

    def test_with_split_b(self):
        o_ring = emp.o_ring(self.pat1, self.cols1, 'y:2', 'b', [0,.1,.2],
                            density=False)
        assert_array_equal(o_ring[0][1]['y'], [])  # Bottom
        assert_array_equal(o_ring[1][1]['y'], [6, 6])  # Top

    def test_density_a(self):
        # First radius is 0.05
        o_ring = emp.o_ring(self.pat1, self.cols1, '', 'a', [0,.10000001])
        assert_array_almost_equal(o_ring[0][1]['y'],
                                  [8 / (1.25*np.pi*(0.1)**2)],
                                  3)

    def test_density_b(self):
        # First radius is 0.05
        o_ring = emp.o_ring(self.pat1, self.cols1, '', 'b', [0,.10000001,.1828427])
        assert_array_almost_equal(o_ring[0][1]['y'],
                                  [6 / (1.25*np.pi*(0.1)**2),
                                   6 / (3/8 * np.pi*(0.1828427**2 - 0.1**2))],
                                  3)


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

