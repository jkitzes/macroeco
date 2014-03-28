'''
Unit tests for empirical.py
'''

from __future__ import division

from numpy.testing import (TestCase, assert_equal, assert_array_equal,
                           assert_almost_equal, assert_array_almost_equal,
                           assert_allclose, assert_, assert_raises)

from macroeco.empirical import *
import numpy as np
import scipy.stats as stats
import numpy.testing as nt

class TestEmpiricalCDF(TestCase):
    """ Unittests for Empirical cdf """

    def test_empirical_cdf_vs_R(self):

        #Test against R's ecdf function

        # Test Case 1
        test_data = [1, 1, 1, 1, 2, 3, 4, 5, 6, 6]
        R_res = [.4, .4, .4, .4, .5, .6, .7, .8, 1, 1]
        res = empirical_cdf(test_data)
        assert_array_equal(R_res, res)

        # Test Case 2
        test_data = [3, 3, 3, 3]
        R_res = [1, 1, 1, 1]
        res = empirical_cdf(test_data)
        assert_array_equal(R_res, res)

# class TestPatch(unittest.TestCase):

#     def setUp(self):
#         self.xyfile5 = open('xyfile5.csv','w')
#         self.xyfile5.write('''spp_code, x, y, count
# grt, .1, .1, 2
# grt, .1, .2, 1
# grt, .1, .3, 1
# rty, .1, .2, 1
# rty, .2, .3, 2''')
#         self.xyfile5.close()
#         self.xymeta5 = {('x', 'maximum'): .2, ('x', 'minimum'): .1, ('x',
#         'precision'): .1, ('x', 'type'): 'interval', ('y', 'maximum'): .3,
#         ('y', 'minimum'): .1, ('y', 'precision'): .1, ('y', 'type'): 'interval',
#         ('spp_code', 'maximum'): None, ('spp_code', 'minimum'): None,
#         ('spp_code', 'precision'): None, ('spp_code', 'type'): 'ordinal',
#         ('count', 'maximum'): None, ('count', 'minimum'): None, ('count',
#         'precision'): None, ('count', 'type'): 'ratio'}

#         self.pat1 = Patch('xyfile5.csv')
#         # Line below sets metadata manually-no metadata file loaded
#         self.pat1.data_table.meta = self.xymeta5

#         self.xyfile6 = open('xyfile6.csv', 'w')
#         self.xyfile6.write('''spp_code, x, y, count
# a, 0, 0, 1
# b, 0, 0, 1
# c, 0, 0, 0
# d, 0, 0, 3
# a, 0, 1, 0
# b, 0, 1, 4
# c, 0, 1, 0
# d, 0, 1, 1
# a, 1, 0, 1
# b, 1, 0, 0
# c, 1, 0, 3
# d, 1, 0, 1
# a, 1, 1, 0
# b, 1, 1, 1
# c, 1, 1, 3
# d, 1, 1, 1''')
#         self.xyfile6.close()
#         self.xymeta6 = {('x', 'maximum'): 1, ('x', 'minimum'): 0, ('x',
#         'precision'): 1, ('x', 'type'): 'interval', ('y', 'maximum'): 1,
#         ('y', 'minimum'): 0, ('y', 'precision'): 1, ('y', 'type'): 'interval',
#         ('spp_code', 'maximum'): None, ('spp_code', 'minimum'): None,
#         ('spp_code', 'precision'): None, ('spp_code', 'type'): 'ordinal',
#         ('count', 'maximum'): None, ('count', 'minimum'): None, ('count',
#         'precision'): None, ('count', 'type'): 'ratio'}
#         self.pat2 = Patch('xyfile6.csv')
#         self.pat2.data_table.meta = self.xymeta6

#         self.xyfile7 = open('xyfile7.csv', 'w')
#         self.xyfile7.write('''spp_code, x, y, count
# tery, 1, 1, 1
# 1, 1, 1, 1
# 2, 1, 1, 0
# 3, 1, 1, 3
# 0, 1, 2, 0
# 1, 1, 2, 4
# 2, 1, 2, 0
# tery, 1, 2, 1
# 0, 2, 1, 1
# 1, 2, 1, 0
# 2, 2, 1, 3
# 3, 2, 1, 1
# tery, 2, 2, 0
# 1, 2, 2, 1
# 2, 2, 2, 3
# 3, 2, 2, 1''')
#         self.xyfile7.close()
#         self.xymeta7 = {('x', 'maximum'): 2, ('x', 'minimum'): 1, ('x',
#         'precision'): 1, ('x', 'type'): 'interval', ('y', 'maximum'): 2,
#         ('y', 'minimum'): 1, ('y', 'precision'): 1, ('y', 'type'): 'interval',
#         ('spp_code', 'maximum'): None, ('spp_code', 'minimum'): None,
#         ('spp_code', 'precision'): None, ('spp_code', 'type'): 'ordinal',
#         ('count', 'maximum'): None, ('count', 'minimum'): None, ('count',
#         'precision'): None, ('count', 'type'): 'ratio'}
#         self.pat3 = Patch('xyfile7.csv')
#         self.pat3.data_table.meta = self.xymeta7

#         self.xyfile8 = open('xyfile8.csv', 'w')
#         self.xyfile8.write('''spp_code, x, y, count
# 0, 0, 0, 1
# 1, 0, 0, 1
# 2, 0, 0, 0
# 3, 0, 0, 3
# 0, 0, 1, 0
# 1, 0, 1, 4
# 2, 0, 1, 0
# 3, 0, 1, 1
# 0, 1, 0, 1
# 1, 1, 0, 0
# 2, 1, 0, 3
# 3, 1, 0, 1
# 0, 1, 1, 0
# 1, 1, 1, 1
# 2, 1, 1, 3
# 3, 1, 1, 1
# 0, 2, 0, 0
# 1, 2, 0, 0
# 2, 2, 0, 2
# 3, 2, 0, 4
# 0, 2, 1, 0
# 1, 2, 1, 0
# 2, 2, 1, 0
# 3, 2, 1, 1''')
#         self.xyfile8.close()
#         self.xymeta8 = {('x', 'maximum'): 2, ('x', 'minimum'): 0, ('x',
#         'precision'): 1, ('x', 'type'): 'interval', ('y', 'maximum'): 1,
#         ('y', 'minimum'): 0, ('y', 'precision'): 1, ('y', 'type'): 'interval',
#         ('spp_code', 'maximum'): None, ('spp_code', 'minimum'): None,
#         ('spp_code', 'precision'): None, ('spp_code', 'type'): 'ordinal',
#         ('count', 'maximum'): None, ('count', 'minimum'): None, ('count',
#         'precision'): None, ('count', 'type'): 'ratio'}
#         self.pat4 = Patch('xyfile8.csv')
#         self.pat4.data_table.meta = self.xymeta8
#         self.xyfile9 = open('xyfile9.csv','w')
#         self.xyfile9.write('''spp_code, x, y, count, energy, mass
# grt, .1, .1, 2, 1, 34
# grt, .1, .2, 1, 2, 12
# grt, .1, .3, 1, 3, 23
# rty, .1, .2, 1, 4, 45
# rty, .2, .3, 1, 5, 110''')
#         self.xyfile9.close()
#         self.xymeta9 = {('x', 'maximum'): .2, ('x', 'minimum'): .1, ('x',
#         'precision'): .1, ('x', 'type'): 'interval', ('y', 'maximum'): .3,
#         ('y', 'minimum'): .1, ('y', 'precision'): .1, ('y', 'type'): 'interval',
#         ('spp_code', 'maximum'): None, ('spp_code', 'minimum'): None,
#         ('spp_code', 'precision'): None, ('spp_code', 'type'): 'ordinal',
#         ('count', 'maximum'): None, ('count', 'minimum'): None, ('count',
#         'precision'): None, ('count', 'type'): 'ratio'}

#         self.pat5 = Patch('xyfile9.csv')
#         self.pat5.data_table.meta = self.xymeta9
#         self.xyfile10 = open('xyfile10.csv', 'w')
#         self.xyfile10.write('''spp_code, x, y, count
# a, 0, 0, 1
# b, 0, 0, 1
# d, 0, 0, 3
# b, 0, 1, 4
# d, 0, 1, 1
# a, 1, 0, 1
# c, 1, 0, 3
# d, 1, 0, 1
# b, 1, 1, 1
# c, 1, 1, 3
# d, 1, 1, 1''')
#         self.xyfile10.close()
#         self.xymeta10 = {('x', 'maximum'): 1, ('x', 'minimum'): 0, ('x',
#         'precision'): 1, ('x', 'type'): 'interval', ('y', 'maximum'): 1,
#         ('y', 'minimum'): 0, ('y', 'precision'): 1, ('y', 'type'): 'interval',
#         ('spp_code', 'maximum'): None, ('spp_code', 'minimum'): None,
#         ('spp_code', 'precision'): None, ('spp_code', 'type'): 'ordinal',
#         ('count', 'maximum'): None, ('count', 'minimum'): None, ('count',
#         'precision'): None, ('count', 'type'): 'ratio'}
#         self.pat6 = Patch('xyfile10.csv')
#         self.pat6.data_table.meta = self.xymeta10
#         self.xyfile11 = open('xyfile11.csv', 'w')
#         self.xyfile11.write('''spp_code, x, y, count, reptile
# a, 0, 0, 1, lizard
# b, 0, 0, 1, lizard
# d, 0, 0, 3, snake
# b, 0, 1, 4, lizard
# d, 0, 1, 1, turtle
# a, 1, 0, 1, snake
# c, 1, 0, 3, lizard
# d, 1, 0, 1, snake
# b, 1, 1, 1, tuatara
# c, 1, 1, 3, turtle
# d, 1, 1, 1, snake''')
#         self.xyfile11.close()
#         self.xymeta11 = {('x', 'maximum'): 1, ('x', 'minimum'): 0, ('x',
#         'precision'): 1, ('x', 'type'): 'interval', ('y', 'maximum'): 1,
#         ('y', 'minimum'): 0, ('y', 'precision'): 1, ('y', 'type'): 'interval',
#         ('spp_code', 'maximum'): None, ('spp_code', 'minimum'): None,
#         ('spp_code', 'precision'): None, ('spp_code', 'type'): 'ordinal',
#         ('count', 'maximum'): None, ('count', 'minimum'): None, ('count',
#         'precision'): None, ('count', 'type'): 'ratio', ('reptile', 'maximum')
#         : None, ('reptile', 'minimum') : None, ('reptile', 'precision'):None,
#         ('reptile', 'type') : 'ordinal'}
#         self.pat7 = Patch('xyfile11.csv')
#         self.pat7.data_table.meta = self.xymeta11

#         self.xyfile12 = open('xyfile12.csv', 'w')
#         self.xyfile12.write('''spp_code, x, y, count
# 3, 0, 0, 3
# 3, 0, 1, 1
# 2, 0, 2, 3
# 1, 0, 3, 8
# 3, 1, 0, 1
# 3, 1, 1, 1
# 0, 1, 2, 5
# 3, 1, 3, 1
# 2, 2, 0, 1
# 1, 2, 1, 3
# 1, 2, 2, 6
# 0, 2, 3, 1
# 1, 3, 0, 9
# 2, 3, 1, 1
# 0, 3, 2, 3
# 3, 3, 3, 1''')
#         self.xyfile12.close()
#         self.xymeta12 = {('x', 'maximum'): 3, ('x', 'minimum'): 0, ('x',
#         'precision'): 1, ('x', 'type'): 'interval', ('y', 'maximum'): 3,
#         ('y', 'minimum'): 0, ('y', 'precision'): 1, ('y', 'type'): 'interval',
#         ('spp_code', 'maximum'): None, ('spp_code', 'minimum'): None,
#         ('spp_code', 'precision'): None, ('spp_code', 'type'): 'ordinal',
#         ('count', 'maximum'): None, ('count', 'minimum'): None, ('count',
#         'precision'): None, ('count', 'type'): 'ratio'}
#         self.pat8 = Patch('xyfile12.csv')
#         self.pat8.data_table.meta = self.xymeta12

#         # Data file with three count colums, unique row for each species
#         self.xyfile13 = open('xyfile13.csv', 'w')
#         self.xyfile13.write('''spp_code, order, plot1, plot2, plot3
# a, pred, 0, 0, 0
# b, pred, 0, 0, 1
# c, pred, 0, 1, 0
# d, pred, 0, 2, 3
# e, scav, 0, 1, 0
# f, scav, 0, 1, 4''')
#         self.xyfile13.close()
#         self.xymeta13 = {('spp_code', 'maximum'): None,
#                          ('spp_code', 'minimum'): None,
#                          ('spp_code', 'precision'): None,
#                          ('spp_code', 'type'): 'ordinal',
#                          ('order', 'maximum'): None,
#                          ('order', 'minimum'): None,
#                          ('order', 'precision'): None,
#                          ('order', 'type'): 'ordinal',
#                          ('plot1', 'maximum'): None,
#                          ('plot1', 'minimum'): None,
#                          ('plot1', 'precision'): None,
#                          ('plot1', 'type'): 'ratio',
#                          ('plot2', 'maximum'): None,
#                          ('plot2', 'minimum'): None,
#                          ('plot2', 'precision'): None,
#                          ('plot2', 'type'): 'ratio',
#                          ('plot3', 'maximum'): None,
#                          ('plot3', 'minimum'): None,
#                          ('plot3', 'precision'): None,
#                          ('plot3', 'type'): 'ratio'}
#         self.pat9 = Patch('xyfile13.csv')
#         self.pat9.data_table.meta = self.xymeta13




#     def tearDown(self):
#         os.remove('xyfile5.csv')
#         os.remove('xyfile6.csv')
#         os.remove('xyfile7.csv')
#         os.remove('xyfile8.csv')
#         os.remove('xyfile9.csv')
#         os.remove('xyfile10.csv')
#         os.remove('xyfile11.csv')
#         os.remove('xyfile12.csv')
#         os.remove('xyfile13.csv')

#     #
#     # init and set_attributes
#     #

#     def test_patch_init(self):

#         # Test entire table is loaded
#         self.assertTrue(len(self.pat1.data_table.table) == 5)
#         self.assertTrue(len(self.pat2.data_table.table) == 16)

#         # Test that subsetting works
#         pat = Patch('xyfile6.csv', {'spp_code': [('!=','a'), ('!=', 'b'),
#                                     ('!=','c')]})
#         self.assertTrue(np.all(pat.data_table.table['spp_code'] == 'd'))
#         pat = Patch('xyfile7.csv', {'spp_code': ('==', "tery")})
#         self.assertTrue(sum(pat.data_table.table['count']) == 2)

#         # Testing that metadata was set correctly
#         self.assertTrue(self.pat1.data_table.meta[('x', 'maximum')] == .2)

#     def test_sad(self):

#         # Test correct result with 'whole' and one division
#         sad = self.pat1.sad({'spp_code': 'species', 'count': 'count',
#                                                                     'x': 1})
#         self.assertTrue(np.array_equal(sad[0][1], np.array([4,3])))
#         sad = self.pat1.sad({'spp_code': 'species', 'count': 'count',
#                                                         'x': 'whole'})
#         self.assertTrue(np.array_equal(sad[0][1], np.array([4,3])))
#         sad = self.pat4.sad({'spp_code': 'species', 'count' :'count', 'x': 1})
#         self.assertTrue(np.array_equal(sad[0][2], np.array([0,1,2,3])))

#         # Test correct result with other divisions
#         sad = self.pat4.sad({'spp_code': 'species', 'count': 'count', 'x': 3,
#         'y': 2})
#         self.assertTrue(np.array_equal(sad[-1][1], np.array([0,0,0,1])))

#         # Test that 'whole' and ignore give the same result
#         sad1 = self.pat4.sad({'spp_code': 'species', 'count': 'count'})
#         sad2 = self.pat4.sad({'spp_code': 'species', 'count': 'count', 'x' :
#         'whole'})
#         self.assertTrue(np.array_equal(sad1[0][1], sad2[0][1]))

#         # Test that 'split' keyword returns the correct results
#         sad = self.pat5.sad({'spp_code' :'species', 'energy':'split', 'count'
#                              : 'count'})
#         self.assertTrue(len(sad) == 5)
#         self.assertTrue(np.array_equal(sad[0][1], np.array([2,0])))

#         # Test split and clean on numeric column
#         sad = self.pat5.sad({'spp_code' :'species', 'energy':'split', 'count'
#                              : 'count'}, clean=True)
#         self.assertTrue(len(sad) == 5)
#         self.assertTrue(np.array_equal(sad[0][1], np.array([2])))

#         # Test that cleaning sad and split works on string
#         sad = self.pat7.sad({'spp_code' : 'species', 'count' : 'count',
#                              'reptile' : 'split'}, clean=True)
#         self.assertTrue(len(sad) == 4)
#         self.assertTrue(np.array_equal(sad[0][1], np.array([1,5,3])))
#         self.assertTrue(np.array_equal(sad[2][1], np.array([1])))
#         self.assertTrue(sad[2][2][0] == 'b')

#     def test_parse_criteria(self):

#         # Checking parse returns what we would expect
#         pars = self.pat4.parse_criteria({'spp_code': 'species', 'count': 'count',
#         'x': 1})
#         self.assertTrue(pars[1] == 'spp_code')
#         self.assertTrue(pars[2] == 'count')

#         # Test that energy, mass and count col are None
#         pars = self.pat4.parse_criteria({'spp_code': 'species',
#                                                 'y': 'whole'})
#         self.assertTrue((pars[2] == None) and (pars[3] == None) and (pars[4] ==
#                         None))

#         # If species is not specified correctly an error is thrown
#         self.assertRaises(ValueError, self.pat3.parse_criteria, {'spp_col'
#                             :'species'})
#         # Make sure if count is not passed, no error is thrown
#         self.pat3.parse_criteria({'spp_code': 'species'})

#         # Check energy and mass returns
#         pars = self.pat5.parse_criteria({'spp_code': 'species', 'count':
#         'count', 'energy': 'energy'})

#         self.assertTrue(pars[3] == 'energy')
#         self.assertTrue(pars[4] == None)

#         # Check that combinations in empty dict if no criteria given
#         pars = self.pat5.parse_criteria({'spp_code': 'species', 'count':
#                                 'count'})
#         self.assertTrue(pars[5] == [{}])

#         # TODO: Test that error is thrown if step < prec

#     def test_sar(self):

#         # Checking that sar function returns correct S0 for full plot
#         sar = self.pat3.sar(('x', 'y'), [(1,1)], {'spp_code': 'species',
#         'count': 'count'})
#         self.assertTrue(sar[0]['items'][0] == 5)

#         # Checking for correct result for sar
#         sar = self.pat3.sar(('x', 'y'), [(1,1), (2,2)], {'spp_code': 'species',
#         'count': 'count'})
#         self.assertTrue(np.array_equal(sar[1][1], np.array([3,3,2,3])))
#         sar = self.pat4.sar(('x', 'y'), [(1,1), (1,2), (3,2)], {'spp_code':
#                 'species', 'count': 'count'}, form='sar')
#         self.assertTrue(np.array_equal(sar[1][2], np.array([3,3,2,2,3,1])))

#         # Checking for correct result for ear
#         ear = self.pat3.sar(('x', 'y'), [(1,1), (2,2)], {'spp_code': 'species',
#         'count': 'count'}, form='ear')
#         self.assertTrue(np.array_equal(ear[1][1], np.array([0,1,0,0])))

#         # Test that returned areas are correct
#         sar = self.pat1.sar(('x', 'y'), [(1,1)], {'spp_code': 'species',
#                             'count': 'count'})
#         self.assertTrue(np.round(sar[0]['area'][0], decimals=2) == 0.06)
#         self.assertTrue(sar[0]['items'][0] == 2)

#     def test_universal_sar(self):

#         # Check that it returns the right length
#         criteria = {'spp_code': 'species', 'count' : 'count'}
#         div_cols = ('x', 'y')
#         vals = self.pat8.universal_sar(div_cols, [(1,1), (1,2), (2,2), (2,4),
#                                             (4,4)], criteria)
#         self.assertTrue(len(vals) == 3)

#         # If (1,1) is not passed in it should have a length of zero
#         vals = self.pat8.universal_sar(div_cols, [(1,2), (2,2)], criteria)
#         self.assertTrue(len(vals) == 0)

#         # If (1,1) is not passed in but include_full == True should have len
#         # equal to 1
#         vals = self.pat8.universal_sar(div_cols, [(1,2), (2,2), (2,4)],
#                                        criteria,
#                                                             include_full=True)
#         self.assertTrue(len(vals) == 2)

#         # Test that I get the correct z-value back
#         vals = self.pat8.universal_sar(div_cols, [(1,1), (1,2), (2,2)],
#                                                                     criteria)
#         self.assertTrue(np.round(vals['z'][0], decimals=4) == 0.3390)

#         # If I pass in something other than a halving I should still get
#         # something back
#         vals = self.pat8.universal_sar(div_cols, [(1,1), (2,2), (2,4), (4,4)],
#                                                                     criteria)
#         self.assertTrue(len(vals) == 2)

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

#     def test_ssad(self):

#         # Check that ssad does not lose any individuals
#         ssad = self.pat2.ssad({'spp_code': 'species', 'count': 'count'})
#         sad = self.pat2.sad({'spp_code': 'species', 'count': 'count'})
#         sum_ssad = np.array([sum(val) for val in ssad[1].itervalues()])
#         self.assertTrue(sum(sad[0][1]) == sum(sum_ssad))

#         ssad = self.pat6.ssad({'spp_code': 'species', 'count': 'count'})
#         sad = self.pat6.sad({'spp_code': 'species', 'count': 'count'})
#         sum_ssad = np.array([sum(val) for val in ssad[1].itervalues()])
#         self.assertTrue(sum(sad[0][1]) == sum(sum_ssad))

#         # Manual checks of correct ssad
#         ssad = self.pat2.ssad({'spp_code': 'species', 'count': 'count', 'x':
#                                             2, 'y': 2})
#         self.assertTrue(set(ssad[1]['a']) == {1, 0, 1, 0})
#         self.assertTrue(set(ssad[1]['b']) == {1, 4, 0, 1})
#         self.assertTrue(set(ssad[1]['c']) == {0, 0, 3, 3})
#         self.assertTrue(set(ssad[1]['d']) == {3, 1, 1, 1})

#         ssad = self.pat6.ssad({'spp_code': 'species', 'count': 'count', 'x' :
#                                             2, 'y': 2})
#         self.assertTrue(set(ssad[1]['a']) == {1, 0, 1, 0})
#         self.assertTrue(set(ssad[1]['b']) == {1, 4, 0, 1})
#         self.assertTrue(set(ssad[1]['c']) == {0, 0, 3, 3})
#         self.assertTrue(set(ssad[1]['d']) == {3, 1, 1, 1})

#     def test_ied(self):

#         # Test correct length of result
#         eng = self.pat5.ied({'spp_code': 'species', 'count': 'count',
#         'energy': 'energy'})
#         self.assertTrue(len(eng[0][1]) == 6)

#         # Test error if energy column is missing
#         self.assertRaises(ValueError, self.pat5.ied,
#                                 {'spp_code': 'species', 'count': 'count'})

#         # Test normalize is working
#         eng = self.pat5.ied({'spp_code': 'species', 'count': 'count',
#                         'energy': 'energy', 'x': 2})
#         self.assertTrue(np.array_equal(eng[1][1], np.array([1])))
#         self.assertTrue(len(eng[0][1]) == 5)

#         # Test mass column and normalize
#         eng = self.pat5.ied({'spp_code': 'species', 'count': 'count',
#                         'mass' : 'mass'}, exponent=1, normalize=False)
#         self.assertTrue(np.array_equal(eng[0][1], np.array([17,17,12,23,45,
#                                     110])))

#         # Test that energy overrides mass
#         eng = self.pat5.ied({'spp_code': 'species', 'count': 'count',
#                         'mass' : 'mass', 'energy' : 'energy'}, normalize=False)
#         self.assertTrue(np.array_equal(eng[0][1], np.array([.5,.5,2,3,4,5])))

#     def test_sed(self):

#         # Check correct result
#         eng = self.pat5.sed({'spp_code': 'species', 'count': 'count',
#                                         'energy': 'energy'})
#         self.assertTrue(np.array_equal(eng[0][1]['grt'],
#                                                     np.array([1,1,4,6])))
#         self.assertTrue(np.array_equal(eng[0][1]['rty'],
#                                                     np.array([8,10])))

#         eng = self.pat5.sed({'spp_code': 'species', 'count': 'count',
#                         'energy': 'energy', 'x': 2})
#         self.assertTrue(np.array_equal(eng[1][1]['rty'], np.array([1])))
#         self.assertTrue(len(eng[1][1]) == 2)

# if __name__ == "__main__":
#     unittest.main()
