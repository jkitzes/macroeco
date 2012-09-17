'''
Unit tests for empirical.py
'''

from __future__ import division
import unittest
import os
gcwd = os.getcwd
pd = os.path.dirname
jp = os.path.join
from macroeco.empirical import *
import numpy as np
import random


class TestPatch(unittest.TestCase):

    def setUp(self):
        self.xyfile5 = open('xyfile5.csv','w')
        self.xyfile5.write('''spp_code, x, y, count
grt, .1, .1, 2
grt, .1, .2, 1
grt, .1, .3, 1
rty, .1, .2, 1
rty, .2, .3, 1''')
        self.xyfile5.close()
        self.xymeta5 = {('x', 'maximum'): .2, ('x', 'minimum'): .1, ('x',
        'precision'): .1, ('x', 'type'): 'interval', ('y', 'maximum'): .3,
        ('y', 'minimum'): .1, ('y', 'precision'): .1, ('y', 'type'): 'interval',
        ('spp_code', 'maximum'): None, ('spp_code', 'minimum'): None,
        ('spp_code', 'precision'): None, ('spp_code', 'type'): 'ordinal',
        ('count', 'maximum'): None, ('count', 'minimum'): None, ('count',
        'precision'): None, ('count', 'type'): 'ratio'}

        self.pat1 = Patch('xyfile5.csv')
        # Line below sets metadata manually-no metadata file loaded
        self.pat1.data_table.meta = self.xymeta5 

        self.xyfile6 = open('xyfile6.csv', 'w')
        self.xyfile6.write('''spp_code, x, y, count
a, 0, 0, 1
b, 0, 0, 1
c, 0, 0, 0
d, 0, 0, 3
a, 0, 1, 0
b, 0, 1, 4
c, 0, 1, 0
d, 0, 1, 1
a, 1, 0, 1
b, 1, 0, 0
c, 1, 0, 3
d, 1, 0, 1
a, 1, 1, 0
b, 1, 1, 1
c, 1, 1, 3
d, 1, 1, 1''')
        self.xyfile6.close()
        self.xymeta6 = {('x', 'maximum'): 1, ('x', 'minimum'): 0, ('x',
        'precision'): 1, ('x', 'type'): 'interval', ('y', 'maximum'): 1,
        ('y', 'minimum'): 0, ('y', 'precision'): 1, ('y', 'type'): 'interval',
        ('spp_code', 'maximum'): None, ('spp_code', 'minimum'): None,
        ('spp_code', 'precision'): None, ('spp_code', 'type'): 'ordinal',
        ('count', 'maximum'): None, ('count', 'minimum'): None, ('count',
        'precision'): None, ('count', 'type'): 'ratio'}
        self.pat2 = Patch('xyfile6.csv')
        self.pat2.data_table.meta = self.xymeta6

        self.xyfile7 = open('xyfile7.csv', 'w')
        self.xyfile7.write('''spp_code, x, y, count
tery, 1, 1, 1
1, 1, 1, 1
2, 1, 1, 0
3, 1, 1, 3
0, 1, 2, 0
1, 1, 2, 4
2, 1, 2, 0
tery, 1, 2, 1
0, 2, 1, 1
1, 2, 1, 0
2, 2, 1, 3
3, 2, 1, 1
tery, 2, 2, 0
1, 2, 2, 1
2, 2, 2, 3
3, 2, 2, 1''')
        self.xyfile7.close()
        self.xymeta7 = {('x', 'maximum'): 2, ('x', 'minimum'): 1, ('x',
        'precision'): 1, ('x', 'type'): 'interval', ('y', 'maximum'): 2,
        ('y', 'minimum'): 1, ('y', 'precision'): 1, ('y', 'type'): 'interval',
        ('spp_code', 'maximum'): None, ('spp_code', 'minimum'): None,
        ('spp_code', 'precision'): None, ('spp_code', 'type'): 'ordinal',
        ('count', 'maximum'): None, ('count', 'minimum'): None, ('count',
        'precision'): None, ('count', 'type'): 'ratio'}
        self.pat3 = Patch('xyfile7.csv')
        self.pat3.data_table.meta = self.xymeta7

        self.xyfile8 = open('xyfile8.csv', 'w')
        self.xyfile8.write('''spp_code, x, y, count
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
3, 1, 1, 1
0, 2, 0, 0
1, 2, 0, 0
2, 2, 0, 2
3, 2, 0, 4
0, 2, 1, 0
1, 2, 1, 0
2, 2, 1, 0
3, 2, 1, 1''')
        self.xyfile8.close()
        self.xymeta8 = {('x', 'maximum'): 2, ('x', 'minimum'): 0, ('x',
        'precision'): 1, ('x', 'type'): 'interval', ('y', 'maximum'): 1,
        ('y', 'minimum'): 0, ('y', 'precision'): 1, ('y', 'type'): 'interval',
        ('spp_code', 'maximum'): None, ('spp_code', 'minimum'): None,
        ('spp_code', 'precision'): None, ('spp_code', 'type'): 'ordinal',
        ('count', 'maximum'): None, ('count', 'minimum'): None, ('count',
        'precision'): None, ('count', 'type'): 'ratio'}
        self.pat4 = Patch('xyfile8.csv')
        self.pat4.data_table.meta = self.xymeta8
        self.xyfile9 = open('xyfile9.csv','w')
        self.xyfile9.write('''spp_code, x, y, count, energy
grt, .1, .1, 2, 1
grt, .1, .2, 1, 2
grt, .1, .3, 1, 3
rty, .1, .2, 1, 4
rty, .2, .3, 1, 5''')
        self.xyfile9.close()
        self.xymeta9 = {('x', 'maximum'): .2, ('x', 'minimum'): .1, ('x',
        'precision'): .1, ('x', 'type'): 'interval', ('y', 'maximum'): .3,
        ('y', 'minimum'): .1, ('y', 'precision'): .1, ('y', 'type'): 'interval',
        ('spp_code', 'maximum'): None, ('spp_code', 'minimum'): None,
        ('spp_code', 'precision'): None, ('spp_code', 'type'): 'ordinal',
        ('count', 'maximum'): None, ('count', 'minimum'): None, ('count',
        'precision'): None, ('count', 'type'): 'ratio'}

        self.pat5 = Patch('xyfile9.csv')
        self.pat5.data_table.meta = self.xymeta9 
        self.xyfile10 = open('xyfile10.csv', 'w')
        self.xyfile10.write('''spp_code, x, y, count
a, 0, 0, 1
b, 0, 0, 1
d, 0, 0, 3
b, 0, 1, 4
d, 0, 1, 1
a, 1, 0, 1
c, 1, 0, 3
d, 1, 0, 1
b, 1, 1, 1
c, 1, 1, 3
d, 1, 1, 1''')
        self.xyfile10.close()
        self.xymeta10 = {('x', 'maximum'): 1, ('x', 'minimum'): 0, ('x',
        'precision'): 1, ('x', 'type'): 'interval', ('y', 'maximum'): 1,
        ('y', 'minimum'): 0, ('y', 'precision'): 1, ('y', 'type'): 'interval',
        ('spp_code', 'maximum'): None, ('spp_code', 'minimum'): None,
        ('spp_code', 'precision'): None, ('spp_code', 'type'): 'ordinal',
        ('count', 'maximum'): None, ('count', 'minimum'): None, ('count',
        'precision'): None, ('count', 'type'): 'ratio'}
        self.pat6 = Patch('xyfile10.csv')
        self.pat6.data_table.meta = self.xymeta10



    def tearDown(self):
        os.remove('xyfile5.csv')
        os.remove('xyfile6.csv')
        os.remove('xyfile7.csv')
        os.remove('xyfile8.csv')
        os.remove('xyfile9.csv')
        os.remove('xyfile10.csv')

    #
    # init and set_attributes
    #

    def test_patch_init(self):

        # Test entire table is loaded
        self.assertTrue(len(self.pat1.data_table.table) == 5)
        self.assertTrue(len(self.pat2.data_table.table) == 16)

        # Test that subsetting works
        pat = Patch('xyfile6.csv', {'spp_code': ("!='a'", "!='b'", "!='c'")})
        self.assertTrue(np.all(pat.data_table.table['spp_code'] == 'd'))
        pat = Patch('xyfile7.csv', {'spp_code': "=='tery'"})
        self.assertTrue(sum(pat.data_table.table['count']) == 2)

        # Testing that metadata was set correctly
        self.assertTrue(self.pat1.data_table.meta[('x', 'maximum')] == .2)
        
    def test_sad(self):
        # TODO: Test 'split'
        
        # Test correct result with 'whole' and one division
        sad = self.pat1.sad({'spp_code': 'species', 'count': 'count', 
                                                                    'x': 1})
        self.assertTrue(np.array_equal(sad[1][0][1], np.array([4,2])))
        sad = self.pat1.sad({'spp_code': 'species', 'count': 'count', 
                                                        'x': 'whole'})
        self.assertTrue(np.array_equal(sad[1][0][1], np.array([4,2])))
        sad = self.pat4.sad({'spp_code': 'species', 'count' :'count', 'x': 1})
        self.assertTrue(np.array_equal(sad[0], np.array([0,1,2,3])))

        # Test correct result with other divisions
        sad = self.pat4.sad({'spp_code': 'species', 'count': 'count', 'x': 3,
        'y': 2})
        self.assertTrue(np.array_equal(sad[1][-1][1], np.array([0,0,0,1])))

        # Test that 'whole' and ignore give the same result
        sad1 = self.pat4.sad({'spp_code': 'species', 'count': 'count'})
        sad2 = self.pat4.sad({'spp_code': 'species', 'count': 'count', 'x' :
        'whole'})
        self.assertTrue(np.array_equal(sad1[1][0][1], sad2[1][0][1]))


    def test_parse_criteria(self):

        # Checking parse returns what we would expect 
        pars = self.pat4.parse_criteria({'spp_code': 'species', 'count': 'count',
        'x': 1})
        self.assertTrue(pars[1] == 'spp_code')
        self.assertTrue(pars[2] == 'count')

        # TODO: Remove energy return
        pars = self.pat4.parse_criteria({'spp_code': 'species', 
                                                'y': 'whole'})
        self.assertTrue((pars[2] == None) and (pars[3] == None))

        # If species is not specified correctly an error is thrown
        self.assertRaises(ValueError, self.pat3.parse_criteria, {'spp_col'
                            :'species'})
        # Make sure if count is not passed, no error is thrown
        self.pat3.parse_criteria({'spp_code': 'species'})

        # Check energy and mass returns 
        pars = self.pat5.parse_criteria({'spp_code': 'species', 'count':
        'count', 'energy': 'energy'})

        self.assertTrue(pars[3] == 'energy')
        self.assertTrue(pars[4] == None)

        # Check that combinations in empty dict if no criteria given
        pars = self.pat5.parse_criteria({'spp_code': 'species', 'count':
                                'count'})
        self.assertTrue(pars[5] == [{}])



        # TODO: Write a test where we parse on a string
        # TODO: Test that error is thrown if step < prec

    def test_sar(self):
        
        # Checking that sar function returns correct S0 for full plot
        sar = self.pat3.sar(('x', 'y'), [(1,1)], {'spp_code': 'species',
        'count': 'count'})
        self.assertTrue(sar[1][0] == 5)

        # Checking for correct result for sar
        sar = self.pat3.sar(('x', 'y'), [(1,1), (2,2)], {'spp_code': 'species',
        'count': 'count'})
        self.assertTrue(np.array_equal(sar[2][1], np.array([3,3,2,3])))
        sar = self.pat4.sar(('x', 'y'), [(1,1), (1,2), (3,2)], {'spp_code':
                'species', 'count': 'count'}, form='sar')
        self.assertTrue(np.array_equal(sar[2][2], np.array([3,3,2,2,3,1])))

        # Checking for correct result for ear
        ear = self.pat3.sar(('x', 'y'), [(1,1), (2,2)], {'spp_code': 'species',
        'count': 'count'}, form='ear')
        self.assertTrue(np.array_equal(ear[2][1], np.array([0,1,0,0])))
        
        # Test that returned areas are correct
        sar = self.pat1.sar(('x', 'y'), [(1,1)], {'spp_code': 'species',
                            'count': 'count'})
        self.assertTrue(np.round(sar[0][0], decimals=2) == 0.06)
        self.assertTrue(sar[1][0] == 2)

    def test_ssad(self):
        
        # Check that ssad does not lose any individuals
        ssad = self.pat2.ssad({'spp_code': 'species', 'count': 'count'})
        sad = self.pat2.sad({'spp_code': 'species', 'count': 'count'})
        sum_ssad = np.array([sum(val) for val in ssad[1].itervalues()])
        self.assertTrue(sum(sad[1][0][1]) == sum(sum_ssad))
        
        ssad = self.pat6.ssad({'spp_code': 'species', 'count': 'count'})
        sad = self.pat6.sad({'spp_code': 'species', 'count': 'count'})
        sum_ssad = np.array([sum(val) for val in ssad[1].itervalues()])
        self.assertTrue(sum(sad[1][0][1]) == sum(sum_ssad))

        # Manual checks of correct ssad
        ssad = self.pat2.ssad({'spp_code': 'species', 'count': 'count', 'x':
                                            2, 'y': 2})
        self.assertTrue(set(ssad[1]['a']) == {1, 0, 1, 0})
        self.assertTrue(set(ssad[1]['b']) == {1, 4, 0, 1})
        self.assertTrue(set(ssad[1]['c']) == {0, 0, 3, 3})
        self.assertTrue(set(ssad[1]['d']) == {3, 1, 1, 1})

        ssad = self.pat6.ssad({'spp_code': 'species', 'count': 'count', 'x' :
                                            2, 'y': 2})
        self.assertTrue(set(ssad[1]['a']) == {1, 0, 1, 0})
        self.assertTrue(set(ssad[1]['b']) == {1, 4, 0, 1})
        self.assertTrue(set(ssad[1]['c']) == {0, 0, 3, 3})
        self.assertTrue(set(ssad[1]['d']) == {3, 1, 1, 1})
    
    # TODO: Test mass columns
    def test_ied(self):
        
        # Test correct length of result
        eng = self.pat5.ied({'spp_code': 'species', 'count': 'count',
        'energy': 'energy'})
        self.assertTrue(len(eng[0][1]) == 6)

        # Test error if energy column is missing
        self.assertRaises(ValueError, self.pat5.ied, 
                                {'spp_code': 'species', 'count': 'count'})

        # Test normalize is working
        eng = self.pat5.ied({'spp_code': 'species', 'count': 'count',
                        'energy': 'energy', 'x': 2}) 
        self.assertTrue(np.array_equal(eng[1][1], np.array([1])))
        self.assertTrue(len(eng[0][1]) == 5)

        # TODO: Test correct result without normalize

    def test_sed(self):

        # Check correct result
        eng = self.pat5.sed({'spp_code': 'species', 'count': 'count',
                                        'energy': 'energy'})
        self.assertTrue(np.array_equal(eng[0][1]['grt'],
                                                    np.array([1,1,4,6])))
        self.assertTrue(np.array_equal(eng[0][1]['rty'],
                                                    np.array([8,10])))

        eng = self.pat5.sed({'spp_code': 'species', 'count': 'count',
                        'energy': 'energy', 'x': 2})
        self.assertTrue(np.array_equal(eng[1][1]['rty'], np.array([1])))
        self.assertTrue(len(eng[1][1]) == 2)

