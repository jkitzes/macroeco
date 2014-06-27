from __future__ import division

from numpy.testing import (TestCase, assert_equal, assert_array_equal,
                           assert_almost_equal, assert_array_almost_equal,
                           assert_allclose, assert_, assert_raises)

import numpy as np
from macroeco.misc import *
import pandas as pd

#TODO: Test data_read_write


class TestFormatData(TestCase):

    def test_simple_stack(self):

        # Test that stack gives the same answer as predicted by hand
        test_data = pd.DataFrame({'row': [1, 2, 1, 2],
            'column': [1, 1, 2, 2], 'labelA': [1, 0, 3, 4],
            'labelB': [3, 2, 1, 4]})

        expected = pd.DataFrame({'row': [1,1,2,2,1,1,2,2], 'column':
            [1,1,1,1,2,2,2,2], 'label': np.tile(['labelA', 'labelB'], 4),
            'count': [1,3,0,2,3,1,4,4]}, columns=['row', 'column', 'label',
            'count'])

        stack = format_dense(test_data, ['row', 'column'])
        assert_equal(np.all(stack == expected), True)

    def test_label_count_col(self):
        # Test whether changing label count col work
        test_data = pd.DataFrame({'year': ['02', '03'], 'spp1': [1, 2],
                        'spp2': [3, 4]})

        expected = pd.DataFrame({'year': np.repeat(['02', '03'], 2), 'spp':
            np.tile(['spp1', 'spp2'], 2), 'ind': [1,3,2,4]}, columns=['year',
            'spp', 'ind'])

        stack = format_dense(test_data, ['year'], label_col="spp",
            count_col="ind")

        print stack
        print expected

        assert_equal(np.all(stack == expected), True)

    def test_drop_nan(self):
        # Test whether dropping nan function works

        test_data = pd.DataFrame({'year': ['02', '03'], 'spp1': [1, np.nan],
                        'spp2': [np.nan, 4]})

        expected = pd.DataFrame({'year': ['02', '03'], 'label':
            ['spp1', 'spp2'], 'count': [1,4]}, columns=['year',
            'label', 'count'])

        stack = format_dense(test_data, ['year'], drop_na=True)

        assert_equal(np.all(stack == expected), True)

    def test_nan_to_zero(self):
        # Test whether setting nan to zero function works

        test_data = pd.DataFrame({'year': ['02', '03'], 'spp1': [1, np.nan],
                        'spp2': [np.nan, 4]})

        expected = pd.DataFrame({'year': np.repeat(['02', '03'], 2), 'label':
            np.tile(['spp1', 'spp2'], 2), 'count': [1,0,0,4]}, columns=['year',
            'label', 'count'])

        stack = format_dense(test_data, ['year'], nan_to_zero=True)

        assert_equal(np.all(stack == expected), True)










