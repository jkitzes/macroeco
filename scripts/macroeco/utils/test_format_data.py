#!/usr/bin/python

'''Testing the classes in format_data.py'''

import unittest
import numpy as np
import format_data as form
import os
from matplotlib.mlab import csv2rec
import glob
gcwd = os.getcwd #get current directory
pd = os.path.dirname #get parent directory
chdir = os.chdir #change directories
jp = os.path.join 

class TestFormatData(unittest.TestCase):
    '''Tests the classes within format_data.py'''

    def setUp(self):

        self.grid1 = open('grid1.csv', 'w')
        self.grid1.write('''Harry-1+Joshua - 6+, hg-4+ty -  67,
                            Harry-3+Joshua-1+y-34+ty - 87, hg-23''')
        self.grid1.close()
        
        # Heavily nested names and blank cell
        self.grid2 = open('grid2.csv', 'w')
        self.grid2.write('''aaaa   - 5&aaaa - 4  &  aaaa - 3, aa - 2&a - 5,
                            aaa - 4& aaaa- 3& a - 1, ''')
        self.grid2.close()
        
        # Grid to be cut
        self.grid3 = open('grid3.csv', 'w')
        self.grid3.write('''aaaa   - 5*&aaaa - 4*  &  aa*aa - *3*$please, aa* -2*&a - 5will I be cut 7658?,
                            aa*a -* 4*& aa*aa- 3*& a* - 1*%maybe, **''')
        self.grid3.close()

        self.dense1 = open('dense1.csv', 'w')
        self.dense1.write('''column, row, fry, the, eggs, well, please
                                0,0,1,2,3,4,5
                                0,1,0,0,,0,23
                                1,0,,,5,45,0
                                1,1,1,1,1,1,1''')
        self.dense1.close()

        self.dense2 = open('dense2.csv', 'w')
        self.dense2.write('''column, row, fry, the, eggs, well, please
                                0,0,1,2,3,4,5
                                0,1,0,0,NA,0,23
                                1,0,NA,NA,5,45,0
                                1,1,1,1,1,1,1''')
        self.dense2.close()


    def tearDown(self):
        os.remove('grid1.csv')
        os.remove('grid2.csv')
        os.remove('grid3.csv')
        os.remove('dense1.csv')
        os.remove('dense2.csv')

    def test_Grid_Data(self):
        grid = form.Grid_Data('grid1.csv', 2, spp_sep='+')
        grid.find_unique_spp_in_grid(spacer='-', spp_sep='+')

        # Does it find the right species?
        spp_list = np.array(['Harry', 'Joshua', 'hg', 'ty', 'y'])
        unq_spp = grid.unq_spp_lists[0]
        self.assertTrue(np.all(spp_list == unq_spp))

        # If I don't truncate '+', it still finds the right species
        grid = form.Grid_Data('grid1.csv', 2)
        grid.find_unique_spp_in_grid(spacer='-', spp_sep='+')

        spp_list = np.array(['Harry', 'Joshua', 'hg', 'ty', 'y'])
        unq_spp = grid.unq_spp_lists[0]
        self.assertTrue(np.all(spp_list == unq_spp))

        # Test that the Dense plot is made correctly
        grid = form.Grid_Data('grid1.csv', 2, spp_sep='+')
        grid.grid_to_dense(spacer='-', spp_sep='+')
        columns = ('cell', 'row',  'column', 'Harry', 'Joshua', 'hg', 'ty',
        'y')
        test_names = grid.Dense_Object.dense_data[0].dtype.names
        self.assertTrue(np.all(test_names == columns))

        # Test that values are correct
        dense_obj = grid.Dense_Object
        pred = np.array([0,0,4,23])
        test = dense_obj.dense_data[0]['hg']
        self.assertTrue(np.all(pred == test))
        pred = np.array([1,3,0,0])
        test = dense_obj.dense_data[0]['Harry']
        self.assertTrue(np.all(pred == test))
        pred = np.array([6,1,0,0])
        test = dense_obj.dense_data[0]['Joshua']
        self.assertTrue(np.all(pred == test))
        pred = np.array([0,34,0,0])
        test = dense_obj.dense_data[0]['y']
        self.assertTrue(np.all(pred == test))
        pred = np.array([0,87,67,0])
        test = dense_obj.dense_data[0]['ty']
        self.assertTrue(np.all(pred == test))

        # Tested heavy nesting and empty cell
        grid = form.Grid_Data('grid2.csv', 2)
        grid.find_unique_spp_in_grid(spacer='-', spp_sep='&')
        unq_spp = np.array(['a', 'aa', 'aaa', 'aaaa']) 
        pred = grid.unq_spp_lists[0]
        self.assertTrue(np.all(unq_spp == pred))
        
        grid.grid_to_dense(spacer='-', spp_sep='&')
        dense_obj = grid.Dense_Object
        pred = np.array([0,1,5, 0])
        test = dense_obj.dense_data[0]['a']
        self.assertTrue(np.all(pred == test))
        pred = np.array([0,0,2, 0])
        test = dense_obj.dense_data[0]['aa']
        self.assertTrue(np.all(pred == test))
        pred = np.array([0,4,0, 0])
        test = dense_obj.dense_data[0]['aaa']
        self.assertTrue(np.all(pred == test))
        pred = np.array([12,3,0, 0])
        test = dense_obj.dense_data[0]['aaaa']
        self.assertTrue(np.all(pred == test))

        # Testing remove, replace, and truncation functions
        grid = form.Grid_Data('grid3.csv', [2], spp_sep='&')
        grid.truncate_grid_cells(['$pl', 'will', '%may'])
        grid.remove_and_replace('*', '')

        grid.find_unique_spp_in_grid(spacer='-', spp_sep='&')
        unq_spp = np.array(['a', 'aa', 'aaa', 'aaaa']) 
        pred = grid.unq_spp_lists[0]
        self.assertTrue(np.all(unq_spp == pred))
        
        grid.grid_to_dense(spacer='-', spp_sep='&')
        dense_obj = grid.Dense_Object
        pred = np.array([0,1,5, 0])
        test = dense_obj.dense_data[0]['a']
        self.assertTrue(np.all(pred == test))
        pred = np.array([0,0,2, 0])
        test = dense_obj.dense_data[0]['aa']
        self.assertTrue(np.all(pred == test))
        pred = np.array([0,4,0, 0])
        test = dense_obj.dense_data[0]['aaa']
        self.assertTrue(np.all(pred == test))
        pred = np.array([12,3,0, 0])
        test = dense_obj.dense_data[0]['aaaa']
        self.assertTrue(np.all(pred == test))

        # Testing reset to archival
        grid.reset_grid_data()
        temp_str = 'aaaa-5*&aaaa-4*&aa*aa-*3*$please'
        data_str = grid.grid_data[0]['0'][0]
        self.assertTrue(temp_str == data_str)


        # Test that multiple data sets work
        self.assertRaises(AssertionError, form.Grid_Data, \
                                    glob.glob('grid*.csv'), [2,2])

        grid = form.Grid_Data(glob.glob('grid*.csv'), 2, archival=False)
        
        # reset_archival should fail in this case
        self.assertRaises(ValueError, grid.reset_grid_data)

        # All the truncation should make the the two data sets equal
        grid.truncate_grid_cells(['$pl', 'will', '%may'])
        grid.remove_and_replace('*', '')
        for col in xrange(grid.cols[0]):
            for row in xrange(len(grid.grid_data[0])):
                self.assertTrue(grid.grid_data[1][col][row] ==\
                grid.grid_data[2][col][row])


    def test_Dense_Data(self):
        
        # Test that the expected values are read in
        dense = form.Dense_Data('dense1.csv', replace=('', 0))
        spp_arr = np.array([1,0,0,1])
        read_in = dense.dense_data[0]['fry']
        self.assertTrue(np.all(spp_arr == read_in))
        spp_arr = np.array([2,0,0,1])
        read_in = dense.dense_data[0]['the']
        self.assertTrue(np.all(spp_arr == read_in))
        spp_arr = np.array([3,0,5,1])
        read_in = dense.dense_data[0]['eggs']
        self.assertTrue(np.all(spp_arr == read_in))
        spp_arr = np.array([4,0,45,1])
        read_in = dense.dense_data[0]['well']
        self.assertTrue(np.all(spp_arr == read_in))
        spp_arr = np.array([5,23,0,1])
        read_in = dense.dense_data[0]['please']
        self.assertTrue(np.all(spp_arr == read_in))

        dense = form.Dense_Data('dense1.csv', replace=('NA', 0))
        spp_arr = np.array([1,0,0,1])
        read_in = dense.dense_data[0]['fry']
        self.assertTrue(np.all(spp_arr == read_in))
        spp_arr = np.array([2,0,0,1])
        read_in = dense.dense_data[0]['the']
        self.assertTrue(np.all(spp_arr == read_in))
        spp_arr = np.array([3,0,5,1])
        read_in = dense.dense_data[0]['eggs']
        self.assertTrue(np.all(spp_arr == read_in))
        spp_arr = np.array([4,0,45,1])
        read_in = dense.dense_data[0]['well']
        self.assertTrue(np.all(spp_arr == read_in))
        spp_arr = np.array([5,23,0,1])
        read_in = dense.dense_data[0]['please']
        self.assertTrue(np.all(spp_arr == read_in))

        # Test Dense to Columnar





        
        


        
        











if __name__ == '__main__':
    unittest.main()

        






        
