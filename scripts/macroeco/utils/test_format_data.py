#!/usr/bin/python

'''Testing the classes in format_data.py'''

import unittest
import numpy as np
import format_data as form
import os
import glob
gcwd = os.getcwd #get current directory
pd = os.path.dirname #get parent directory
chdir = os.chdir #change directories
jp = os.path.join 

class TestFormatData(unittest.TestCase):
    '''Tests the classes within format_data.py'''

    def setUp(self):

        self.grid1 = open('grid1.csv', 'w')
        self.grid1.write('''Harry-1+Joshua - 6+, hg-4+ty -  67,\nHarry-3+Joshua-1+y-34+ty - 87, hg-23''')
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

        self.dense3 = open('dense3.csv', 'w')
        self.dense3.write('''column, row, fry, the, eggs, well, please, j
                                0,0,1,2,3,4,5,2
                                0,1,0,0,NA,0,23,5
                                1,0,NA,NA,5,45,0,6
                                1,1,1,1,1,1,1,7''')
        self.dense3.close()

        self.dense4 = open('dense4.csv', 'w')
        self.dense4.write('''column, row, fry, the, eggs, well, please, j,h
                                0,0,1,2,3,4,5,2,t
                                0,1,0,0,0,0,23,5,u
                                1,0,1,0,5,45,0,6,k
                                1,1,1,1,1,1,1,7,m''')
        self.dense4.close()

        self.trans1 = open('trans1.csv', 'w')
        self.trans1.write(
'''spp, island, tree, b1, b2, b3, b4, b5, nm, fun
h,Marta,1,1,2,3,4,5,j,56
t,Marta,2,1,1,1,1,0,k,78
h,Garry,1,2,3,4,5,6,j,123
t,Garry,2,0,1,2,0,5,u,456''')
        self.trans1.close()
        
        self.col1 = open('col1.csv', 'w')
        self.col1.write('''spp, x, y, dbh1, dbh2, john
l,1,1,34,38,g
y,2,1,100,10,g
h,1,2,1,1,g
y,2,2,300,2,f''')
        self.col1.close()

        self.col2 = open('col2.csv', 'w')
        self.col2.write('''spp, x, y, dbh1, dbh2, john
l,1,,34,38,g
y,2,1,100,10,g
h,,2,1,1,NA
y,2,1,300,2,f''')
        self.col2.close()



    def tearDown(self):
        os.remove('grid1.csv')
        os.remove('grid2.csv')
        os.remove('grid3.csv')
        os.remove('dense1.csv')
        os.remove('dense2.csv')
        os.remove('dense3.csv')
        os.remove('dense4.csv')
        os.remove('trans1.csv')
        os.remove('col1.csv')
        os.remove('col2.csv')

    def test_Grid_Data(self):
        grid = form.Grid_Data('grid1.csv', spp_sep='+')
        grid.find_unique_spp_in_grid(spacer='-', spp_sep='+')

        # Does it find the right species?
        spp_list = np.array(['Harry', 'Joshua', 'hg', 'ty', 'y'])
        unq_spp = grid.unq_spp_lists[0]
        self.assertTrue(np.all(spp_list == unq_spp))

        # If I don't truncate '+', it still finds the right species
        grid = form.Grid_Data('grid1.csv')
        grid.find_unique_spp_in_grid(spacer='-', spp_sep='+')

        spp_list = np.array(['Harry', 'Joshua', 'hg', 'ty', 'y'])
        unq_spp = grid.unq_spp_lists[0]
        self.assertTrue(np.all(spp_list == unq_spp))

        # Test that the Dense plot is made correctly
        grid = form.Grid_Data('grid1.csv', spp_sep='+')
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
        grid = form.Grid_Data('grid3.csv', spp_sep='&')
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

        grid = form.Grid_Data(glob.glob('grid*.csv'), archival=False)
        
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
        
        # NAs should all be turned to 0's
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

        # Test dense_to_columnar
        dense = form.Dense_Data(['dense2.csv', 'dense3.csv'], replace=('NA',0))
        dense.dense_to_columnar(2, (5,6))
        col = dense.Columnar_Object
        col.merge_data()
        unq_spp = np.unique(['eggs', 'fry', 'the', 'well', 'please', 'j'])
        pred_unq_spp = np.unique(col.merged_data['spp'])
        self.assertTrue(np.all(unq_spp == pred_unq_spp))
        count = [1,2,3,4,5]
        self.assertTrue(np.all(count == col.merged_data['count'][:5]))
        count = [1,1,1,1,1,7]        
        self.assertTrue(np.all(count == col.merged_data['count'][-6:]))
        self.assertTrue(len(col.merged_data) == 30)
        self.assertTrue(col.merged_data['count'][5] == 23)
        self.assertTrue(col.merged_data['spp'][5] == 'please')

        # Test correct extension of num_spp
        dense = form.Dense_Data(['dense2.csv', 'dense2.csv'], replace=('NA',0))
        self.assertRaises(TypeError, dense.dense_to_columnar, 2, (5,6,7))
        dense.dense_to_columnar(2, 5)
        col = dense.Columnar_Object
        count = np.array([1,2,3,4,5])
        self.assertTrue(np.all(col.columnar_data[0]['count'][:5] == count))
        self.assertTrue(np.all(col.columnar_data[1]['count'][:5] == count))

        # Test trailing column after species
        dense = form.Dense_Data(['dense4.csv'])
        dense.dense_to_columnar(2, 5)
        col = dense.Columnar_Object
        comp = np.array([2,2,2,2,2,5,6,6,6,7,7,7,7,7])
        self.assertTrue(np.all(comp == col.columnar_data[0]['j']))
        comp = np.array(['t', 't', 't', 't', 't', 'u', 'k', 'k', 'k', 'm', 'm',
                            'm', 'm', 'm'])
        self.assertTrue(np.all(comp == col.columnar_data[0]['h']))

    def test_Transect_Data(self):
        

        # Already tested replace_vals test_Dense_Data
        trans = form.Transect_Data('trans1.csv', replace=('0', 1))
        trans.transect_to_columnar(3, 5)
        col = trans.Columnar_Object
        count = np.array([1,2,3,4,5,1,1,1,1,1,2,3,4,5,6,1,1,2,1,5])
        self.assertTrue(np.all(count == col.columnar_data[0]['count']))

        # Test that transect data reads in correctly and converts to columnar
        trans = form.Transect_Data('trans1.csv')
        trans.transect_to_columnar(3, 5)
        col = trans.Columnar_Object
        count = np.array([1,2,3,4,5,1,1,1,1,2,3,4,5,6,1,2,5])
        self.assertTrue(np.all(count == col.columnar_data[0]['count']))

        # Test multiple datasets are converted correctly
        trans = form.Transect_Data(['trans1.csv', 'trans1.csv'])
        trans.transect_to_columnar(3,5)
        col = trans.Columnar_Object
        col.merge_data()
        self.assertTrue(np.all(np.concatenate((count, count)) ==
                                                    col.merged_data['count']))
    def test_Columnar_Data(self):

        # Testing missing values
        col = form.Columnar_Data('col2.csv', missingd={'y' : '', 'x' : '',
                                            'john' : 'NA'}, delete_missing=True)
        self.assertTrue(len(col.columnar_data[0]) == 2)
        self.assertTrue(np.all(col.columnar_data[0]['spp'] == np.array(['y',
                                                                        'y'])))
        self.assertTrue(np.all(col.columnar_data[0]['dbh1'] == np.array([100,
                                                                        300])))

        # No missing values; Test subsetting
        col = form.Columnar_Data('col1.csv')
        col.subset_data({'john' : ('!=', 'f')})
        self.assertTrue(np.all(col.columnar_data[0]['john'] == np.array(['g',
                                                                   'g', 'g'])))
        # Test reset
        col.reset_columnar_data()
        check = np.array(['g','g','g','f'])
        self.assertTrue(np.all(col.columnar_data[0]['john'] == check))

        # Test splitting
        col.split_up_data_by_field([('dbh1',), ('dbh2',)])
        self.assertTrue(len(col.columnar_data) == 2)
        dbh1 = np.array([34,100,1,300])
        dbh2 = np.array([38,10,1,2])
        try:
            col.columnar_data[0]['dbh2']
        except ValueError:
            pass

        try:
            col.columnar_data[1]['dbh1']
        except ValueError:
            pass

        self.assertTrue(np.all(col.columnar_data[0]['dbh1'] == dbh1))
        self.assertTrue(np.all(col.columnar_data[1]['dbh2'] == dbh2))

        col.reset_columnar_data()

        col.split_up_data_by_field([('spp', 'x'), ('y',)])
        self.assertTrue(len(col.columnar_data) == 2)
        td1 = np.array(['spp', 'x', 'dbh1', 'dbh2', 'john'])
        td2 = np.array(['y', 'dbh1', 'dbh2', 'john'])
        d1 = np.array(col.columnar_data[0].dtype.names)
        d2 = np.array(col.columnar_data[1].dtype.names)
        self.assertTrue(np.all(d1 == td1))
        self.assertTrue(np.all(d2 == td2))
        
        # Test change column names
        col.reset_columnar_data()
        col.split_up_data_by_field([('dbh1',), ('dbh2',)])
        self.assertRaises(ValueError,col.change_column_names, [('x', 'y')], 
                            ['hello'])

        col.reset_columnar_data()
        col.split_up_data_by_field([('dbh1',), ('dbh2',)])
        col.change_column_names([('dbh1', 'dbh2'), ('john',)], ['dbh', 'harry'])
        nms = np.array(['spp', 'x', 'y', 'dbh', 'harry'])
        dtnms1 = np.array(col.columnar_data[0].dtype.names)
        dtnms2 = np.array(col.columnar_data[1].dtype.names)
        self.assertTrue(np.all(nms == dtnms1))
        self.assertTrue(np.all(nms == dtnms2))

        # Test if long names added
        col.reset_columnar_data()
        col.split_up_data_by_field([('dbh1',), ('dbh2',)])
        col.change_column_names([('dbh1', 'dbh2')], ['goofy_chew'])
        nms = np.array(['spp', 'x', 'y', 'goofy_chew', 'john'])
        dtnms1 = np.array(col.columnar_data[0].dtype.names)
        dtnms2 = np.array(col.columnar_data[1].dtype.names)
        self.assertTrue(np.all(nms == dtnms1))
        self.assertTrue(np.all(nms == dtnms2))

        # Test adding fields to data
        col.reset_columnar_data()
        col.split_up_data_by_field([('dbh1',), ('dbh2',)])
        col.change_column_names([('dbh1', 'dbh2')], ['dbh'])
        col.add_fields_to_data_list({'year' : (1998, 2001), 'body' : ('large',
                                                                     'small')})
        year1 = np.repeat('1998', 4)
        year2 = np.repeat('2001', 4)
        body1 = np.repeat('large', 4)
        body2 = np.repeat('small', 4)
        self.assertTrue(np.all(year1 == col.columnar_data[0]['year']))
        self.assertTrue(np.all(year2 == col.columnar_data[1]['year']))
        self.assertTrue(np.all(body1 == col.columnar_data[0]['body']))
        self.assertTrue(np.all(body2 == col.columnar_data[1]['body']))

        # Test adding different dtypes
        col.reset_columnar_data()
        col.split_up_data_by_field([('dbh1',), ('dbh2',)])
        col.change_column_names([('dbh1', 'dbh2')], ['dbh'])
        col.add_fields_to_data_list({'year' : (1998, 2001), 'body' : ('large',
                                            'small')}, descr={'year': np.int, 
                                                              'body': 'S20'})

        year1 = np.repeat(1998, 4)
        year2 = np.repeat(2001, 4)
        body1 = np.repeat('large', 4)
        body2 = np.repeat('small', 4)
        self.assertTrue(np.all(year1 == col.columnar_data[0]['year']))
        self.assertTrue(np.all(year2 == col.columnar_data[1]['year']))
        self.assertTrue(np.all(body1 == col.columnar_data[0]['body']))
        self.assertTrue(np.all(body2 == col.columnar_data[1]['body']))

        # Test remove columns
        col = form.Columnar_Data(['col1.csv', 'col2.csv'],  missingd={'y' : '',
                        'x' : '', 'john' : 'NA'}, delete_missing=True)
        self.assertTrue(len(col.columnar_data[0]) == 4)
        self.assertTrue(len(col.columnar_data[1]) == 2)
        col.remove_columns('john')
        test_nm = set(['x','y', 'spp', 'dbh1', 'dbh2'])
        self.assertTrue(test_nm == set(col.columnar_data[0].dtype.names))
        self.assertTrue(test_nm == set(col.columnar_data[1].dtype.names))

        col.remove_columns(['x', 'y'])
        test_nm = set(['spp', 'dbh1', 'dbh2'])
        self.assertTrue(test_nm == set(col.columnar_data[0].dtype.names))
        self.assertTrue(test_nm == set(col.columnar_data[1].dtype.names))
        
        # Try removing row that is not there, no error is thrown
        col.remove_columns(['x'])
        self.assertTrue(test_nm == set(col.columnar_data[0].dtype.names))
        self.assertTrue(test_nm == set(col.columnar_data[1].dtype.names))

        # Fractionate is tested in test_form_func.py
        col.reset_columnar_data()
        col.fractionate_data((1,1), (.5,.5), ('x', 'y'))
        self.assertTrue(np.all(np.array([0,.5,0,.5]) ==
                                                    col.columnar_data[0]['x']))

        # Test merge data
        col.merge_data()
        self.assertTrue(len(col.merged_data) == 6)
        spp = np.array(['l','y','h','y','y','y'])
        self.assertTrue(np.all(col.merged_data['spp'] == spp))
        dbh2 = np.array([34,100,1,300,100,300])
        self.assertTrue(np.all(col.merged_data['dbh1'] == dbh2))

        # Try to break merge data
        col.columnar_data = [col.merged_data]
        col.merge_data()



if __name__ == '__main__':
    unittest.main() 
