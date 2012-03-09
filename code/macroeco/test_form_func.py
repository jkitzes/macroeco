#Testing form_func.py



import unittest
import os
from form_func import *
import shutil
import numpy as np

gcwd = os.getcwd
jp = os.path.join
pd = os.path.dirname

class TestFormFunc(unittest.TestCase):
    '''Tests the functions with in form_func.py'''

    def setUp(self):
        self.spp_array1 = np.array(['as', 'as', 'as', 'as', 'as'])
        self.spp_array2 = np.array([2,2,3,5,6,3,4,5,7,8])
        self.spp_array3 = np.array(['as','bn', 'as', 'ty', 'bn'])

    def test_make_spec_dict(self):
        spec_dict = make_spec_dict(self.spp_array1)
        self.assertTrue(len(spec_dict) == 1)
        self.assertTrue(spec_dict['spp_code'][0] == 0)
        self.assertTrue(spec_dict['spp'][0] == 'as')
        spec_dict = make_spec_dict(self.spp_array2)
        self.assertTrue(len(spec_dict) == 7)
        self.assertTrue(sum(np.equal(spec_dict['spp_code'], \
                            np.array([0,1,2,3,4,5,6]))) == 7)
        spec_dict = make_spec_dict(self.spp_array3)
        self.assertTrue(len(spec_dict) == 3)
        
        

if __name__ == '__main__':
    unittest.main()

