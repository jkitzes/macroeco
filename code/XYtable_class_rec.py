'''
Python program to examine SAR and EAR
'''

import numpy as np
import random
import matplotlib as plt

# CONVERSION UTILITIES
def xytable_to_sparse(xytable):
    ''' Convert xytable data into sparse plot data '''
    # TODO: Loop through rows, if unique row add to sparse w/ a 1 in last col, 
    # if matches a row already in sparse, increment last col by 1.
    pass


def dense_to_sparse(dense, unit):
    ''' Convert dense plot data into sparse plot data '''
    sp_y = dense.shape[0]
    sp_x = dense.shape[1]
    nspp = dense.shape[2]

    sparse = []
    for spp in xrange(0, nspp):
        for x in xrange(0, sp_x):
            for y in xrange(0, sp_y):
                if dense[y, x, spp] > 0:
                    sparse.append([spp, x * unit, y * unit, dense[y, x, spp]])

    return np.array(sparse)


def sparse_to_dense(sparse, x_minmax, y_minmax, unit, nspp):
    ''' Convert dense plot data into sparse plot data '''
    x_min = x_minmax[0]
    y_min = y_minmax[0]
    nx = (x_minmax[1] - x_minmax[0] + unit) / float(unit)
    ny = (y_minmax[1] - y_minmax[0] + unit) / float(unit)
    
    dense = np.zeros((ny, nx, nspp))
    for row in sparse:
        dense[(row[2]-y_min)/unit, (row[1]-x_min)/unit, row[0]] += row[3]
    
    return dense


# Temporarily define sample sparse and dense plots for testing
test_dense = np.array((0,0,0,1,1,3,0,0,4,1,0,0,0,2,2,1)*2).reshape(4,4,2)
test_sparse = dense_to_sparse(test_dense, 0.1)


# LOAD_PLOT - Load plot data from file
def load_plot(filename):
    '''
    Load plot and metadata from file. Data should be returned as a 3D numpy 
    array if dense, or a 4 col numpy array if sparse. If sparse, first col must 
    be an int index of species number, and function should also load the 
    recarray that shows the match between index and spp name/code as given in 
    the original data.
    '''
    # TODO: Add loading and parsing metadata
    # TODO: Currently only takes in comma delimited, allow other options
    # TODO: Currently only takes in sparse - may be harder to read grid data in 
    # file - what format will it be?
    try:
        data = np.loadtxt(filename, skiprows = 1, delimiter = ',')
    except IOError as detail:
        print detail
    
    return data


# CLASS PLOT - Define plot class
class Plot:
    '''
    Plot class to store abundance data for multiple species in contiguous 
    subplots. Multiple non-contiguous plots, each of which may contain 
    subplots, should be declared as class Network.
    '''

    def __init__(self, data, x_minmax, y_minmax, nspp, sparse = True,
                 unit = 0.1):
        '''
        x_minmax and y_minmax are tuples giving the lowest and highest possible 
        subplot centroid, respectively, nspp is number of species 
        '''
        self.data = data

        self.x_min = x_minmax[0]
        self.x_max = x_minmax[1]
        self.p_width = self.x_max - self.x_min + unit  # Width in num units

        self.y_min = y_minmax[0]
        self.y_max = y_minmax[1]
        self.p_height = self.y_max - self.y_min + unit  # Height in num units

        self.nspp = nspp
        self.sparse = sparse
        self.unit = unit

        if sparse:
            self.total_abund = self.sparse_abund(data)
        else:
            self.total_abund = self.data.sum(axis=0).sum(axis=0)

        # TODO: Error checking for correct plot type
        # TODO: Change name Plot to something more unique - grid, tract, area?
        # TODO: Convert arguments after data to metadata input from 


    def SAD_sample(self, width, height, full = False):
        '''
        Calculate a sampled SAR with subplots with subplots of specified width 
        and height (given in units).
        '''
        pass



    def SAD_grid(self, div_list, summary = ''):
        '''
        summary SAD, SAR, EAR
        Calculate a gridded SAR from subplots created by evenly dividing a plot 
        along the vertical and horizontal dimensions, with numbers of divisions 
        given by tuples div in div_list, wher div = (x_divs, y_divs). A 
        division of one corresponds to no cut (ie, the whole plot).

        Because of rounding, x_ and y_ divs that do not evenly split the plot 
        in integer number of units will lead to subplots with different numbers 
        of sample points, artificially inflating the variance of the species 
        abundance or counts of the subplots. To avoid this, the values of x_ 
        and y_divs must divide each dimension into an even number of units.
        '''
        # TODO: Error check that div must be >= 1

        out_area = []
        out_result = []

        for div in div_list:  # Loop each division tuple
            sp_width = self.p_width / float(div[0])
            sp_height = self.p_height / float(div[1])
            # TODO: Error check that x and y strips divide dimensions evenly - 
            # use remainder function and ensure 0.

            for x_div_count in xrange(0,div[0]):  # Loop x_divs
                x_st = self.x_min + x_div_count * sp_width
                x_en = x_st + sp_width

                for y_div_count in xrange(0,div[1]):  # Loop y_divs
                    y_st = self.y_min + y_div_count * sp_height
                    y_en = y_st + sp_height

                    out_area.append(sp_width * sp_height)
                    out_result.append(self.subSAD(x_st, x_en, y_st, y_en, 
                                                  summary))

        return np.vstack((out_area, 
                          np.array(out_result).transpose())).transpose()


    def subSAD(self, x_st, x_en, y_st, y_en, summary):
        '''
        Calculates a SAD for a subplot of known starting and ending 
        coordinates. Only works for rectangles. If full, returns a ndarray of 
        all species with abundances of each, else returns an int count of the 
        species in the subplot.
        '''
        
        if self.sparse:
            in_sp = np.all([self.data[:,1] >= x_st, self.data[:,1] < x_en,
                           self.data[:,2] >= y_st, self.data[:,2] < y_en],
                           axis = 0)
            sp_data = self.data[in_sp]
            sub_abund = self.sparse_abund(sp_data)
                
        else:
            sub_abund = self.data[y_st:y_en, x_st:x_en, 
                                  :].sum(axis=0).sum(axis=0)

        if summary is 'SAR':
            return sum(sub_abund > 0)
        elif summary is 'EAR':
            return sum(sub_abund == self.total_abund)
        else:
            return sub_abund


    def sparse_abund(self, data):
        ''' Calculate abundance of each species in sparse plot data '''
        abund = np.zeros(self.nspp)
        for row in data:
            abund[row[0]] += row[3]
        return abund


# CLASS NETWORK - Network of plots
class Network():
    pass

