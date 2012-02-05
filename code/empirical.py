'''
Python program to examine SAR and EAR
'''

import numpy as np
import numpy.random as rand

# UTILITY METHODS
def sample(array):
    '''
    Returns random element from 1D numpy array.
    Equivalent to np.random.choice in Numpy v1.7
    '''
    return array[rand.randint(0,len(array))]

def distance(pt1, pt2):
    ''' Calculate Euclidean distance between two points '''
    return np.sqrt((pt1[0] - pt2[0])**2 + (pt1[1] - pt2[1])**2)


# CONVERSION UTILITIES
def xytable_to_sparse(xytable):
    ''' Convert xytable data into sparse plot data '''
    sparse_row = [np.array((None, None, None))]  # Dummy row, removed later,
    sparse_count = [None]                        # so '==' always works
    
    for row in xytable:
        if ~np.all(sparse_row == row, axis = 1).any():
            sparse_row.append(row); sparse_count.append(1)
        else:
            row_ind = np.where(np.all(sparse_row == row, axis = 1))[0][0]
            sparse_count[row_ind] += 1

    return np.vstack((np.array(sparse_row[1:]).transpose(), 
                      np.array(sparse_count[1:]))).transpose()


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
test_dense[1,0,0] = 0
test_sparse = dense_to_sparse(test_dense, 0.1)
test_xy = test_sparse[:,0:3]
test_xy = np.vstack((test_xy,test_xy[0:8,:]))


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

        if self.sparse:
            self.total_abund = self.sparse_abund(data)
        else:
            self.total_abund = self.data.sum(axis=0).sum(axis=0)

        # TODO: Error checking for correct plot type
        # TODO: Change name Plot to something more unique - grid, tract, area?
        # TODO: Convert arguments after data to metadata input 


    def SAD_grid(self, div_list, summary = ''):
        '''
        Calculate gridded SADs from subplots created by evenly dividing a plot 
        along the vertical and horizontal dimensions, with numbers of divisions 
        given by tuples div in div_list, wher div = (x_divs, y_divs). A 
        division of one corresponds to no cut (ie, the whole plot).

        summary takes three arguments: '' (full SAD), 'SAR', and 'EAR'

        Because of rounding, x_ and y_ divs that do not evenly split the plot 
        in integer number of units will lead to subplots with different numbers 
        of sample points, artificially inflating the variance of the species 
        abundance or counts of the subplots. To avoid this, the values of x_ 
        and y_divs must divide each dimension into an even number of units.
        '''
        # TODO: Error check that div must be >= 1

        result = []

        for div in div_list:  # Loop each division tuple
            div_result = []
            sp_width = self.p_width / float(div[0])
            sp_height = self.p_height / float(div[1])
            # TODO: Error check that x and y strips divide dimensions evenly - 
            # use remainder function on *_max + unit and ensure 0.

            for x_div_count in xrange(0,div[0]):  # Loop x_divs
                x_st = self.x_min + x_div_count * sp_width
                x_en = x_st + sp_width

                for y_div_count in xrange(0,div[1]):  # Loop y_divs
                    y_st = self.y_min + y_div_count * sp_height
                    y_en = y_st + sp_height

                    div_result.append(self.subSAD(x_st, x_en, y_st, y_en, 
                                                  summary))

            result.append(np.array(div_result))

        return result


    def SAD_sample(self, wh_list, samples, summary = ''):
        '''
        Calculate a sampled SAR with subplots with subplots of specified width 
        and height (given in units).
        '''
        result = []

        for wh in wh_list:  # Loop each width-height tuple
            wh_result = []
            sp_width = wh[0]
            sp_height = wh[1]
            x_origins = np.arange(self.x_min, self.x_max - sp_width + 
                                  2*self.unit, self.unit)
            y_origins = np.arange(self.y_min, self.y_max - sp_width + 
                                  2*self.unit, self.unit)

            for s in xrange(0,samples):  # Loop each sample
                # TODO: Currently fails for sp_width = whole plot
                x_st = sample(x_origins)
                y_st = sample(y_origins)

                x_en = x_st + sp_width
                y_en = y_st + sp_height

                wh_result.append(self.subSAD(x_st, x_en, y_st, y_en, 
                                             summary))

            result.append(np.array(wh_result))

        return result


    def QS_grid(self, div_list):
        '''
        Calculates commonality between pairs of subplots in a gridded plot as 
        the Sorensen index (equivalent to Chi in Harte et al. papers, and 1 - 
        Bray-Curtis dissimilarity).

        In result, the first two cols are the plot indices, indexed in single 
        sequence rowwise.
        '''
        # TODO: Make sure that divs divide plot evenly and that divs are of
        # width at least one unit.

        result = []

        # SAD and sp_cent have the same number of list elements and row
        # structure within arrays
        SAD = self.SAD_grid(div_list, summary = '')
        sp_cent = self.get_sp_centers(div_list)

        for ind_div, div in enumerate(div_list):
            div_result = []
            nsp = div[0] * div[1]

            div_SAD = SAD[ind_div]
            div_sp_cent = sp_cent[ind_div]

            for ind_a in xrange(0, nsp - 1):
                spp_in_a = (div_SAD[ind_a,:] > 0)
                for ind_b in xrange(ind_a + 1, nsp):
                    spp_in_b = (div_SAD[ind_b,:] > 0)
                    dist = distance(div_sp_cent[ind_a,:], div_sp_cent[ind_b,:])
                    QS = sum(spp_in_a * spp_in_b) / (0.5 * (sum(spp_in_a) + 
                                                            sum(spp_in_b)))
                    div_result.append((ind_a, ind_b, dist, QS))

            result.append(np.array(div_result))

        return result

    def get_sp_centers(self, div_list):
        '''
        Get coordinate of center of plots in landscape gridded according to 
        divisions in div_list
        '''
        # TODO: Did not confirm that this works for x_min and y_min > 0
        sp_centers = []
        for div in div_list:
            sp_width = self.p_width / float(div[0])
            sp_height = self.p_height / float(div[1])

            div_sp_cent = []
            for sp_x in xrange(0, div[0]):  # Same sp order as SAD_grid
                x_origin = (self.x_min + sp_x * sp_width)
                x_cent = x_origin + 0.5 * (sp_width - self.unit)
                for sp_y in xrange(0, div[1]):
                    y_origin = (self.y_min + sp_y * sp_height)
                    y_cent = y_origin + 0.5 * (sp_height - self.unit)

                    div_sp_cent.append((x_cent, y_cent))

            sp_centers.append(np.array(div_sp_cent))

        return sp_centers


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


# CLASS NETWORK - Network of plots - NOT YET IMPLEMENTED
class Network():
    pass

