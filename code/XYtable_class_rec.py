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


def dense_to_sparse(dense):
    ''' Convert dense plot data into sparse plot data '''
    sp_y = dense.shape[0]
    sp_x = dense.shape[1]
    nspp = dense.shape[2]

    sparse = []
    for spp in xrange(0, nspp):
        for x in xrange(0, sp_x):
            for y in xrange(0, sp_y):
                if dense[y, x, spp] > 0:
                    sparse.append([spp, x, y, dense[y, x, spp]])

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
test_sparse = dense_to_sparse(test_dense)


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



    #assuming for now that I am only passing one type of shape
    #Will need to adjust parameters later  
    #Going to implement this for only handling rectanglar plots
    def sar_method(self, shape_name='rectangle', metric=(2,1), 
                    area_range=range(0,10000, 100),table='Yes'):
        '''This function calculates the SAR curve for a given data set
        and returns an array with species counts for a given area.
        
        If one sets the keyword table to 'Yes', this function will 
        return a list such that each entry in the list is a 2D
        array in which the number of each individual species
        occuring within a sampling grid is displayed for each samples
        grid of that specific area.
        
        shape_name is a string that takes in the shape of the 
        sampling grid.
        
        metric is a tuple that takes in the ratio of a rectangle 
        sampling grid
        
        area_range is a list that contains the areas to be sampled
        
        '''
        #Constants
        SAMPLE_SIZE = 10
        X_MAX = 200
        Y_MAX = 100
        SCALE = 1 
        DIFF = 0.1 # space between two samples in linear space
        
        if shape_name == 'rectangle':
            #Define the variables
            area_array = (np.array(area_range)) * (SCALE**2)
            x_max_convert = X_MAX * SCALE 
            y_max_convert = Y_MAX * SCALE
            buffer = (SCALE * DIFF)/2.0
            width_ratio = metric[0]
            length_ratio = metric[1]
            species_counts = []
            sad_table_list = []
            #Make an array of all possible widths and lengths
            wid_len_array = calculate_width_length(width_ratio, length_ratio, area_array)
            
            #Make two arrays of all possible corner points
            start = 0 - buffer
            stop_x = x_max_convert + buffer
            stop_y = y_max_convert + buffer
            x_corner_array = np.linspace(start, stop_x, ((stop_x - start)/(DIFF * SCALE)) + 1)
            y_corner_array = np.linspace(start, stop_y, ((stop_y - start)/(DIFF * SCALE)) + 1)
            
            #This for loop iterates through each different sample area
            for i in xrange(len(area_range)):
                x_boundary = round((x_max_convert + buffer) - wid_len_array[i,0], 1)
                y_boundary = round((y_max_convert + buffer) - wid_len_array[i,1], 1)
                #Possible values to sample from
                x_sample_values = x_corner_array[(x_corner_array <= x_boundary)]
                y_sample_values = y_corner_array[(y_corner_array <= y_boundary)]
                grid_samples = np.zeros(SAMPLE_SIZE) 
                sad = set_up_rec_array(self.xy_table, SAMPLE_SIZE)
               
                #This for loop takes SAMPLE_SIZE samples and determines
                #how many species are in the grid    
                for s in xrange(SAMPLE_SIZE):
                    rand_x = random.sample(x_sample_values, 1)
                    rand_y = random.sample(y_sample_values, 1)
                    in_cell = test_if_in_cell(rand_x, rand_y, wid_len_array[i,0],
                                              wid_len_array[i,1], self.xy_table,
                                              SCALE)
                    species_in_cell = self.xy_table.species[in_cell]
                    
                    if table == 'No':
                        unique_species = set(species_in_cell)
                        grid_samples[s] = len(unique_species)
                    if table == 'Yes':
                        for m in xrange(len(sad.species)):
                            
                            if area_array[i] == 0:
                                sad['sample_' + str(s + 1)][m] = 0
                            else:
                                sad['sample_' + str(s + 1)][m] = sum((sad.species[m] 
                                                                    == species_in_cell))
                if table == 'No':        
                    species_counts.append([area_range[i], np.mean(grid_samples)])
                if table == 'Yes':
                    sad_table_list.append(sad)
            
            if table == 'No':
                return np.array(species_counts) 
            if table == 'Yes':
                return sad_table_list
                
                
    def ear_method(self, shape='rectangle', ratio=(2,1), 
                    areas=range(0,20001,1000)):
        '''This function calculates the endemic species in 
        each sampling area.  It uses the SAR identity to 
        calculate. 
        
        Returns an 2D sample array with area and ear counts
        '''
        TOTAL_SPECIES = len(set(self.xy_table['species']))
        X_MAX = 200
        Y_MAX = 100
        TOTAL_AREA = X_MAX * Y_MAX
        
        if shape == 'rectangle':
            area_convert = list(TOTAL_AREA - np.array(areas)) 
            sar_counts = self.sar_method(shape_name=shape, metric=ratio, 
                                        area_range=area_convert, table='No')
            ear_counts = np.zeros((len(areas),2))
            ear_counts[:,0] = np.array(areas)
            ear_counts[:,1] = TOTAL_SPECIES - sar_counts[:,1]
            
            return ear_counts
                            
                
                
                
                
def calculate_width_length(wr, lr, area):
    '''This function takes in a width:length
    ratio and a array of areas and returns an
    array consisting of the appropriate widths
    and lengths for those areas.
    
    '''
    width_length_array = np.zeros((len(area), 2))
    scal_fact = (lambda x: (x/(wr*lr))**0.5)(area)
    width_length_array[:,0] = wr * scal_fact
    width_length_array[:,1] = lr * scal_fact
    return width_length_array
        
def test_if_in_cell(rx, ry, wid, len, point_array, scale):
    '''This function takes in two points, a width
    a length, and a list of points and returns a 
    boolean vector which says which of the points
    in the list are in the rectangle formed by 
    the two points and the length and width
    
    '''
    is_in_cell = ((point_array['x'] * scale) >= rx)
    is_in_cell *= ((point_array['x'] * scale) <= (rx + wid))
    is_in_cell *= ((point_array['y'] * scale) >= ry)
    is_in_cell *= ((point_array['y'] * scale) <= (ry + len))
    return is_in_cell
    
def set_up_rec_array(rec_array, sample):
    '''This function sets up the sad rec array.
    It takes in a rec_array, performs a few 
    manipulations on the rec array and returns another
    rec array that is an sad grid with samples + 1
    columns.
    
    sample is an int and rec_array is a rec array
    
    '''
    
    species_dtype = np.dtype([('species', 'S6')])
    species_array = np.array(np.unique(rec_array.species))
    full = np.array(species_array, species_dtype).view(np.recarray)
    num_array = np.zeros(len(full.species))
    for i in xrange(sample):
        full = plt.mlab.rec_append_fields(full, 'sample_' + str(i + 1), 
                                           num_array)
    return full
    
    
        
                
                
                                                             
                    
                
                
            
                    
            

            
            
            
            
                    
            
        
            
        
           


    
