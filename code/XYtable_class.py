'''Python program to examine SAR and EAR 

Creates an XYtable class which generates SAR and EAR
curves

Date of last update:1/26/2012

'''

import numpy as np
import exceptions
import random
import math
import matplotlib as plt




class XYtable:
    '''The XYtable class
    
    This class will perform SAR and EAR operations
    on the species occurence data
    
    '''
    def __init__(self, filename):
        '''This method takes in a filename and turns this file into 
        a numpy array.  This method EXPECTS that the file will be formatted
        similar to cocoli.txt and sherman.txt (**see files for formatting**).
        Two arrays are created: a 2D array with the (x,y) coordinates of the 
        trees called xy_table and a 1D array with the names of the species 
        called species_list.
        
        In xy_table, col[0] contains x values and col[1] contains y values 
        
        
        '''
        try:
            #This array contains the (x,y) values 
            self.full_xy_table = np.loadtxt(filename, dtype = float, skiprows = 1, 
                                       usecols = (2,3))
            
            #This array contains the species list                        
            self.full_species_list = np.loadtxt(filename, dtype = str, skiprows = 1,
                                            usecols = (1,))
            
            #This array contains the species species_list represented
            #as ints
            self.full_spec_list_ints = make_names_to_ints(self.full_species_list)
                                             
                        
            #These lines are just because we are cutting off cocoli
            self.in_plot = test_if_in_cell(0, 0, 200, 100, self.full_xy_table)
            self.xy_table = self.full_xy_table[self.in_plot]
            self.species_list = self.full_species_list[self.in_plot]
            self.spec_list_ints = self.full_spec_list_ints[self.in_plot]
            
            
        except IOError as detail:
                print 'That file does not exist...Details: ' , detail
        
    #assuming for now that I am only passing one type of shape
    #Going to implement this for only handling rectanglar plots
    def sar_method(self, shape_name='rectangle', metric=(2,1), 
                    area_range=range(0,10000, 100),table='No'):
        '''This function calculates the SAR curve for a given data set
        and returns an array with species counts for a given area.
        
        If one sets the keyword table to 'Yes', this function will 
        return a list such that each entry in the array is a 2D
        array in which the number of each individual species
        occuring within a sampling grid is displayed for each sample
        grid of that specific area.
        
        shape_name is a string that takes in the shape of the 
        sampling grid.
        
        metric is a tuple that takes in the ratio of of a rectangle 
        sampling grid
        
        area_range is a list that contains the areas to be sampled
        
        '''
        SAMPLE_SIZE = 10
        X_MAX = 200
        Y_MAX = 100
        SCALE = 1
        DIFF = 0.1 # space between two samples in linear space
        
        if shape_name == 'rectangle':
            area_array = (np.array(area_range)) * (SCALE**2)
            x_max_convert = X_MAX * SCALE 
            y_max_convert = Y_MAX * SCALE
            buffer = (SCALE * DIFF)/2.0
            width_ratio = metric[0]
            length_ratio = metric[1]
            species_counts = []
            sad_tables_list = []
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
                sad_table = set_up_sad_table(self.spec_list_ints, SAMPLE_SIZE) 
                #This for loop takes SAMPLE_SIZE samples and determines
                #how many species are in the grid
                for s in xrange(SAMPLE_SIZE):
                    rand_x = random.sample(x_sample_values, 1)
                    rand_y = random.sample(y_sample_values, 1)
                    in_cell = test_if_in_cell(rand_x, rand_y, wid_len_array[i,0],
                                                wid_len_array[i,1], (self.xy_table
                                                * SCALE))
                    species_in_cell = self.spec_list_ints[in_cell] 
                    if table == 'No':
                        unique_species = set(species_in_cell)
                        grid_samples[s] = len(unique_species)
                    if table == 'Yes':
                        for m in xrange(len(sad_table[:,0])):
                            if area_array[i] == 0:
                                sad_table[m, s + 1] = 0
                            else:
                                sad_table[m, s + 1] = sum((sad_table[m,0] 
                                                            == species_in_cell))
                    
                if table == 'No':
                    species_counts.append([area_range[i], np.mean(grid_samples)])
                if table == 'Yes':
                    sad_tables_list.append(sad_table)
            
            if table == 'No':
                return np.array(species_counts) 
            if table == 'Yes':
                return (area_range, sad_tables_list)
    
    def ear_method(self, shape='rectangle', ratio=(2,1), 
                    areas=range(1000,20001,1000)):
        '''This function calculates the endemic species in 
        each sampling area.  It uses the SAR identity to 
        calculate. 
        
        Returns an 2D sample array with area and ear counts
        '''
        TOTAL_SPECIES = len(np.unique(self.spec_list_ints))
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
        
def test_if_in_cell(rx, ry, wid, len, point_array):
    '''This function takes in two points, a width
    a length, and a list of points and returns a 
    boolean vector which says which of the points
    in the list are in the rectangle formed by 
    the two points and the length and width
    
    '''
    is_in_cell = (point_array[:,0] >= rx)
    is_in_cell *= (point_array[:,0] <= rx + wid)
    is_in_cell *= (point_array[:,1] >= ry)
    is_in_cell *= (point_array[:,1] <= ry + len)
    return is_in_cell
    
def make_names_to_ints(s_list):
    '''This function takes in the 
    complete species array for the xy_table
    and converts the species strings to
    ints.  The function return the an
    equivalent species array with all strings
    converted to ints
    
    '''
    
    s_unique = np.unique(s_list)
    s_ints = np.linspace(0, len(s_unique) - 1, 
                         num=len(s_unique))
    
    spec_to_ints = np.zeros(len(s_list))
            
    for s in xrange(len(s_unique)):
        check = (s_unique[s] == s_list)
        for i in xrange(len(check)):
            if check[i]:
                spec_to_ints[i] = s_ints[s]
    
    return spec_to_ints 

def set_up_sad_table(spec_int, size):
    '''This fucntion sets takes in a species 
    array and a sample size (int) and makes and
    returns a SAD table.
    
    '''
    
    s_unique = np.unique(spec_int)
    sad_tab = np.zeros((len(s_unique), size + 1))
    sad_tab[:,0] = s_unique
    return sad_tab
    
                
                
                                                             
                    
                
                
            
                    
            

            
            
            
            
                    
            
        
            
        
           


    