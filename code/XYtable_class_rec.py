'''Python program to examine SAR and EAR 

Creates an XYtable class which generates SAR and EAR
curves

Date of last update:1/31/2012

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
        a numpy rec array.  This method EXPECTS that the file will be formatted
        similar to cocoli.txt and sherman.txt (**see files for formatting**).
        One rec array is created with three columns (names=(species, x ,y)).  
        
         
        '''
        try:
            #Careful with 'S6' datatype. Refers to a six character string
            dt = dict(names = ('species','x','y'),
                                formats = ('S6',np.float,np.float,))
            self.full_xy_table = np.loadtxt(filename, dt, skiprows = 1, 
                                            usecols = (1,2,3)).view(np.recarray)
            #Making cocoli just a rectangle 
            in_plot = test_if_in_cell(0, 0, 200, 100, self.full_xy_table, 1)
            self.xy_table = self.full_xy_table[in_plot]
            
        except IOError as detail:
                print 'That file does not exist...Details: ' , detail
        
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
    
    
        
                
                
                                                             
                    
                
                
            
                    
            

            
            
            
            
                    
            
        
            
        
           


    