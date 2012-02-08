'''
Routines for loading and converting data.
'''

import numpy as np

#
# LOADING
#

def load_data(filename):
    ''' Load plot and metadata from file. NOT YET COMPLETE.
    
    Data should be returned as a 3D numpy array if dense, or a 4 col numpy 
    array if sparse. If sparse, first col must be an int index of species 
    number. Should also read metadata file so that it can be passed to class 
    Plot later.
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


#
# CONVERSIONS
#

def xytable_to_sparse(xytable):
    ''' Convert xytable data into sparse plot data
    
    Note that an xytable with a fourth column of all '1' values can also be 
    used in place of sparse data matrix, although calculations may be slower.
    '''
    # TODO: Add flag not to compress

    # Add dummy first row to sparse so that equality check sparse_row == row
    # works for first real row.
    sparse_row = [np.array((None, None, None))]
    sparse_count = [None]

    # Loop through rows in xytable and add new row to sparse if this spp and
    # coordinate is not already in sparse. If this combination already exists,
    # increment the associated count.
    for row in xytable:
        # TODO: Return row_ind as below, use implicit false
        if ~np.all(sparse_row == row, axis=1).any():
            sparse_row.append(row)
            sparse_count.append(1)
        else:
            # TODO: sparse_row[np.where(sparse_row == row)]?
            row_ind = np.where(np.all(sparse_row == row, axis=1))[0][0]
            sparse_count[row_ind] += 1

    # TODO: Change this to use .T without extra transpose
    return np.vstack((np.array(sparse_row[1:]).transpose(),
                      np.array(sparse_count[1:]))).transpose()


def dense_to_sparse(dense, unit):
    ''' Convert dense plot data into sparse plot data
    
    Dense data must already be numpy array or a similar data type with three 
    dimensions.
    '''

    # Record number of subplots in y and x direction, and num of species
    sp_y = dense.shape[0]
    sp_x = dense.shape[1]
    nspp = dense.shape[2]

    # Loop through each cell of dense and add row to sparse
    sparse = []
    for spp in xrange(0, nspp):
        for x in xrange(0, sp_x):
            for y in xrange(0, sp_y):
                if dense[y, x, spp] > 0:
                    sparse.append([spp, x * unit, y * unit, dense[y, x, spp]])

    return np.array(sparse)


def sparse_to_dense(sparse, x_minmax, y_minmax, unit, nspp):
    ''' Convert sparse plot data into dense plot data
    
    Note that the resulting dense data can be VERY large is the extent of x and 
    y is large and unit is small. Use with caution. If the desired operation is 
    to take a sparse table and calculate the SAD (third dimension of dense) at 
    some coarser scale, use SAD_grid instead.
    '''

    # Calculate dimension of dense in x and y direction
    x_min = x_minmax[0]
    y_min = y_minmax[0]
    nx = (x_minmax[1] - x_minmax[0] + unit) / float(unit)
    ny = (y_minmax[1] - y_minmax[0] + unit) / float(unit)

    # Fill cells in dense
    dense = np.zeros((ny, nx, nspp))
    for row in sparse:
        dense[(row[2] - y_min) / unit, (row[1] - x_min) / unit,
              row[0]] += row[3]

    return dense


