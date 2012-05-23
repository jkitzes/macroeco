'''
Initial script for examining commonality of multiple data sets
'''

import numpy as np
import matplotlib.pyplot as plt
from data import *
from empirical import *

# Generate data for BCI. Note that there are 22 spp w/ 0 abund.
div_list_bci = [(10,10), (20,20)]
bci_xy = load_data('../data/archival/BCIS_temp/bci5.csv')
bci_sparse = xytable_to_sparse(bci_xy, compress = False)
bci_patch = Patch(bci_sparse, (0,999.9), (0,499.9), 0.1)
bci_QS = bci_patch.QS_grid(div_list_bci)

# Generate data for serpentine
div_list_lbr = [(4,4), (16,16)]
lbr_xy = load_data('../data/formatted/LBRI_1998/LBRI1998_all_xy.csv')
lbr_sparse = xytable_to_sparse(lbr_xy, compress = False)
lbr_patch = Patch(lbr_sparse, (0,15), (0,15), 1)
lbr_QS = lbr_patch.QS_grid(div_list_lbr)

# Generate data for UCSC
div_list_ucsc = [(4,4), (20,20)]
ucsc_xy = load_data('../data/archival/UCSC2007_temp/FERP07data_xy.csv')
ucsc_sparse = xytable_to_sparse(ucsc_xy, compress = False)
ucsc_patch = Patch(ucsc_sparse, (0,199.9), (0,299.9), .1)
ucsc_QS = ucsc_patch.QS_grid(div_list_ucsc)

# Define function to calculate mean commonality at each unique distance
def mean_comm(data):
    all_comms = data[:,3]
    unique_dists = np.unique(data[:,2])
    unique_mean_comms = []
    for dist in unique_dists:
        unique_mean_comms.append(np.mean(all_comms[dist == data[:,2]]))
    plt.scatter(unique_dists, unique_mean_comms, c = 'red')

# Plot BCI
plt.figure()
for i, div in enumerate(div_list_bci):
    plt.subplot(2,1,i+1)
    data = bci_QS[i]
    plt.scatter(data[:,2], data[:,3], c = 'gray')
    plt.xlim(0,1200)
    mean_comm(data)

# Plot LBR
plt.figure()
for i, div in enumerate(div_list_lbr):
    plt.subplot(2,1,i+1)
    data = lbr_QS[i]
    plt.scatter(data[:,2], data[:,3], c = 'gray')
    mean_comm(data)

# Plot UCSC
plt.figure()
for i, div in enumerate(div_list_ucsc):
    plt.subplot(2,1,i+1)
    data = ucsc_QS[i]
    plt.scatter(data[:,2], data[:,3], c = 'gray')
    mean_comm(data)
