import sys
from matplotlib.mlab import csv2rec
import matplotlib.pyplot as plt


data = csv2rec(sys.argv[1],names=None)
outputID = sys.argv[2]
   
x = data.dtype.names[1]
y = data.dtype.names[2]
phenom = data.dtype.names[0]

fig= plt.figure()
ax = fig.add_subplot(111)
ax.scatter(data[x],data[y],s=data[phenom]*10)
ax.set_xlabel(x)
ax.set_ylabel(y)
ax.set_title(phenom)
fig.savefig(outputID+'.png')
plt.show()


