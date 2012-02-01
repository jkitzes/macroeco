import sys
from matplotlib.mlab import csv2rec
import matplotlib.pyplot as plt
from datetime import datetime

dt = datetime.utcnow()
with open("logfile.txt","a") as log:
    log.write( dt.strftime("%Y %I:%M%p UTC")+" :" + sys.argv[1] + " using " + sys.argv[0]+'\n')
    
data = csv2rec(sys.argv[1],names=None)

x = data.dtype.names[1]
y = data.dtype.names[2]
phenom = data.dtype.names[0]

fig= plt.figure()
ax = fig.add_subplot(111)
ax.scatter(data[x],data[y],s=data[phenom]*10)
ax.set_xlabel(x)
ax.set_ylabel(y)
ax.set_title(phenom)
plt.show()


