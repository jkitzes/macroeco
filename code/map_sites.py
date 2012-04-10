'''Scans all ../data/formatted/*.xml and maps their locations.'''

import os, glob
import macroeco.metadata as metadata
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt

format_root = os.path.join(os.path.split(os.path.split(os.path.abspath(__file__))[0])[0],
                           'data/formatted/*/*.xml')
files = glob.glob(format_root)

lats = []
lons = []
names = []
for f in files:
    meta = metadata.Metadata(f)
    bounds =  meta.get_coverage_region()

    lats.append(bounds[0])
    lons.append(bounds[1])
    #    names.append(meta.get_title()) #Too long.
    fname, fext = os.path.splitext(os.path.split(f)[-1])
    names.append(fname)


m = Basemap(projection='cyl', lat_0=50, lon_0=-100,
            llcrnrlon=min(lons)-5, llcrnrlat=min(lats)-5,
            urcrnrlon=max(lons)+5, urcrnrlat=max(lats)+5,
            resolution='l')

m.bluemarble()
m.drawcoastlines()
m.drawcountries()
#m.fillcontinents(color='beige') #Not with bluemarble
m.drawmapboundary()

x, y = m(lons, lats)
m.plot(x, y, 'wo')
for n, xpt, ypt in zip(names,x,y):
    if n == 'BCIS': ypt += 1 #Cleanup for crowded areas 
    if n == 'SHER': ypt += 2
    plt.text(xpt+.5,ypt+.5,n,color='white')
plt.title('Field sites of testable data')
plt.savefig('sites_map.png')
plt.show()
    
