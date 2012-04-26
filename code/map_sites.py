'''Scans all ../data/formatted/*.xml and maps their locations.'''

import os, glob
import macroeco.metadata as metadata
from macroeco.workflow import make_map
#from mpl_toolkits.basemap import Basemap
#import matplotlib.pyplot as plt

format_root = os.path.join(os.path.split(os.path.split(os.path.abspath(__file__))[0])[0],
                           'data/formatted/*/*.xml')
files = glob.glob(format_root)

make_map(files, 'sites_map', True)
