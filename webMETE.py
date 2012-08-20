#!/usr/bin/python

'''
webMETE runs a local webserver to allow a browser interface to the Macroeco scripts.

Classes
-------

'''
from wsgiref.simple_server import make_server, demo_app
#import xml.etree.ElementTree as etree
import sys, os, logging
#from matplotlib.mlab import csv2rec
#from mpl_toolkits.basemap import Basemap
#import matplotlib.pyplot as plt
#import metadata as metadata
import subprocess
import glob
import os
from datetime import datetime


__author__ = "Chloe Lewis"
__copyright__ = "Copyright 2012, Regents of the University of California"
__credits__ = []
__license__ = None
__version__ = "0.5"
__maintainer__ = "Chloe Lewis"
__email__ = "chlewis@berkeley.edu"
__status__ = "Development"

localport = 8000
from wsgiref.util import setup_testing_defaults
from wsgiref.simple_server import make_server

# A relatively simple WSGI application. It's going to print out the
# environment dictionary after being updated by setup_testing_defaults
def simple_app(environ, start_response):
    setup_testing_defaults(environ)

    status = '200 OK'
    headers = [('Content-type', 'text/plain')]

    start_response(status, headers)

    ret = ["%s: %s\n" % (key, value)
           for key, value in environ.iteritems()]
    return ret

def call_script(data='data/formatted/LBRI/LBRI_1998.csv',analysis='code/analyze_my_sad.py',output='projects/BOUNDS'):
    '''Runs the chosen data,analysis,output triple.
     There can be more than one dataset and more than one analysis.
     All analyses will be run with all datasets.

     Note that the GUI stays open: user can choose another triple to run.'''
    METEbase = os.path.dirname(os.path.abspath(__file__)) #Using METE file structure
    later ='''    for afile in self.alist.curselection():
        script = self.alist.realcontent[int(afile)]
        spath = os.path.join(METEbase,script)

        output = self.structure['output']

        dfiles = []
        for dfile in self.dlist.curselection():
            data = self.dlist.realcontent[int(dfile)]
            dpath = os.path.join(METEbase, data)
            dfiles.append(dpath)
'''
    
    dt = datetime.utcnow()
    dfiles = [os.path.join(METEbase,data)]
    spath = os.path.join(METEbase,analysis)
    with open("logfile.txt","a") as log:
        log.write( dt.strftime("%Y %I:%M%p UTC")+" :\t"
                   + spath + "\t" + str(dfiles) +'\n')
    subprocess.Popen(["python", spath] + dfiles, cwd=output,
                     shell=False,stdin=None,stdout=None,close_fds=True)


def AutoMETE_app(environ, start_response):
    '''
    Serves an interface to AutoMETE data and analysis on port %d of the localhost.
    '''%localport
    status = '200 OK'
    headers = [('Content-type','text/html')]
    start_response(status, headers)
    call_script()
    return ["Hello World"]

httpd = make_server('', localport, AutoMETE_app)
print "Running analysis engine. Access through: \n\n http://localhost:%d"%localport
httpd.serve_forever()
httpd.handle_request()
