#!/usr/bin/python

'''
CoolName1 runs a local webserver to allow a browser interface to the Macroeco scripts.

Classes
-------

'''
from wsgiref.simple_server import make_server, demo_app
#import xml.etree.ElementTree as etree
import sys, os, logging, re
#from matplotlib.mlab import csv2rec
#from mpl_toolkits.basemap import Basemap
#import matplotlib.pyplot as plt
#import metadata as metadata
import subprocess
import glob, cgi
from datetime import datetime


__author__ = "Chloe Lewis"
__copyright__ = "Copyright 2012, Regents of the University of California"
__credits__ = []
__license__ = None
__version__ = "0.5"
__maintainer__ = "Chloe Lewis"
__email__ = "chlewis@berkeley.edu"
__status__ = "Development"


from wsgiref.util import setup_testing_defaults
from wsgiref.simple_server import make_server
from string import Template

# Could subclass BaseWSGI instead of running as a script,
# but it doesn't gain anything (and is self-ish). 

localport = 8000
layout = Template(open('CoolName1.html').read())



def CoolName1_app(environ, start_response):
    '''
    Runs CoolName1 on a local web server; access through
    http://localhost:%d


    Dispatch the current request to
    the functions from above and store the regular expression
    captures in the WSGI environment as  `myapp.url_args` so that
    the functions from above can access the url placeholders.

    If nothing matches call the `not_found` function.
    '''%localport
    path = environ.get('PATH_INFO', '').lstrip('/')
    for regex, callback in urls:
        match = re.search(regex, path)
        if match is not None:
            environ['myapp.url_args'] = match.groups()
            return callback(environ, start_response)

    return NIY(environ, start_response)    

def make_pretty_lists():
    '''Makes a list of script GUI names and summaries
    from the .py files in the ./scripts directory.
        Also makes up the forms describing each script, while we're there.'''
    # formatting for each line
    entry = Template('<p><a href="setup/$link">$GUI_name</a> $summary</p>')
    scriptlist = ''
    
    # list of .py modules
    for module in os.listdir("scripts"): #structure-dependent
        if module[:-11] == '__init__.py' or module[-3:] != '.py':
            break
        modname = module[:-3]
        print 'modname:', modname
        current_module = __import__('scripts.%s'%modname) 
        print current_module
        # script descriptors get prettied up and stored
        # in the home-page list 
        each = entry(modname, modname.GUI_name, modname.summary)
        scriptlist.append(each)

        #and in the description & form page for each script

        #save out separate d & f page
        del module
        
    #save pretty list
        
def css(environ, start_response):
    '''Read in the CSS file that styles our GUI.'''
    c = open('CoolName1.css').read()
    start_response('200 OK', [('Content-Type', 'text/css')])
    return c

def index(environ, start_response):
    '''Returns a HTML page indexing the scripts and
    giving their summaries.'''
    start_response('200 OK', [('Content-Type', 'text/html')])
    summaries = open('script_list.txt').read()
    return layout.safe_substitute(maincontent = summaries, localport = localport)

def setup_script(environ, start_response):
    '''Presents the long explanation for a chosen script
    and provides a form for setting up a run:
    output directory
    data
    parameter values, if any
    '''
    form = environ['myapp.url_args'][0]
    fill = open('scripts/'+form+'.txt').read()
    start_response('200 OK', [('Content-Type', 'text/html')])
    return layout.safe_substitute(maincontent = fill, localport = localport)

def results(environ, start_response):
    print environ
    fields = cgi.parse_qs(environ['QUERY_STRING'])
    scriptname = environ['HTTP_REFERER'].split('/')[-1]+".py"
    start_response('200 OK', [('Content-Type', 'text/html')])
    #    with open("logfile.txt","a") as log:
    #    log.write( dt.strftime("%Y %I:%M%p UTC")+" :\t"
    #               + spath + "\t" + str(dfiles) +'\n') #TODO: use logger
    spath = os.path.dirname(os.path.abspath(__file__))
    print 'spath:', spath
    spath = os.path.join(spath, 'scripts', scriptname)

    dpath = ' '.join(fields['data'])
    callstring = ["python", spath, dpath]
    print callstring
    subprocess.Popen(callstring, cwd = fields['output'][0], shell=False, stdin=None, stdout=None, close_fds=True)

    # TUrn this into an iterator! (stackoverflow); yield Starting... , then status?
    # possibly list of files written to output?
    
    return layout.safe_substitute(maincontent = "<h2>Running %s</h2><p>Parameters:  %s</p>"%(scriptname,str(fields)), localport = localport)

def call_script(data='',analysis='scripts/sample_script.py',output='scripts/samples/'):
    '''Runs the chosen data,analysis,output triple.
     There can be more than one dataset and more than one analysis.
     All analyses will be run with all datasets.

     Note that the GUI stays open: user can choose another triple to run.'''
    METEbase = os.path.dirname(os.path.abspath(__file__)) #Using METE file structure
    
    dt = datetime.utcnow()
    dfiles = [os.path.join(METEbase,data)]
    spath = os.path.join(METEbase,analysis)
    opath = os.path.join(METEbase, output)

    

def NIY(environ, start_response):
    start_response('404 NOT FOUND', [('Content-Type', 'text/html')])
    return ['Not Found or Not Implemented Yet']

# WSGI parses URLs, using them like arguments even if
# they look like paths. Nb: needs to be after function defs. 
urls = [
    (r'^$', index), # Index of scripts.
    (r'css', css),  # CSS style file.
    (r'setup/(.+)$',setup_script), # Set up a script
    (r'data/(.+)$', NIY), # Choose a dataset
    (r'documentation/?$', NIY), # Local documentation
    (r'results/?$', results) # Point to or display results (if any).
]


def simple_app(environ, start_response):
    '''A simple Web application. Prints the environment dictionary after being updated by setup_testing_defaults.
    '''

    setup_testing_defaults(environ)

    status = '200 OK'
    headers = [('Content-type', 'text/plain')]

    start_response(status, headers)

    ret = ["%s: %s\n" % (key, value)
           for key, value in environ.iteritems()]
    return ret


# Try this if the app fails and it might be the webservice itself:
#httpd = make_server('', localport, simple_app)
#make_pretty_lists()
httpd = make_server('', localport, CoolName1_app)
print "Running CoolName1 analysis. Use your browser to interact: \n\n http://localhost:%d"%localport
httpd.serve_forever()
httpd.handle_request()


