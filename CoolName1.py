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
import webbrowser
from numpy.random import random

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

localport = 8000 # WARNING:  hardcoded into CoolName1.html
layout = Template(open('CoolName1.txt').read())
css_file = open('CoolName1.css').read()
output_path = '.'
processes = {}

def CoolName1_app(environ, start_response):
    '''
    Runs CoolName1 on a local web server; access through
    http://localhost:%d


    Dispatch the current request to
    the functions from above and store the regular expression
    captures in the WSGI environment as  `myapp.url_args` so that
    the functions from above can access the url placeholders.

    If nothing matches call NIY().
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

        #save out separate d & f page as scripts/modname.txt
        del module
        
    #save pretty list as script_list.txt
        
def css(environ, start_response):
    '''Read in the CSS file that styles our GUI.'''
    start_response('200 OK', [('Content-Type', 'text/css')])
    return css_file

def index(environ, start_response):
    '''Home page: explains what to do with CoolName1.'''
    start_response('200 OK', [('Content-Type', 'text/html')])
    explanation = '''<p>CoolName1 can analyze ecological data in several ways,
                   using scripts we have written and tested. Input parameters and
                   output files are stored together in a Project directory so you
                   always have a record of how each result was generated. </p>
                   <p>The Setup page lists scripts, project directories, and
                   formatted data for you to choose. Each time you start a script,
                   another tab will open to remind you where the results are.</p>
                   <p>This application runs on your home computer, and your data
                   and results remain private. </p>
                   <p>You can also run any of these scripts from the command-line,
                   or write your own scripts using our mathematical and helper
                   functions.</p>
                   <p>Contact the CoolName1 team at <a href=
                   "mailto:CoolName1@berkeley.edu">CoolName1</a>.</p>'''
    return layout.safe_substitute(maincontent = explanation,
                                  localport = localport)

def scripts(environ, start_response):
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

def new_procID():
    '''This can't be the simplest way. Need a unique number for each subprocess.'''
    later = '''pID = 1; old = 1
    yield old
    while 1:
        pID = old + 1
        yield pID
        old = pID'''
    return 1

def setup_all(environ, start_response):
    '''Present scripts and recent project directories to user;
    user chooses script, project, and data.'''
    start_response('200 OK', [('Content-Type', 'text/html')])
    scripts = open('script_list.txt').read()
    RecentList = open('CoolName1.recent').readlines()
    projects = ' '.join(map(dir_to_radio,RecentList))
    proc_number = new_procID()    
    fill = '''<form name="setup_all" action="results" target="_blank{!s}" method="get">
           <h1>Available scripts:</h1>
           <p>{!s}</p>
           <h1>Recent project directories:</h1>
           <p>{!s}</p>
           <p>Demo projects already specify data.</p>
           <input type="submit" value="Run" />
           <input type="hidden" name="procID" value={!s} />
           </form>'''.format(proc_number, scripts, projects, proc_number)
    return layout.safe_substitute(maincontent = fill, localport = localport)

def run(environ, start_response):
    '''Parse run elements out of parameters.xml;
    offer a list of runs;
    execute chosen runs as subprocesses.'''
    later = '''           <h1>Data:</h1>
           <p>Enter the path relative to the project directory.</p>
           <p><input type=text name="data">
'''
    fields = cgi.parse_qs(environ['QUERY_STRING'])
    print fields
    try:
        pf = open('parameters.xml', 'r')
        print pf.read()
    except:
        start_response('404 NOT FOUND', [('Content-Type', 'text/html')])
        return layout.safe_substitute(maincontent = '<h1>No parameters.xml file</h1><p>We need a parameters.xml file in the project directory to set up and run the analyses. Examples are in the /projects/demos subdirectories of CoolName1. Project directory:<code>%s</code></p>'%output_path, localport=localport)
 

def results(environ, start_response):
    '''Run a script. Project directory and script passed in HTML request;
    if script has required parameters, including data,
    they must be defined in parameters.xml in the project directory.'''
    #print environ
    fields = cgi.parse_qs(environ['QUERY_STRING'])
    print fields
    #scriptname = environ['HTTP_REFERER'].split('/')[-1]+".py" #TODO: ok on Windows?
    try:
        scriptname = fields['script'][0]
    except:
        logging.error('No script specified')
    try:
        output_path = fields['project'][0]
    except:
        logging.error('No project directory specified')
    try:
        procID = fields['procID'][0]
    except:
        logging.warning('No process ID')
    start_response('200 OK', [('Content-Type', 'text/html')])

    spath = os.path.dirname(os.path.abspath(__file__))
    spath = os.path.join(spath, 'scripts', scriptname+'.py')

    #process_output = file('results.tmp','w') 
    p = subprocess.Popen(['python',spath],
                      cwd=output_path, shell=False, stdin=None,
                      stdout=None, close_fds=True)
    print type(p)
    processes[procID] = p    
    # Should this be an iterator? (see StackOverflow); yield Starting... , then status?

    return layout.safe_substitute(maincontent = '''<h2>Running %s</h2>
                    <p>Results stored in: %s</p>
                    <p>Stop process: <form action="kill/%s">
                    <input type="submit" value="Cancel" />
                    </form> </p>
                    '''%(scriptname,output_path,procID), localport = localport)

def kill(environ, start_response):
    '''Terminates a specific subprocess.'''
    args = environ['myapp.url_args']
    if args:
        procID = cgi.escape(args[0])
        print 'procID = ', procID
    else:
        return layout.safe_substitute(maincontent = "Process ID lost",
                                      localport = localport)
    processes[procID].kill()
    #subprocess.Popen.terminate()
    
def dir_to_radio(dirstring):
    o = Template('''<input type=radio name="project" value=$value />
                     $pretty<br />''')
    p = dirstring.rstrip()
    mpath = os.path.dirname(os.path.abspath(__file__))
    v = os.path.abspath(os.path.join(mpath, p))
    return o.safe_substitute(value=v, pretty = p)

def dir_to_link(dirstring):
    o = Template('''<a href="results?output=$value">
                     $pretty</a><br />''')
    p = dirstring.rstrip()
    mpath = os.path.dirname(os.path.abspath(__file__))
    v = os.path.abspath(os.path.join(mpath, p))
    return o.safe_substitute(value=v, pretty = p)

def project(environ, start_response):
    '''Directories containing a parameters.xml file;
    TODO ship the list with the contents of CoolName1.recent.'''
    RecentList = open("CoolName1.recent").readlines()
    start_response('200 OK', [('Content-Type', 'text/html')])
    if len(RecentList) > 5:
        control = ' '.join(map(dir_to_link,RecentList[-5:]))
    else:
        control = ' '.join(map(dir_to_link,RecentList))
    return layout.safe_substitute(maincontent =
                                  "<h1>Choose project:</h1>%s"%control,
                                   localport = localport)

def explanation(environ, start_response):
    '''Display the long explanation of the script.'''
    start_response('200 OK', [('Content-Type', 'text/html')])
    scriptname = cgi.parse_qs(environ['QUERY_STRING'])['script'][0]
    spath = os.path.dirname(os.path.abspath(__file__))
    spath = os.path.join(spath, 'scripts', scriptname+'.txt')
    fill = open(spath).read()
    return layout.safe_substitute(maincontent=fill,
                                  localport = localport)

    
    

    
    
def NIY(environ, start_response):
    start_response('404 NOT FOUND', [('Content-Type', 'text/html')])
    return ['Not Found or Not Implemented Yet']

# WSGI parses URLs, using them like arguments even if
# they look like filesystem paths. Nb: this declaration needs
# to be after function defs. 
urls = [
    (r'^$', index), # Homepage.
    (r'scripts', scripts), #Index of scripts
    (r'css', css),  # CSS style file.
    (r'setup',setup_all), # Choose script, directory, data
    (r'project', project), # Choose a parameters.xml file
                       # (and therefore an output directory)
    (r'run', run), # Parse available runs out of the parameters.xml,
    (r'kill/(.+)$', kill), # kill the subprocess specified in URL
    (r'data/(.+)$', NIY), # Choose a dataset
    (r'explanation/?$', explanation), 
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
# make_pretty_lists()
# print('Trying to open entry page in browser...')
# subprocess.Popen("CoolName1.html", cwd='.', shell=False, stdin=None, stdout=None, close_fds=True)
httpd = make_server('', localport, CoolName1_app)
print '''

The CoolName1 engine is running now.

Open this URL in your browser:

http://localhost:{!s}

'''.format(localport)
webbrowser.open('http://localhost:{!s}'.format(localport))
httpd.serve_forever()
httpd.handle_request()


