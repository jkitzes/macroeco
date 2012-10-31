#!/usr/bin/python

'''
Manages the details of a reproducible workflow within macroeco. Main Workflow 
class is called with one argument, required_params, and the surrounding script 
must be called with a single sys.argv with the output directory.

Classes
-------
- `Workflow` -- tracks the analysis, data requested, and parameters; maps sites
- `Parameters` -- finds/asks for and stores run names and parameters
'''

import xml.etree.ElementTree as etree
import sys, os, logging
import matplotlib.pyplot as plt
from macroeco.data import Metadata

__author__ = "Chloe Lewis"
__copyright__ = "Copyright 2012, Regents of the University of California"
__credits__ = []
__license__ = None
__version__ = "0.5"
__maintainer__ = "Chloe Lewis"
__email__ = "chlewis@berkeley.edu"
__status__ = "Development"

paramfile = 'parameters.xml'  # Parameter file found in output dir
logfile   = 'logfile.txt'  # Logfile to save in output dir


class Workflow:
    '''
    Manages the details of a reproducible workflow with macroeco scripts.

    Arguments
    ---------
    required_params : dictionary
        Parameters needed for analysis, in form of 
        'parameter_name':'short_description'. All of these parameters must be 
        present in params file in output directory, or analysis will not run. 
        This argument is empty only when no data or parameters are required for 
        a script to run.
    clog : bool
        Whether to log to console in addition to file, False by default
        
    Attributes
    ----------
    script_name : string
        Name of script originating the workflow
    output_path : string
        Path to output directory
    interactive : bool
        Whether the script can pause for user interaction
    runs : dict
        If parameters are needed, sets of parameter values are named runs
    '''

    def __init__(self, required_params={}, clog=False, svers=None):

        # Store script name from command line call
        script_path, script_extension = os.path.splitext(sys.argv[0])
        self.script_name = os.path.split(script_path)[-1]
        self.script_vers = svers

        # Store output directory path - contains params file, log, results
        # TODO: Make more robust to non-absolute path entries
        output_path = os.getcwd()
        self.output_path = output_path

        # Prepare logger
        logging.basicConfig(filename=logfile,   # Add file logging
                            level=logging.INFO, format=('%(asctime)s | '
                            '%(levelname)s | %(filename)s:%(lineno)d | '
                            '%(message)s'), datefmt='%H:%M:%S')

        if clog:  # Add console logging
            console = logging.StreamHandler()
            console.setLevel(logging.INFO)
            format = logging.Formatter('%(levelname)-8s %(message)s')
            console.setFormatter(format)
            logging.getLogger('').addHandler(console)
       
            def excepthook(*args):  # Catch errors to log
                logging.getLogger().error('Error', exc_info=args)
        else:
            def excepthook(*args):  # Catch errors to log + stderr
                logging.getLogger().error('Error', exc_info=args)
                sys.__excepthook__(*args)  # Show err in console if clog False

        sys.excepthook = excepthook  # Define error handler as above

        logging.captureWarnings(True)  # Catch warnings
        
        logging.debug('Creating workflow object')

        # Get parameters from file, including data paths
        assert type(required_params) == type({}), ('Required params must be a' 
                                                   ' dict.')
        self.parameters = Parameters(self.script_name, self.script_vers,
                                     required_params)
        self.interactive = self.parameters.interactive

        
    def single_datasets(self):
        '''
        Generator that yields data files and descriptive parameters.

        Special parameter 'data_path' is a list of locations of data files to 
        use for analysis - if present, map of sites will be generated for each 
        run.

        Yields
        ------
        data_path : string
            Path to data to analyze, relative to current working directory
        output_ID : string
            Concatenates script, run, and dataset identifiers
        run_params : dict
            Dictionary of parameters for each script_name and run
        '''

        def clean_name(fp):  # Extract file name from path
            return os.path.splitext(os.path.split(fp)[-1])[0]

        # Run script on all runs (parameter sets), and data sets
        for run_name in self.parameters.params.keys():
        # TODO: Check for output_ID conflicts (must be unique)

            # Check if data_paths in params. If not, add one empty data_path
            # for the loop below. If so, make a map.
            if len(self.parameters.data_path) == 0:
                logging.debug(('No data paths given for run %s, no map of '
                              'sites created') % run_name)
                self.parameters.data_path[run_name] = ['']
            else:
                make_map(self.parameters.data_path[run_name], run_name)
                
            # Loop through each dataset and yield values for dataset and run
            for data_path in self.parameters.data_path[run_name]:
                abs_data_path = os.path.abspath(os.path.join(self.output_path, 
                                                             data_path))
                output_ID = '_'.join([self.script_name, 
                                      run_name, clean_name(data_path)])
                logging.info('Beginning %s' % output_ID)
                yield (abs_data_path, output_ID, 
                       self.parameters.params[run_name])

        
class Parameters:
    '''
    Load parameters from parameter file in current working directory
    and make available as self.params.
    Checks that all required_params are present and loaded.

    Arguments
    ---------
    script_name : string
        Name of script originating the workflow
    required_params : dictionary
        Parameters needed for analysis, in form of 
        'parameter_name':'short_description'. All of these parameters must be 
        present in params file in output directory for this script_name and 
        run, or analysis will not run. This argument is empty only when no data 
        or parameters are required for a script to run.

    Attributes
    ----------
    script_name : string
        Name of script originating the workflow
    script_vers : string
        Version of script originating the workflow
    interactive : bool
        Whether the script can pause for user interaction
    params : dict
        Dictionary of dictionaries, with each outer key a run name and each 
        outer value a dictionary of parameter names and values for each run.
    data_path : dict
        Dictonarity where keys are run names and values are lists of data paths 
        associated with each run.
        
    '''
    
    def __init__(self, script_name, script_vers, required_params, 
                 output_path=False):

        # Store initial attributes
        self.script_name = script_name
        self.script_vers = script_vers
        self.interactive = False
        self.params = {}
        self.data_path = {}
        if not output_path:
            output_path = os.getcwd()

        # Check that parameter file exists, if not use default values
        try:
            pf = open(paramfile, 'r')
            pf.close()
        except IOError:
            logging.info(('No parameter file found at %s, proceeding without '
                          'parameters') % output_path)
            self.params[''] = {}
            self.data_path[''] = {}
            self.interactive = False
            return

        # Read parameter file
        logging.info('Reading parameters from %s' % os.path.join(output_path, 
                     paramfile))
        self.read_from_xml()

        # Check that all required parameters present in all runs
        if not self.required_params_present(required_params):
            raise IOError('Required parameters missing')
        logging.info('Parameters: %s' % str(self.params))
        logging.info('Data: %s' % str(self.data_path))

        # Evaluate param values into appropriate types
        self.eval_params()

         
    def read_from_xml(self):
        ''' Read parameters from xml file into self.params dictionary. '''

        # Define class for checking keys
        class AllEntities:
            def __getitem__(self, key):
                return key

        # Declare parser object
        # TODO: Without next line, works in iPython, console, not script ??
        parser = etree.XMLParser()
        parser.parser.UseForeignDTD(True)
        parser.entity = AllEntities()

        # Try to open paramfile from output_path
        # TODO: Integration test
        try:
            pml = etree.parse(paramfile, parser=parser).getroot() 
        except etree.ParseError:
            raise IOError('ParseError trying to read %s' % paramfile)
        except:
            raise
        
        # Create params dictionary
        if len(pml) == 0:  # Error if no analyses in param file
            raise IOError('Parameter file %s contains no valid analyses' % 
                          paramfile)
        for analysis in pml:  # Loop analyses looking for script_name
            if analysis.get('script_name') == self.script_name:

                if 'version' in analysis.attrib:  # Set version
                    vers = analysis.get('version')
                    if self.script_vers:  # If got script_vers, check
                        if float(vers) != float(self.script_vers):
                            logging.warning(('Script version does not match ' 
                                             'version in parameters. '
                                             'Continuing, but may fail.'))

                if 'interactive' in analysis.attrib:  # Set interactive
                    ia = analysis.get('interactive')
                    if ia in ['T', 'True', 't', 'true']:
                        self.interactive = True
                    else:
                        self.interactive = False
                else:
                    self.interactive = False

                if len(analysis) == 0:  # Error if no runs in analysis
                    raise IOError(('Analysis found for this script, but no ' 
                                   'valid runs found'))

                run_counter = 1
                for run in analysis.getchildren():  # Loop runs
                    run_name = run.get('name')
                    if run_name is None:
                        run_name = 'run' + str(run_counter)
                        run_counter += 1
                    self.params[run_name] = {}
                    self.data_path[run_name] = []
                    for elt in run.getchildren():  # Loop params in run
                        if elt.tag == 'param':
                            param = elt.get('name')
                            value = elt.get('value')
                            self.params[run_name][param] = value
                        if elt.tag == 'data':
                            data_type = elt.get('type')
                            data_location = elt.get('location')
                            if data_location == 'system':
                                # User responsible for sys paths, security, etc
                                prepend = ''
                            else:
                                prepend = os.path.join('..','..','data',
                                                       'formatted')
                            if data_type == '' or data_type == None:
                                logging.warning(('No data type specified,'
                                                'assuming .csv'))
                                data_type = 'csv'
                            if data_type == 'csv':
                                directory = elt.find('directory').text
                                data_file = os.path.extsep.join((elt.find('file').text,
                                                                'csv'))
                                data_path = os.path.join(prepend,
                                                         directory, data_file)
                                self.data_path[run_name].append(data_path)
                            else:
                                logging.error('Data type {!s} not yet handled; '
                                              'not using this data.'.format(
                                                  data_type))
                            
    def required_params_present(self, req_params):
        ''' Check if any required parameters missing from any runs. '''

        status = 1
        for run_name in self.params.keys():
            run_params = self.params[run_name]
            if not set(req_params.keys()).issubset(set(run_params.keys())):
		logging.error('In run {!s}, missing parameters {!s}'.format(
			       run_name, set(req_params.keys()).difference(set(run_params.keys()))))
                status = 0
        return status


    def eval_params(self):
        '''
        Attempts to evaluate parameters to appropriate types.
        
        If eval() fails, parameter will stay a string, possibly leading to 
        cryptic errors later if there is a typo in a param value.
        '''

        for run_name in self.params.keys():
            for param_name in self.params[run_name].keys():
                try:
                    value = eval(self.params[run_name][param_name])
                    self.params[run_name][param_name] = value
                    value_type = str(type(value)).split("'")[1]
                    logging.debug('In run %s, parameter %s evaluated to %s' % 
                                  (run_name, param_name, value_type))
                except:
                    logging.debug('In run %s, parameter %s left as string' % 
                                  (run_name, param_name))

            
def make_map(data_paths, run_name, whole_globe=False):
    '''
    Makes a map of all sites in run.

    Parameter
    ---------
    data_paths : list
        Paths to data files (csv's). Data location will be extracted from 
        corresponding xml metadata file.
    run_name : str
        Name of run, used as name of map file.
    whole_globe : bool
        If True, map is entire globe. If False, map is "zoomed in" on data 
        locations.

    Returns
    -------
    map_created : bool
        True if file was created, False if a file already existed and none was 
        created.

    Notes
    -----
    Map will be given the name of a run. If multiple runs have the same name, 
    only the map associated with the first run of that name will be saved.

    The label for each site will be the data file base name 
    (e.g., LBRI_2000.csv and LBRI.csv will be LBRI_2000 and LBRI respectively).
    '''

    # Check if Basemap present - if not, log and return
    try:
        from mpl_toolkits.basemap import Basemap
    except:
        logging.debug('Basemap module is not available, no map of data ' + 
                      'locations can be created')
        return False

    # Set map_name
    map_name = 'map_' + run_name + '.png'

    # TODO: Check if run_name is unique
    # Check if map with this run_name already exists
    if os.path.isfile(map_name):
        logging.debug('Map with this run name already exists. New map ' + 
                      'overwriting old one.')

    # Get lat, long, and name of each data set
    lats = []
    lons = []
    names = []

    for path in data_paths:
        x = path[:-3]+'xml'
        
        try:
            meta = Metadata(x, {})
            bounds = meta.get_physical_coverage()
            lats.append(bounds[0])
            lons.append(bounds[1])
            
            fname, fext = os.path.splitext(os.path.split(path)[-1])
            names.append(fname)  # First 4 letters of data set name
        except:
            logging.info('No location data found in %s, no map point '
                         'added.' % x)

    # If no valid location data, return without making a map
    if len(names) == 0:
        return False

    # Set up map
    logging.debug('Creating map for run %s' % run_name)
    if whole_globe:
        m = Basemap(projection='cyl', resolution='i')
    else:
        # 10 degree buffer around min and max lat/long
        m = Basemap(projection='cyl', lat_0=50, lon_0=-100,
            llcrnrlon=min(lons)-10, llcrnrlat=min(lats)-10,
            urcrnrlon=max(lons)+10, urcrnrlat=max(lats)+10,
            resolution='l')

    # Draw features
    m.bluemarble()
    m.drawcoastlines()
    m.drawcountries()
    m.drawmapboundary()

    # Add sites
    x, y = m(lons, lats)
    m.plot(x, y, 'yo')
    for n, xpt, ypt in zip(names,x,y):
        if n == 'BCIS': ypt += 1 # Manual Cleanup for crowded areas
        if n == 'SHER': ypt += 2
        plt.text(xpt+.5,ypt+.5,n,color='yellow')

    plt.savefig(map_name)
    plt.close()
    return True
