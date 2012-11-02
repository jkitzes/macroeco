#!/usr/bin/python

'''
Routines for loading census data and metadata.

Classes
-------
- `DataTable` -- data and metadata for a single censused area
- `Metadata` -- load and parse EML metadata for data file
'''

from __future__ import division
import os
import logging
import numpy as np
import xml.etree.ElementTree as etree
from matplotlib.mlab import csv2rec


class DataTable:
    '''
    Class to hold data table and metadata.

    Parameters
    ----------
    data_path : str
        Path to data - location of metadata determined from this path.

    Attributes
    ----------
    asklist : list
        A list of tuples of column name and attribute, e.g., [('x', 
        'precision'), ('y', 'maximum')], that defines the columns and 
        parameters that are needed for analysis. Defined in data_load method.
    table : recarray
        Census data table.
    meta : dict
        Dictionary of metadata needed for analysis. Needed variables for each 
        column are defined in asklist 
    '''

    def __init__(self, data_path):
        '''Initialize DataTable object. See class docstring.'''

        self.table, self.meta = self.data_load(data_path)


    def data_load(self, data_path):
        '''
        Load data and metadata from files.
        
        Parameters
        ----------
        data_path : str
            Path to data table file.
            
        Returns
        -------
        table : recarray
            Census data table.
        meta : dict
            Dictionary of metadata associated with table.
        '''
        
        # Check that file is csv
        assert data_path[-4:] == '.csv', 'File must be .csv'

        # Load main table - dtype detected automatically
        table = csv2rec(data_path)

        # Store asklist defining columns and fields needed for analysis.
        # asklist is
        self.asklist = []
        for name in table.dtype.names:
            self.asklist.append((name, 'minimum'))
            self.asklist.append((name, 'maximum'))
            self.asklist.append((name, 'precision'))
            self.asklist.append((name, 'type'))  
        
        # Load metadata from file
        meta = Metadata(data_path, self.asklist).meta_dict

        return table, meta


    def get_subtable(self, subset):
        '''
        Return subtable matching all conditions in subset.

        Parameters
        ----------
        subset : dict
            Dictionary of conditions for subsetting data (see description in 
            Patch Class docstring).

        Returns
        -------
        subtable : ndarray
            Subtable with records from table meeting requirements in subset.

        '''

        # If no subset, return original table
        if subset == {}:
            return self.table
        
        # Declare array to track valid rows of table
        valid = np.ones(len(self.table), dtype=bool)

        # TODO: Add ability to do logical or - and is just multiple subsets on 
        # same column.
        for key, value in subset.iteritems():
            if type(value) is not type(['a']):  # Make all iterables
                value = [value]
            
            # Merge tuples into a string
            merged_values = []
            for val in value:
                try: # check if val[1] is a string
                    eval(str(val[1]))
                    merged_values.append(val[0] + str(val[1]))
                except:
                    merged_values.append(val[0]  + "'" +  val[1] + "'")
                    
            for this_value in merged_values:
                if this_value != "=='whole'":
                    this_valid = eval("self.table[key]" + this_value)
                    valid = np.logical_and(valid, this_valid)

        subtable = self.table[valid]
        return subtable


class Metadata:
    '''
    Metadata values for any analysis stored using Ecological Metadata Language.

    Parameters
    ----------
    data_path : str
        Path to csv data file. Metadata file must be in same dir, with same 
        filename, but with .xml extension.
    
    Attributes
    ----------
    valid_file : bool
        Whether valid metadata file was found.
    root : object
        Root of Element Tree representation of metadata xml file.
    meta_dict : dict
        Dictionary of metadata with values given by asklist.

    '''

    def __init__(self, data_path, asklist):
        '''Initialize Metadata object. See class docstring.'''
       
        # Get path to metadata file
        data_path, data_extension = os.path.splitext(data_path)
        xml_path = os.path.abspath(os.path.join(data_path + '.xml'))
   
        # Determine if metadata file is valid and if so store self.root
        self.valid_file = True

        try:
            open(xml_path)
        except:
            logging.info('Missing or invalid metadata file at %s' % xml_path)
            self.valid_file = False
        
        try:
            self.root = etree.ElementTree(file=xml_path).getroot()
        except:
            logging.info('Error parsing metadata file at %s' % xml_path)
            self.root = None
            self.valid_file = False
        
        # Check if metadata file is missing or invalid, if so return None
        if self.valid_file == False:
            self.meta_dict = None
        else:
            self.meta_dict = self.get_meta_dict(asklist)


    def get_meta_dict(self, asklist):
        '''
        Parse metadata dictionary from xml file.
        
        Parameters
        ----------
        asklist : list
            A list of tuples of column name and attribute, e.g., [('x', 
            'precision'), ('y', 'maximum')], that defines the columns and 
            parameters that are needed for analysis.

        Returns
        -------
        meta_dict : dict
            Dictionary of metadata values for each item in asklist, in form 
            {('column_name', 'element'): value}. column_name in data table is 
            equivalent to attribute in xml.
        '''

        # TODO: Column attribute will be None if either column entry does not
        # exist in metadata or if column entry exists but attribute is missing. 
        # We may want to distinguish these, perhaps just with logging.

        # Populate dictionary of metadata values for asklist items
        meta_dict = {}

        for item in asklist:
            # Get list of all elements for this attribute
            all_elements = self.get_all_elements(item[0])
            #import pdb; pdb.set_trace()

            # Get value of element for this attribute if it exists
            if all_elements is None:
                value = None
            else:
                value = self.get_element_value(all_elements, item[1], item[0])
            
            # Eval value if possible and log outcome
            try:
                value = eval(value)
                value_type = str(type(value)).split("'")[1]
                logging.debug('Metadata value %s, %s evaluated to %s' % 
                              (item[0], item[1], value_type))
            except:
                logging.debug('Metadata value %s, %s left as string' % 
                              (item[0], item[1]))

            # Store value for this item
            meta_dict[item] = value

        return meta_dict


    def get_all_elements(self, attribute):
        '''Returns list of XML elements of type attribute for attribute.'''

        attributes = self.root.findall('.//dataTable/attributeList/attribute')
        for a in attributes:
            if a.find('.//attributeName').text == attribute:
                return a


    def get_element_value(self, all_elements, element_name, col_name):
        '''Returns value of attribute_name from all_attributes list.'''
        if element_name == 'type':
            if len(all_elements.findall('.//dateTime')) == 1:
                return 'ordinal'
            elif len(all_elements.findall('.//interval')) == 1:
                return 'interval'
            elif len(all_elements.findall('.//ordinal')) == 1:
                return 'ordinal'
            elif len(all_elements.findall('.//nominal')) == 1:
                return 'nominal'
            elif len(all_elements.findall('.//ratio')) == 1:
                return 'ratio'
            else:
                logging.warning("Could not find recognizable column type. " +\
                             "Setting type of column name '%s' to ordinal." %\
                             col_name)
                return 'ordinal'
        else:
            try:
                value = all_elements.find('.//%s' % element_name).text
                return value
            except AttributeError:
                return None


    def get_physical_coverage(self):
        '''Returns a tuple of physical limits of the dataset (NESW).'''
        coords = self.root.find('.//coverage/geographicCoverage/' + 
                                'boundingCoordinates')
        bounds = []
        for d in ('north','east','south','west'):
            bounds.append(float(coords.find('%sBoundingCoordinate'%d).text))
        return bounds


    def get_title(self):
        '''Extracts the title of the dataset. Not currently used.'''
        return self.root.find('.//dataset/title').text
