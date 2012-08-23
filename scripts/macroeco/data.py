#!/usr/bin/python

'''
Routines for loading census data and metadata.

Classes
-------
- `DataTable` -- data and metadata for a single censused area
- `Metadata` -- load and parse EML metadata for data file
'''

from __future__ import division
import os, sys
import logging
import numpy as np
import xml.etree.ElementTree as etree
from matplotlib.mlab import csv2rec


class DataTable():
    '''
    Class to hold data table and metadata.

    Parameters
    ----------
    datapath : str
        Path to data. Location of metadata determined from this path.

    Attributes
    ----------
    table : recarray
        Census data loaded from csv
    meta_dict : dict
        Dictionary of metadata associated with table.
    '''

    def __init__(self, datapath):
        '''Initialize XYTable object. See class docstring.'''
        self.table, self.meta = self.dataload(datapath)

    def dataload(self, datapath):
        '''
        Load table and metadata from files.
        
        Parameters
        ----------
        datapath : str
            Path to data table file.
            
        Returns
        -------
        table : recarray
            Census data loaded from csv.
        meta : dict
            Dictionary of metadata associated with table.
        '''
        
        # Check that file is .csv.
        assert datapath[-3:] == 'csv', 'File must be .csv'

        # Load main table
        table = csv2rec(datapath)

        # Load metadata from file, printing warning and returning dict of zeros 
        # if file does not exist.
        metadatapath = datapath[:-3] + 'xml'

        asklist = []
        for name in table.dtype.names:
            asklist.append((name, 'minimum'))
            asklist.append((name, 'maximum'))
            asklist.append((name, 'precision'))

        meta_dict = self.get_metadata(asklist, metadatapath)

        return table, meta_dict

    def get_metadata(self, asklist, metadatapath):
        '''
        Gets the relavent metadata from metadata file.

        Parameters
        ----------
        asklist : list
            List of tuples of column name and attribute, e.g., [('x', 
            'precision'), ('y', 'maximum')].

        metadatapath : str
            Path to the metadata file.
        '''

        meta_raw = Metadata(metadatapath)

        # Check if metadata is missing, if so return None
        if meta_raw.valid_file == False:
            return None

        meta_raw.get_dataTable_values(asklist)
        meta_raw_dict = meta_raw.TableDescr

        meta_dict = {}
        for key, value in meta_raw_dict.iteritems():
            if value is None:
                new_value = value
            else:
                new_value = float(value)

            meta_dict[key] = new_value

        return meta_dict

    def get_sub_table(self, x_st, x_en, y_st, y_en):
        '''
        Returns subtable including only species within rectangle defined by 
        x_st, x_en, y_st, and y_en.
        
        Rectangle is inclusive only on the lower and left sides.
        '''

        x_col = self.head.index('x')
        y_col = self.head.index('y')

        in_sub = np.all([self.table[:,x_col] >= x_st, self.table[:,x_col] < 
                         x_en, self.table[:,y_col] >= y_st, self.table[:,y_col] 
                         < y_en], axis = 0)
        sub_table = self.table[in_sub]
        return sub_table


class Metadata:
    '''
    Metadata values for any analysis. Metadata is intrinsic to a dataset, 
    and is stored using Ecological Markup Language. %citation TODO.
    
    '''

    def __init__(self, datapath):
        ''' Parses the metadata from the EML file sibling to datapath (same 
        location, same filename, .XML extension).

        datapath is expected to be relative (from cwd).
        
        '''#TODO: handle absolute.
        
        dPath, dExtension = os.path.splitext(datapath)
        xmlpath = os.path.abspath(os.path.join(dPath + '.xml'))
    
        self.valid_file = True

        try:
            open(xmlpath)
        except:
            logging.info('Missing or invalid metadata file at %s' % xmlpath)
            self.valid_file = False
        
        try:
            self.doc = etree.ElementTree(file=xmlpath)
            self.root = self.doc.getroot()
        except:
            logging.info('Error parsing metadata file at %s' % xmlpath)
            self.valid_file = False

	    self.TableDescr = None


    def get_dataTable_values(self, asklist):
        ''' The asklist is pairs of column names and sub-attributes of that 
        column, e.g. ('gx','precision').

        Metadata.TableDescr = {(columnName, subAttribute):subAttributeValue}'''
        self.TableDescr = {}
        for request in asklist:
            print 'looking for %s in %s'%(request[1],request[0])
            a = self.get_columnAtt_by_name(request[0])
            if a is None:
                v = None
            else:
                v = self.get_subattribute_value(a, request[1])
            self.TableDescr[request] = v


    def get_columnAtt_by_name(self, attribName):
        '''
        Returns list of XML elements of type attribute with the requested 
        attributeName.
        '''

        attributes = self.root.findall('.//dataTable/attributeList/attribute')
        for a in attributes:
            if a.find('.//attributeName').text == attribName:
                return a
            #print 'No match found for ', attribName 

    def get_subattribute_value(self, att, subatt):
        '''
        Returns the value of the subatt of the given att
        '''
        try:
            value = att.find('.//%s'%subatt).text
            return value
        except AttributeError:
            return None

    def get_coverage_region(self):
        '''
        Returns a tuple of the limits of the physical area of the dataset:  
        NESW.
        '''
        coords = self.root.find('.//coverage/geographicCoverage/' + 
                                'boundingCoordinates')
        bounds = []
        for d in ('north','east','south','west'):
            bounds.append(float(coords.find('%sBoundingCoordinate'%d).text))
        return bounds

    def get_title(self):
        '''Extracts the title of the dataset'''
        return self.root.find('.//dataset/title').text
