#!/usr/bin/env python

'''
This module contains a minimal metadata writer class for quickly making
metadata

'''


import xml.etree.ElementTree as ET
import os

sub = ET.SubElement


__author__ = "Mark Wilber"
__copyright__ = "Copyright 2012, Regents of University of California"
__credits__ = "John Harte"
__license__ = None
__version__ = "0.1"
__maintainer__ = "Mark Wilber"
__email__ = "mqw@berkeley.edu"
__status__ = "Development"


class MetaWriter:
    '''
    Writes a metadata file based on the given filename and user input
    
    '''

    def __init__(self, datapath):
        '''
        Class takes in a datafile path name and creates an xml tree using the
        column heading of the recarray generated from the csv file.

        Parameters
        ----------
        datapath : string
            Datafile name

        '''
        assert datapath[-4:] == '.csv', "%s must end in .csv" % (datapath)
        self.filename = datapath.split('.')[0]
        fin = open(datapath, 'r')
        self.column_names = fin.readline().strip().split(',')
        fin.close()
        self.root = ET.Element('eml:eml')
        self.root.attrib = {'packageId' : self.filename, 'system' : 'knb',
            "xmlns:eml" : "eml://ecoinformatics.org/eml-2.1.0", 'xmlns:xsi': 
            "http://www.w3.org/2001/XMLSchema-instance", "xsi:schemaLocation"
            : "eml://ecoinformatics.org/eml-2.1.0 eml.xsd"}
        self.dataset = sub(self.root, 'dataset')
        self.title = sub(self.dataset, 'title')
        self.title.text = "Data set " + os.path.split(datapath)[1]

        self.creator = sub(self.dataset, 'creator')
        self.individualName = sub(self.creator, 'individualName')
        self.surName = sub(self.individualName, 'surName')
        self.surName.text = "None"

        self.contact = sub(self.dataset, 'contact')
        self.individualName2 = sub(self.contact, 'individualName')
        self.surName2 = sub(self.individualName2, 'surName')
        self.surName2.text = "None"

        self.dataTable = sub(self.dataset, 'dataTable')

        self.entityName = sub(self.dataTable, 'entityName')
        self.entityName.text = os.path.split(datapath)[1]

        self.physical = sub(self.dataTable, 'physical')
        self.objectName = sub(self.physical, 'objectName')
        self.objectName.text = os.path.split(datapath)[1]
        self.size = sub(self.physical, 'size')
        self.size.attrib = {'unit' : "byte"}
        self.size.text = str(os.path.getsize(datapath))
            
        # Nested in physical
        self.dataFormat = sub(self.physical, 'dataFormat')
        self.textFormat = sub(self.dataFormat, 'textFormat')
        self.numHeaderLines = sub(self.textFormat, 'numHeaderLines')
        self.numHeaderLines.text = "1"
        self.recordDelimiter = sub(self.textFormat, 'recordDelimiter')
        self.recordDelimiter.text = "#x0A"
        self.attributeOrientation = sub(self.textFormat, 'attributeOrientation')
        self.attributeOrientation.text = "column"
        self.simpleDelimited = sub(self.textFormat, 'simpleDelimited')
        self.fieldDelimiter = sub(self.simpleDelimited, 'fieldDelimiter')
        self.fieldDelimiter.text = ","
        
        self.distribution = sub(self.physical, 'distribution')
        self.online = sub(self.distribution, 'online')
        self.url = sub(self.online, 'url')
        self.url.text = "macroeco://" + os.path.split(datapath)[1]
        

        self.attributeList = sub(self.dataTable, 'attributeList')
        self.attributes = []
        self.attributeTypes = []
        for i, name in enumerate(self.column_names):
            attribute = sub(self.attributeList, 'attribute')
            attributeName = sub(attribute, 'attributeName')
            attributeDefinition = sub(attribute, 'attributeDefinition')
            attributeDefinition.text = "None"
            measurementScale = sub(attribute, 'measurementScale')

            # Default Ordinal
            attributeType = sub(measurementScale, 'ordinal')
            nonNumericDomain = sub(attributeType,'nonNumericDomain')
            textDomain = sub(nonNumericDomain, 'textDomain')
            definition = sub(textDomain, 'definition')
            definition.text = "None"

            attributeName.text = name
            self.attributes.append(attribute)
            self.attributeTypes.append(attributeType)

        self.numberOfRecords = sub(self.dataTable, 'numberOfRecords')
        self.numberOfRecords.text = "Unknown"

    def add_attribute_types(self, typelist):
        '''
        Sets the type of the attribute to either ordinal (categorical) or
        interval (categorical). Initialized in constructor as ordinal. 

        Parameters
        ----------
        typelist : list
            A list of tuples.  Each tuple contains 2 elements: a string and a
            dict. The dict must contain the keyword cat (categorical) or a 
            KeyError will be thrown.

        Example of typelist:

            [('x', {'cat' : True}), ('y' : {'cat' : True}), ('year',
            {'cat' : False}]

        '''

        for item in typelist:
            for attribute in self.attributes:
                tree = ET.ElementTree(attribute)
                att = tree.findall('attributeName')[0]
                if (att.text == item[0]):
                    measure = tree.findall('measurementScale')[0]
                    if item[1]['cat'] == True:
                        if len(measure.findall('interval')) == 1:
                            measure.remove(measure.find('interval'))
                            att_type = sub(measure, 'ordinal')
                            nonNumericDomain = sub(att_type,'nonNumericDomain')
                            textDomain = sub(nonNumericDomain, 'textDomain')
                            definition = sub(textDomain, 'definition')
                            definition.text = "None"

                        elif len(measure.findall('ordinal')) == 1:
                            measure.remove(measure.find('ordinal'))
                            att_type = sub(measure, 'ordinal')
                            nonNumericDomain = sub(att_type,'nonNumericDomain')
                            textDomain = sub(nonNumericDomain, 'textDomain')
                            definition = sub(textDomain, 'definition')
                            definition.text = "None"

                    elif item[1]['cat'] == False:

                        if len(measure.findall('ordinal')) == 1:
                            measure.remove(measure.find('ordinal'))
                            att_type = sub(measure, 'interval')
                            unit = sub(att_type, 'unit')
                            standardUnit = sub(unit, 'standardUnit')
                            standardUnit.text = "dimensionless"
                            precision = sub(att_type, 'precision')
                            precision.text = "0"
                            numericDomain = sub(att_type, 'numericDomain')
                            numberType = sub(numericDomain, 'numberType')
                            numberType.text = 'natural'


                        elif len(measure.findall('interval')) == 1:
                            measure.remove(measure.find('interval'))
                            att_type = sub(measure, 'interval')
                            unit = sub(att_type, 'unit')
                            standardUnit = sub(unit, 'standardUnit')
                            standardUnit.text = "dimensionless"
                            precision = sub(att_type, 'precision')
                            precision.text = "0"
                            numericDomain = sub(att_type, 'numericDomain')
                            numberType = sub(numericDomain, 'numberType')
                            numberType.text = 'natural'

    def add_attribute_traits(self, traitlist):
        '''
        Adds traits to the attributes contained in self.attributes as specified
        by the traitlist.  Traitlist is a list of tuples with each tuple
        containting two elements: the attribute name (string) and a dictionary
        of traits to be added to the attribute.  If the type of the trait
        ordinal, nothing will be changed.  Only traits with type interval will
        be appened too. 
        
        Parameters
        ----------
        traitlist : list
            A list of 2 element tuples where the first element contains a
            string and the second element conatins a dict. See example in
            docstring. The only keywords that are recognized are maximum,
            minimum, and precision.

        Example of traitlist:

            [('x', {'minimum' : '0', 'maximum' : '100'}), ('y', {'precision' :
            '0.1'})]
        
        '''

        for item in traitlist:
            for attribute in self.attributes:
                tree = ET.ElementTree(attribute)
                child = tree.findall('attributeName')[0]
                if child.text == item[0]:
                    #TODO:Cleaner way to do this than with if?
                    measure = tree.findall('measurementScale')[0]
                    if len(measure.findall('interval')) == 1:
                        interval = measure.findall('interval')[0]
                        for key in item[1].iterkeys():
                            if key == 'precision':
                                prec = interval.findall('precision')
                                if len(prec) == 0:
                                    precision = sub(interval, 'precision')
                                    precision.text = str(item[1][key])
                                elif len(prec) == 1:
                                    prec[0].text = str(item[1][key])
                            elif key == 'minimum':
                                numericDomain =\
                                           interval.findall('numericDomain')[0]
                                bnd = numericDomain.findall('bounds')
                                if len(bnd) == 0:
                                    bounds = sub(numericDomain, 'bounds')
                                    minimum = sub(bounds, 'minimum')
                                    minimum.attrib = {'exclusive' :
                                                                   'false'}
                                    minimum.text = str(item[1][key])
                                elif len(bnd) == 1:
                                    mins = bnd[0].findall('minimum')
                                    if len(mins) == 0:
                                        minimum = sub(bnd[0], 'minimum')
                                        minimum = sub(bnd[0], 'minimum')
                                        minimum.attrib = {'exclusive' :
                                                                       'false'}
                                        minimum.text = str(item[1][key])
                                    elif len(mins) == 1:
                                        bnd[0].remove(mins[0])
                                        minimum = sub(bnd[0], 'minimum')
                                        minimum.attrib = {'exclusive' :
                                                                   'false'}
                                        minimum.text = str(item[1][key])
                            elif key == 'maximum':
                                numericDomain =\
                                    interval.findall('numericDomain')[0]
                                bnd = numericDomain.findall('bounds')
                                if len(bnd) == 0:
                                    bounds = sub(numericDomain, 'bounds')
                                    maximum = sub(bounds, 'maximum')
                                    maximum.attrib = {'exclusive' :
                                                                   'false'}
                                    maximum.text = str(item[1][key])
                                elif len(bnd) == 1:
                                    maxs = bnd[0].findall('maximum')
                                    if len(maxs) == 0:
                                        maximum = sub(bnd[0], 'maximum')
                                        maximum.attrib = {'exclusive' :
                                                                   'false'}
                                        maximum.text = str(item[1][key])
                                    elif len(maxs) == 1:
                                        bnd[0].remove(maxs[0])
                                        maximum = sub(bnd[0], 'maximum')
                                        maximum.attrib = {'exclusive' :
                                                                   'false'}
                                        maximum.text = str(item[1][key])



    def write_meta_data(self, name=None):
        '''
        Writes out the xml tree that is contained in self.root and saves and
        .xml file in the currect working directory under the given filename. If
        no name is given save the xml as the same name as the input file.

        
        '''
        
        tree = ET.ElementTree(self.root)
        if name == None:
            tree.write(self.filename + '.xml')
        else:
            tree.write(name + '.xml')

                






