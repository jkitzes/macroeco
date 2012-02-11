'''Python Script which takes in metadata.xml
and gets required values from the file

Date: 2/3/2011

'''

import os
import datetime

class NoAdditionalMetadata(Exception):
    '''This exception is thrown when the metadata
    file that the user included does not contain 
    the mdoule 'additionalMetadata'
    
    '''
    
    def __init__(self,value=None):
        Exception.__init__(self)
        self.value = value
        print "exception with %s at time %s" % (self.value,datetime.datetime.now())
    def __str__(self):
        return self.value


def open_metadata_xml(filename):
    '''This function will take in a file name as an arguement
    and save the data from the file to a variable and 
    return it
    '''
    
    with open(filename, 'r') as fin:
        metadata = fin.read()
    
    return metadata

def get_value_from_metadata(value_name,filename):
    '''This function accepts two strings value and
    filename.  value is the name of the required
    metadata value that is to be retrieved.
    filename is the name of the metadata file. 
    The function returns a string.
    
    More error checking needs to be done!
    
    '''
    
    metadata = open_metadata_xml(filename)
    additional_meta_start = metadata.find('<additionalMetadata>')
    if additional_meta_start == -1:
        raise NoAdditionalMetadata('%s does not contain the module \
                                    additionalMetaData' % (filename))
    start = metadata.find('!!' + value_name, additional_meta_start)
    if(start == -1):
        raise Exception('%s does not exist or is formatted improperly \
                         in %s' % (value_name,filename))
    end = metadata.find('@@', start + 1)
    value = metadata[start + 2 : end].strip().split('=')[1]
    return value
    
    
    
        
    
    
