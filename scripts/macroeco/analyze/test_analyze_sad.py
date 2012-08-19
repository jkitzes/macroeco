#!/usr/bin/python

'''
Unit tests for sad_analysis.py
'''

from __future__ import division
import unittest
import numpy as np
from sad_analysis import *
from empirical import *
from predict_sad import *
import os

class TestSadAnalysis(unittest.TestCase):
    '''Test functions within sad_analysis.py'''

    def setUp(self):
        self.xyfile1 = open('xyfile1.csv', 'w')
        self.xyfile1.write('''spp_code, x, y, count
                        0, 0, 0, 1
                        1, 0, 0, 1
                        2, 0, 0, 0
                        3, 0, 0, 3
                        0, 0, 1, 0
                        1, 0, 1, 4
                        2, 0, 1, 0
                        3, 0, 1, 1
                        0, 1, 0, 1
                        1, 1, 0, 0
                        2, 1, 0, 3
                        3, 1, 0, 1
                        0, 1, 1, 0
                        1, 1, 1, 1
                        2, 1, 1, 3
                        3, 1, 1, 1
                        0, 2, 0, 0
                        1, 2, 0, 0
                        2, 2, 0, 2
                        3, 2, 0, 4
                        0, 2, 1, 0
                        1, 2, 1, 0
                        2, 2, 1, 0
                        3, 2, 1, 1''')
        self.xyfile1.close()
        self.xymeta = open('xyfile1.xml', 'w')

        #THIS METADATA IS JUST FOR TESTING, THE INFORMATION MEANS NOTHING
        self.xymeta.write('''<?xml version="1.0" encoding="UTF-8"?>
<eml:eml packageId="macroeco.113.1" system="knb" xmlns:eml="eml://ecoinformatics.org/eml-2.1.0" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:schemaLocation="eml://ecoinformatics.org/eml-2.1.0 eml.xsd"> <dataset> <title>Crosstimbers data set</title>
 <creator id="1336688248976"> <individualName><givenName>Ethan</givenName>
<surName>White</surName>
</individualName>
</creator>
 <abstract><para>This data contains a tree census on a 200 m by 200 m plot.  Collected in Oklahoma?</para>
</abstract>
<keywordSet><keyword>macroecology</keyword>
<keyword>crosstimbers</keyword>
</keywordSet>
<intellectualRights><para>This data is not to be distributed</para>
</intellectualRights>
<coverage><geographicCoverage><geographicDescription>Oklahoma</geographicDescription>
<boundingCoordinates><westBoundingCoordinate>-96.8461</westBoundingCoordinate>
<eastBoundingCoordinate>-96.8461</eastBoundingCoordinate>
<northBoundingCoordinate>35.5608</northBoundingCoordinate>
<southBoundingCoordinate>35.5608</southBoundingCoordinate>
</boundingCoordinates>
</geographicCoverage>
<temporalCoverage><singleDateTime><calendarDate>1998</calendarDate>
</singleDateTime>
</temporalCoverage>
</coverage>
<contact id="1336688260131"><individualName><givenName>Justin</givenName>
<surName>Kitzes</surName>
</individualName>
</contact>
<methods><methodStep><description><section><para>All trees were censused in a 200m by 200m plot. Not much is known about this dataset as of May, 10th 2012</para>
</section>
</description>
</methodStep>
</methods>
<dataTable id="1336688731568"><entityName>CRTI_1998.csv</entityName>
<physical id="1336688572976"><objectName>CRTI_1998.csv</objectName>
<size unit="byte">105393</size>
<dataFormat><textFormat><numHeaderLines>1</numHeaderLines>
<recordDelimiter>#x0A</recordDelimiter>
<attributeOrientation>column</attributeOrientation>
<simpleDelimited><fieldDelimiter>,</fieldDelimiter>
</simpleDelimited>
</textFormat>
</dataFormat>
<distribution><online><url>ecogrid://knb/macroeco.112.1</url>
</online>
</distribution>
</physical>
<attributeList><attribute id="1336688731569"><attributeName>spp_code</attributeName>
<attributeDefinition>An integer code that is unique to each species.  Integer codes are defined in CRTI_species_codes.csv</attributeDefinition>
<measurementScale><nominal><nonNumericDomain><textDomain><definition>See definition</definition>
</textDomain>
</nonNumericDomain>
</nominal>
</measurementScale>
</attribute>
<attribute id="1336688731570"><attributeName>x</attributeName>
<attributeDefinition>X coordinate within the plot, ranges from 0 to 199.9m.</attributeDefinition>
<measurementScale><interval><unit><standardUnit>meter</standardUnit>
</unit>
<precision>1</precision>
<numericDomain><numberType>real</numberType>
<bounds><minimum exclusive="false">0</minimum>
<maximum exclusive="false">2</maximum>
</bounds>
</numericDomain>
</interval>
</measurementScale>
</attribute>
<attribute id="1336688731571"><attributeName>y</attributeName>
<attributeDefinition>Y coordinate within the plot, ranges from 0 to 199.9m.</attributeDefinition>
<measurementScale><interval><unit><standardUnit>meter</standardUnit>
</unit>
<precision>1</precision>
<numericDomain><numberType>real</numberType>
<bounds><minimum exclusive="false">0</minimum>
<maximum exclusive="false">1</maximum>
</bounds>
</numericDomain>
</interval>
</measurementScale>
</attribute>
</attributeList>
<numberOfRecords>7626</numberOfRecords>
</dataTable>
</dataset>
 </eml:eml>

        ''')
        
        self.xymeta.close()
        self.xyfile2 = open('xyfile2.csv', 'w')
        self.xyfile2.write('''spp_code, x, y, count
                        0, 0, 0, 1
                        1, 0, 0, 1
                        2, 0, 0, 0
                        3, 0, 0, 3
                        4, 0, 0, 0
                        0, 0, 1, 0
                        1, 0, 1, 4
                        2, 0, 1, 0
                        3, 0, 1, 1
                        4, 0, 1, 0
                        0, 1, 0, 1
                        1, 1, 0, 0
                        2, 1, 0, 3
                        3, 1, 0, 1
                        4, 1, 0, 0
                        0, 1, 1, 0
                        1, 1, 1, 1
                        2, 1, 1, 3
                        3, 1, 1, 1
                        4, 1, 1, 0
                        0, 2, 0, 0
                        1, 2, 0, 0
                        2, 2, 0, 2
                        3, 2, 0, 4
                        4, 2, 0, 0
                        0, 2, 1, 0
                        1, 2, 1, 0
                        2, 2, 1, 0
                        3, 2, 1, 1
                        4, 2, 1, 0''')
        self.xyfile2.close()
        self.xymeta2 = open('xyfile2.xml', 'w')

        #THIS METADATA IS JUST FOR TESTING, THE INFORMATION MEANS NOTHING
        self.xymeta2.write('''<?xml version="1.0" encoding="UTF-8"?>
<eml:eml packageId="macroeco.113.1" system="knb" xmlns:eml="eml://ecoinformatics.org/eml-2.1.0" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:schemaLocation="eml://ecoinformatics.org/eml-2.1.0 eml.xsd"> <dataset> <title>Crosstimbers data set</title>
 <creator id="1336688248976"> <individualName><givenName>Ethan</givenName>
<surName>White</surName>
</individualName>
</creator>
 <abstract><para>This data contains a tree census on a 200 m by 200 m plot.  Collected in Oklahoma?</para>
</abstract>
<keywordSet><keyword>macroecology</keyword>
<keyword>crosstimbers</keyword>
</keywordSet>
<intellectualRights><para>This data is not to be distributed</para>
</intellectualRights>
<coverage><geographicCoverage><geographicDescription>Oklahoma</geographicDescription>
<boundingCoordinates><westBoundingCoordinate>-96.8461</westBoundingCoordinate>
<eastBoundingCoordinate>-96.8461</eastBoundingCoordinate>
<northBoundingCoordinate>35.5608</northBoundingCoordinate>
<southBoundingCoordinate>35.5608</southBoundingCoordinate>
</boundingCoordinates>
</geographicCoverage>
<temporalCoverage><singleDateTime><calendarDate>1998</calendarDate>
</singleDateTime>
</temporalCoverage>
</coverage>
<contact id="1336688260131"><individualName><givenName>Justin</givenName>
<surName>Kitzes</surName>
</individualName>
</contact>
<methods><methodStep><description><section><para>All trees were censused in a 200m by 200m plot. Not much is known about this dataset as of May, 10th 2012</para>
</section>
</description>
</methodStep>
</methods>
<dataTable id="1336688731568"><entityName>CRTI_1998.csv</entityName>
<physical id="1336688572976"><objectName>CRTI_1998.csv</objectName>
<size unit="byte">105393</size>
<dataFormat><textFormat><numHeaderLines>1</numHeaderLines>
<recordDelimiter>#x0A</recordDelimiter>
<attributeOrientation>column</attributeOrientation>
<simpleDelimited><fieldDelimiter>,</fieldDelimiter>
</simpleDelimited>
</textFormat>
</dataFormat>
<distribution><online><url>ecogrid://knb/macroeco.112.1</url>
</online>
</distribution>
</physical>
<attributeList><attribute id="1336688731569"><attributeName>spp_code</attributeName>
<attributeDefinition>An integer code that is unique to each species.  Integer codes are defined in CRTI_species_codes.csv</attributeDefinition>
<measurementScale><nominal><nonNumericDomain><textDomain><definition>See definition</definition>
</textDomain>
</nonNumericDomain>
</nominal>
</measurementScale>
</attribute>
<attribute id="1336688731570"><attributeName>x</attributeName>
<attributeDefinition>X coordinate within the plot, ranges from 0 to 199.9m.</attributeDefinition>
<measurementScale><interval><unit><standardUnit>meter</standardUnit>
</unit>
<precision>1</precision>
<numericDomain><numberType>real</numberType>
<bounds><minimum exclusive="false">0</minimum>
<maximum exclusive="false">2</maximum>
</bounds>
</numericDomain>
</interval>
</measurementScale>
</attribute>
<attribute id="1336688731571"><attributeName>y</attributeName>
<attributeDefinition>Y coordinate within the plot, ranges from 0 to 199.9m.</attributeDefinition>
<measurementScale><interval><unit><standardUnit>meter</standardUnit>
</unit>
<precision>1</precision>
<numericDomain><numberType>real</numberType>
<bounds><minimum exclusive="false">0</minimum>
<maximum exclusive="false">1</maximum>
</bounds>
</numericDomain>
</interval>
</measurementScale>
</attribute>
</attributeList>
<numberOfRecords>7626</numberOfRecords>
</dataTable>
</dataset>
 </eml:eml>

        ''')
        self.xymeta2.close()
        self.patch1 = Patch('xyfile1.csv')
                
    
    def tearDown(self):
        os.remove('xyfile1.csv')
        os.remove('xyfile1.xml')
        os.remove('xyfile2.csv')
        os.remove('xyfile2.xml')

    def test_get_gridded_sad_list(self):
        sad_list = get_gridded_sad_list('xyfile1.csv', [(1,1)])
        self.assertTrue(len(sad_list) == 1)
        self.assertTrue(len(sad_list[0]) == 1)
        self.assertTrue(len(sad_list[0][0]) == 4)
        self.assertTrue(np.array_equal(sad_list[0][0], np.array([2,6,8,11])))
        sad_list = get_gridded_sad_list('xyfile2.csv', [(1,1)])
        self.assertTrue(np.array_equal(sad_list[0][0], np.array([2,6,8,11,0])))
        sad_list = get_gridded_sad_list('xyfile2.csv', [(1,1)], clean=True)
        self.assertTrue(np.array_equal(sad_list[0][0], np.array([2,6,8,11])))
        sad_list = get_gridded_sad_list('xyfile2.csv', [(1,1), (3,1), (3,2)])
        self.assertTrue(len(sad_list) == 3)
        self.assertTrue(len(sad_list[0]) == 1 and len(sad_list[1]) == 3 and \
                        len(sad_list[2]) == 6)
        sd1 = np.array([1,1,0,3,0])
        sd2 = np.array([1,0,3,1,0])
        sd3 = np.array([0,0,2,4,0])
        sd4 = np.array([0,4,0,1,0])
        sd5 = np.array([0,1,3,1,0])
        sd6 = np.array([0,0,0,1,0])
        for sad in sad_list[2]:
            self.assertTrue(np.array_equal(sad, sd1) or\
                            np.array_equal(sad, sd2) or\
                            np.array_equal(sad, sd3) or\
                            np.array_equal(sad, sd4) or\
                            np.array_equal(sad, sd5) or\
                            np.array_equal(sad, sd6))
        sd1 = np.array([1,5,0,4,0])
        sd2 = np.array([1,1,6,2,0])
        sd3 = np.array([0,0,2,5,0])
        for sad in sad_list[1]:
            self.assertTrue(np.array_equal(sad, sd1) or\
                            np.array_equal(sad, sd2) or\
                            np.array_equal(sad, sd3))
        sad_list = get_gridded_sad_list('xyfile2.csv', [(1,1), (3,1), (3,2)], clean=True)
        self.assertTrue(len(sad_list) == 3)
        self.assertTrue(len(sad_list[0]) == 1 and len(sad_list[1]) == 3 and \
                        len(sad_list[2]) == 6)
        sd1 = np.array([1,1,3])
        sd2 = np.array([1,3,1])
        sd3 = np.array([2,4])
        sd4 = np.array([4,1])
        sd5 = np.array([1,3,1])
        sd6 = np.array([1])
        for sad in sad_list[2]:
            self.assertTrue(np.array_equal(sad, sd1) or\
                            np.array_equal(sad, sd2) or\
                            np.array_equal(sad, sd3) or\
                            np.array_equal(sad, sd4) or\
                            np.array_equal(sad, sd5) or\
                            np.array_equal(sad, sd6))
        sd1 = np.array([1,5,4])
        sd2 = np.array([1,1,6,2])
        sd3 = np.array([2,5])
        for sad in sad_list[1]:
            self.assertTrue(np.array_equal(sad, sd1) or\
                            np.array_equal(sad, sd2) or\
                            np.array_equal(sad, sd3))

    def test_get_obs_cdf_values(self):
        sad = get_gridded_sad_list('xyfile1.csv', [(1,1)])[0][0]
        cdf = get_obs_cdf_values(sad)
        self.assertTrue(np.array_equal(cdf, np.array([1/4, 2/4, 3/4, 1])))
        cdf = get_obs_cdf_values([1,1,2,2,5])
        self.assertTrue(np.array_equal(np.array([2/5, 2/5, 4/5, 4/5, 5/5]), cdf))
        self.assertRaises(ValueError, get_obs_cdf_values, [0,0,1,1,2])
        cdf = get_obs_cdf_values([5,2,1,1,2])
        self.assertTrue(np.array_equal(np.array([2/5, 2/5, 4/5, 4/5, 5/5]), cdf))


    #This function uses two already tested functions
    def test_get_obs_vs_pred_cdf(self):
        obs_pred = get_obs_vs_pred_cdf([1,1,2,2,5], 'mete')
        self.assertTrue(len(obs_pred) == 5)
        self.assertTrue(np.array_equal(obs_pred['n'], np.array([1,1,2,2,5])))
        self.assertRaises(ValueError, get_obs_vs_pred_cdf, [0,0,1,2,3], 'mete')
        
    def test_get_obs_and_pred_rarity(self):
        rarity = get_obs_and_pred_rarity([1,1,1,1,5,6,23,45], 'mete')
        self.assertTrue(rarity[0] == 6)
        self.assertRaises(ValueError, get_obs_and_pred_rarity, [0,0,1,2,4], 'mete')

    def test_get_obs_and_pred_Nmax(self):
        Nmax = get_obs_and_pred_Nmax([1,1,1,1,15], 'mete')
        self.assertTrue(Nmax[0] == 15)
        Nmax = get_obs_and_pred_Nmax([1,1,2,4,45,102], 'mete')
        self.assertTrue(Nmax[0] == 102)
        self.assertRaises(ValueError, get_obs_and_pred_Nmax, [0,1,2,4], 'mete')

    #No new functions are being used here
    def test_get_obs_pred_abund(self):
        pass

    def self_get_values_for_sad(self):
        self.asertRaises(ValueError, get_values_for_sad, [1,1,1,1,0,4,6,7,9,23,45,67], 'mete')
        

if __name__ == '__main__':
    unittest.main()
        
        



