#!/usr/bin/python

import os
import unittest
import workflow


# Cases for future testing:
# params file has one interactive run, user changes values.
# params file has plural interactive runs (complicated dialog ahoy).
# No params file. Dialog, write, reload.
# Params file doesn't match ask. Dialog, write, reload, check against ask.
## workflow.xml proper subset of ask
## Neither a proper subset of the other
# Types of param: string, int, float, lists of those; mixed-type list (ick).


class ParamfileTestCase(unittest.TestCase):

    def setUp(self):
        self.pf = open('parameters.xml', 'w')
        self.pf.write("""
<?xml version='1.0'?>
<macroeco>
    <analysis script_name='RunExists' version='0.5' interactive='F'>
        <run name='ParamfileTestCase'>
            <param name='size' value='4.4' />
            <param name='species' value='E. coli' />
            <param name='layers' value='[0,3,12]' />
        </run>
    </analysis>
    <analysis script_name='ManyNIRuns' version='0.5' interactive='f'>
        <run name='FirstCase' >
            <param name='size' value='4.4' />
            <param name='species' value='E. coli' />
            <param name='layers' value='[0,3,12]' />
        </run>
        <run name='SecondCase' >
            <param name='size' value='2.2' />
            <param name='species' value='H. sapiens' />
            <param name='layers' value='[5]' />
        </run>
    </analysis>
    <analysis script_name='Unnamed' version='0.5'>
        <run>
            <param name='size' value='4.4' />
            <param name='species' value='E. coli' />
            <param name='layers' value='[0,3,12]' />
        </run>
        <run>
            <param name='size' value='2.2' />
            <param name='species' value='H. sapiens' />
            <param name='layers' value='[5]' />
        </run>
    </analysis>
    <analysis script_name='Interactive' version='0.5' interactive='T'>
        <run>
            <param name='size' value='4.4' />
            <param name='species' value='E. coli' />
            <param name='layers' value='[0,3,12]' />
        </run>
        <run>
            <param name='size' value='2.2' />
            <param name='species' value='H. sapiens' />
            <param name='layers' value='[5]' />
        </run>
    </analysis>
    
</macroeco>
                      """)
        self.pf.close()

    def tearDown(self):
        os.remove(workflow.paramfile)

    def test_emptyask(self):
        pa = workflow.Parameters('unittest', {})
        self.assertEqual(pa.params, {})
        self.assertEqual(pa.interactive, None)

    def test_NIrunExists(self):
        pa = workflow.Parameters('RunExists', self.asklist)
        self.assertTrue(len(pa.params) == 1)
        self.assertTrue(set(self.asklist.keys()).issubset(set(pa.params['ParamfileTestCase'].keys())))
        self.assertTrue(pa.interactive == False)
        run = pa.params['ParamfileTestCase']
        convertlist = map(int, run['layers'][1:-1].split(','))
        self.assertTrue(float(run['size'])*convertlist[1] == 3*4.4)

    def test_MultipleNIRunsExist(self):
        pa = workflow.Parameters('ManyNIRuns', self.asklist)
        self.assertTrue(len(pa.params) == 2)
        self.assertTrue(pa.params['FirstCase']['size'] == '4.4')
        self.assertTrue(pa.params['FirstCase']['species'] == 'E. coli')
        self.assertTrue(pa.params['FirstCase']['layers'] == '[0,3,12]')
        self.assertTrue(pa.params['SecondCase']['size'] == '2.2')
        self.assertTrue(pa.params['SecondCase']['species'] == 'H. sapiens')
        self.assertTrue(pa.params['SecondCase']['layers'] == '[5]')

    def test_MultipleUnnamedRunsExist(self):
        pa = workflow.Parameters('Unnamed', self.asklist)
        self.assertTrue(len(pa.params) == 2)
        self.assertTrue(pa.params['autoname0']['size'] == '4.4')
        self.assertTrue(pa.params['autoname0']['species'] == 'E. coli')
        self.assertTrue(pa.params['autoname0']['layers'] == '[0,3,12]')
        self.assertTrue(pa.params['autoname1']['size'] == '2.2')
        self.assertTrue(pa.params['autoname1']['species'] == 'H. sapiens')
        self.assertTrue(pa.params['autoname1']['layers'] == '[5]')

    def test_InteractiveRun(self):
        pa = workflow.Parameters('Interactive', self.asklist)
        self.assertTrue(pa.interactive == True)

if __name__ == '__main__':
    unittest.main()
