#!/usr/bin/python

import os
import unittest
from macroeco.utils import workflow


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
        self.cwd = os.getcwd() + '/'

        self.pf = open('parameters.xml', 'w')
        self.pf.write("""<?xml version='1.0'?>
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
        pass
        os.remove(workflow.paramfile)

    def test_emptyask(self):
        pa = workflow.Parameters('nonexistantrun', None, {})
        self.assertEqual(pa.params, {})
        self.assertEqual(pa.interactive, False)

    def test_NIrunExists(self):
        req_params = {'size': 'descripsize', 'species': 'descripspp'}
        pa = workflow.Parameters('RunExists', None, req_params)
        self.assertTrue(len(pa.params) == 1)
        self.assertTrue(set(req_params.keys()).issubset(\
            set(pa.params['ParamfileTestCase'].keys())))
        self.assertTrue(pa.interactive == False)

        run = pa.params['ParamfileTestCase']
        self.assertTrue(run['size']*run['layers'][1] == 3*4.4)

    def test_MultipleNIRunsExist(self):
        pa = workflow.Parameters('ManyNIRuns', None, {})
        self.assertEqual(len(pa.params), 2)
        self.assertEqual(pa.params['FirstCase']['size'], 4.4)
        self.assertEqual(pa.params['FirstCase']['species'], 'E. coli')
        self.assertEqual(pa.params['FirstCase']['layers'], [0,3,12])
        self.assertEqual(pa.params['SecondCase']['size'], 2.2)
        self.assertEqual(pa.params['SecondCase']['species'], 'H. sapiens')
        self.assertEqual(pa.params['SecondCase']['layers'], [5])

    def test_UnnamedRunErrors(self):
        pa = workflow.Parameters('Unnamed', None, {})
        self.assertEqual(len(pa.params), 2)
        self.assertEqual(pa.params['run1']['size'], 4.4)
        self.assertEqual(pa.params['run1']['species'], 'E. coli')
        self.assertEqual(pa.params['run1']['layers'], [0,3,12])
        self.assertEqual(pa.params['run2']['size'], 2.2)
        self.assertEqual(pa.params['run2']['species'], 'H. sapiens')
        self.assertEqual(pa.params['run2']['layers'], [5])

    def test_InteractiveRun(self):
        pa = workflow.Parameters('Interactive', None, {})
        self.assertTrue(pa.interactive == True)

