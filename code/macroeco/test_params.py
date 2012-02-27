import unittest
import params


# Cases for future testing:
# params file has one interactive run, user changes values.
# params file has plural interactive runs (complicated dialog ahoy).
# No params file. Dialog, write, reload.
# Params file doesn't match ask. Dialog, write, reload, check against ask.
## params.xml proper subset of ask
## Neither a proper subset of the other
# Types of param: string, int, float, lists of any of those; mixed-type list (ick).


class ParamfileTestCase(unittest.TestCase):

    def setUp(self):
        self.asklist = {'size':'positive number', 'species':'string','layers':'list of positive integers'}
        self.pf = open(params.paramfile, 'w')
        self.pf.write("""<?xml version='1.0'?>
<METE>
    <analysis scriptname='RunExists' version='0.5'>
        <run name='ParamfileTestCase'>
            <interactive>F</interactive>
            <param name='size' value='4.4' />
            <param name='species' value='E. coli' />
            <param name='layers' value='[0,3,12]' />
        </run>
    </analysis>
    <analysis scriptname='ManyNIRuns' version='0.5'>
        <run name='FirstCase' >
            <interactive>F</interactive>
            <param name='size' value='4.4' />
            <param name='species' value='E. coli' />
            <param name='layers' value='[0,3,12]' />
        </run>
        <run name='SecondCase' >
            <interactive>F</interactive>
            <param name='size' value='2.2' />
            <param name='species' value='H. sapiens' />
            <param name='layers' value='[5]' />
        </run>
    </analysis>
    <analysis scriptname='Unnamed' version='0.5'>
        <run>
            <interactive>F</interactive>
            <param name='size' value='4.4' />
            <param name='species' value='E. coli' />
            <param name='layers' value='[0,3,12]' />
        </run>
        <run>
            <interactive>F</interactive>
            <param name='size' value='2.2' />
            <param name='species' value='H. sapiens' />
            <param name='layers' value='[5]' />
        </run>
    </analysis>
</METE>""")
        self.pf.close()

    def test_emptyask(self):
        pa = params.Parameters('unittest', {})
        self.assertEqual(pa.params, {})
        self.assertEqual(pa.interactive, None)

    def test_NIrunExists(self):
        pa = params.Parameters('RunExists', self.asklist)
        self.assertTrue(len(pa.params) == 1)
        self.assertTrue(set(self.asklist.keys()).issubset(set(pa.params['ParamfileTestCase'].keys())))
        self.assertTrue(pa.interactive == False)
        run = pa.params['ParamfileTestCase']
        convertlist = map(int, run['layers'][1:-1].split(','))
        self.assertTrue(float(run['size'])*convertlist[1] == 3*4.4)

    def test_MultipleNIRunsExist(self):
        pa = params.Parameters('ManyNIRuns', self.asklist)
        self.assertTrue(len(pa.params) == 2)
        self.assertTrue(pa.params['FirstCase']['size'] == '4.4')
        self.assertTrue(pa.params['FirstCase']['species'] == 'E. coli')
        self.assertTrue(pa.params['FirstCase']['layers'] == '[0,3,12]')
        self.assertTrue(pa.params['SecondCase']['size'] == '2.2')
        self.assertTrue(pa.params['SecondCase']['species'] == 'H. sapiens')
        self.assertTrue(pa.params['SecondCase']['layers'] == '[5]')

    def test_MultipleUnnamedRunsExist(self):
        pa = params.Parameters('Unnamed', self.asklist)
        self.assertTrue(len(pa.params) == 2)
        self.assertTrue(pa.params['autoname0']['size'] == '4.4')
        self.assertTrue(pa.params['autoname0']['species'] == 'E. coli')
        self.assertTrue(pa.params['autoname0']['layers'] == '[0,3,12]')
        self.assertTrue(pa.params['autoname1']['size'] == '2.2')
        self.assertTrue(pa.params['autoname1']['species'] == 'H. sapiens')
        self.assertTrue(pa.params['autoname1']['layers'] == '[5]')



if __name__ == '__main__':
    unittest.main()
