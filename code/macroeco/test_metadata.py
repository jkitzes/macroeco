import unittest
import metadata
from matplotlib.mlab import csv

later = '''class DummyMetadataString(unittest.TestCase):

    def setUp(self):
        self.asklist = {'gx/measurementScale/interval/precision':'1.1',
                        'gx/measurementScale/interval/numericDomain/bounds/minimum':'0',
                        'gx/measurementScale/interval/numericDomain/bounds/maximum':'399',
                        'gy/measurementScale/interval/precision':'10.5',
                        'gy/measurementScale/interval/numericDomain/bounds/minimum':'100',
                        'gy/measurementScale/interval/numericDomain/bounds/maximum':'4000'}
        self.xf = open('test_metadata.xml','w')
        self.xf.write('TODO: For metadata practice, First line will be ignored.\n')
        for key in self.asklist.keys():
            self.xf.write(key+ ','+str(self.asklist[key])+'\n')

        self.xf.close()

    def test_getFromComment(self):
        meta = metadata.Metadata('test_metadata.csv',('gx','gy'))
        for ask in self.asklist:
            assert meta.data[ask] == self.asklist[ask]'''

class MorphoMetadata(unittest.TestCase):

    def test_getfromBCIS_Arch(self):
        asklist = [('gx','precision'),
                        ('gy','precision'),
                        ('gx','minimum'),
                        ('gx','maximum'),
                        ('gy','maximum')]
        values = ['0.1', '0.1','0','999.9','499.9']
        expect = dict(zip(asklist, values))
        meta = metadata.Metadata('../../data/archival/BCIS/BCIS_1995.csv')
        meta.get_dataTable_values(asklist)
        for ask in asklist:
            assert meta.TableDescr[ask] == expect[ask]

if __name__=="__main__":
    unittest.main()
