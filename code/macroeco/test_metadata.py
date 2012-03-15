import unittest
import metadata
from matplotlib.mlab import csv


class MorphoMetadata(unittest.TestCase):
    def setUp(self):
        self.meta = metadata.Metadata('../../data/archival/BCIS/BCIS_1995.csv')        

    def test_tablefromBCIS_Archival(self):
        asklist = [('gx','precision'),
                        ('gy','precision'),
                        ('gx','minimum'),
                        ('gx','maximum'),
                        ('gy','maximum')]
        values = ['0.1', '0.1','0','999.9','499.9']
        expect = dict(zip(asklist, values))

        self.meta.get_dataTable_values(asklist)
        for ask in asklist:
            assert self.meta.TableDescr[ask] == expect[ask]

    def test_region(self):
        edges = self.meta.get_coverage_region()
        assert edges == [9.15,-79.85,9.15,-79.85]


if __name__=="__main__":
    unittest.main()
