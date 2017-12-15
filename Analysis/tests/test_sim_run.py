import PolyLibScan.Analysis as Analysis
import pathlib2 as pl
import mock
import numpy as np

import unittest as ut

local_path = pl.Path(__file__).absolute().parent

class TestRun(ut.TestCase):

    def __init__(self, *args, **kwargs):
        super(TestRun, self).__init__(*args, **kwargs)
        self.project = Analysis.Project(local_path.joinpath('data/').as_posix(), 
                experimental_data=local_path.joinpath('data/static', 'ic50.h5').as_posix(),
                parameters=local_path.joinpath('data/static', 'parameters_hp.yml').as_posix(),
                protein_path='asd.pdb')
        self.run_ = self.project.jobs[0][0]

    def test_distance(self):
        self.assertTrue(self.run_.distance>0)

    def test_meta(self):
        results = self.run_.meta()
        for key, value in results.iteritems():
            self.assertTrue(isinstance(value, basestring))

    def test_weights(self):
        self.assertTrue(self.run_.weights, dict)
        for key,val in self.run_.weights().items():
            self.assertTrue(isinstance(key, basestring))

    def test_sequence(self):
        self.assertTrue(len(self.run_.sequence()) > 0)

    def test_coorinates(self):
        results = self.run_.coordinates()
        self.assertTrue('start' in results)
        self.assertTrue('end' in results)

    def test_distance_to_active_site(self):
        results = self.run_._distance_to_active_site()
        self.assertTrue(isinstance(results, float))
        self.assertTrue(results > 0)

    def test_dist(self):
        type_ = [('x', np.float), ('y', np.float), ('z', np.float)]
        a = np.array((2.0, 3.0, 4.0), dtype=type_)
        b = np.array((0.0, 0.0, 0.0), dtype=type_)
        results = self.run_._dist(a, b)
        self.assertTrue(results > 0)

    def test_distance_time_series(self):
        results = self.run_.distance_time_series()
        self.assertTrue(len(results) > 0)
        self.assertEqual(results[0][0], 0)
        self.assertAlmostEqual(results[0][1], 11.44684301)
        self.assertEqual(results[-1][0], 200000)
        self.assertAlmostEqual(results[-1][1], 8.68289875388)

    def test_polymer_trajectory(self):
        results = self.run_.trajectory(molecule='polymer')
        self.assertTrue(len(results.next()) > 0)
        #self.assertTrue(results.next()

if __name__ == '__main__':
    ut.main(verbosity=2)
