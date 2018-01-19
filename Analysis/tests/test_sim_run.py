import PolyLibScan.Analysis as Analysis
from PolyLibScan.Analysis.sim_run import AtomFilter
import pathlib2 as pl
import numpy as np

import unittest as ut

local_path = pl.Path(__file__).absolute().parent


def generate_project():
    return Analysis.Project(local_path.joinpath('data/').as_posix(),
                            experimental_data=local_path.joinpath('data/static', 'ic50.h5').as_posix(),
                            parameters=local_path.joinpath('data/static', 'parameters_hp.yml').as_posix(),
                            protein_path='asd.pdb')


class TestRun(ut.TestCase):

    def __init__(self, *args, **kwargs):
        super(TestRun, self).__init__(*args, **kwargs)
        self.project = generate_project()
        self.run_ = self.project.jobs[0][0]

    def test_distance(self):
        self.assertTrue(self.run_.distance > 0)

    def test_meta(self):
        results = self.run_.meta()
        for key, value in results.iteritems():
            self.assertTrue(isinstance(value, basestring))

    def test_weights(self):
        self.assertTrue(self.run_.weights, dict)
        for key, val in self.run_.weights().items():
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

    def test_trajectory(self):
        protein_traj = self.run_.trajectory(molecule='protein').next()
        polymer_traj = self.run_.trajectory(molecule='polymer').next()
        full_traj = self.run_.trajectory(molecule='full').next()

        self.assertTrue(len(polymer_traj) == 20)    # hard coded for given jobdata.h5

        # only works as long as there's only one protein and one polymer in jobdata.h5
        protein_and_polymer_traj = np.concatenate((protein_traj, polymer_traj))
        self.assertTrue(np.alltrue(protein_and_polymer_traj == full_traj))


class TestAtomFilter(ut.TestCase):

    def __init__(self, *args, **kwargs):
        super(TestAtomFilter, self).__init__(*args, **kwargs)
        self.project = generate_project()
        self.run_ = self.project.jobs[0][0]

    def test_mask(self):
        polymer_monomers = AtomFilter(self.run_.job.trajectory_order, sim=self.run_.job,
                                      monomer_id=self.run_.job.particle_ids['polymer'], molecule='polymer')
        protein_monomers = AtomFilter(self.run_.job.trajectory_order, sim=self.run_.job,
                                      monomer_id=self.run_.job.particle_ids['protein'], molecule='protein')
        all_monomers = AtomFilter(self.run_.job.trajectory_order, sim=self.run_.job,
                                  monomer_id=np.concatenate((self.run_.job.particle_ids['protein'],self.run_.job.particle_ids['polymer'])), molecule='full')
        self.assertEqual(polymer_monomers.mask.sum(), len(self.run_.sequence()))
        self.assertEqual(all_monomers.mask.sum(), polymer_monomers.mask.sum() + protein_monomers.mask.sum())

        # hard coded check if correct number of protein beads is filtered
        if 'PIS2' in self.project.protein_path.upper():
            self.assertEqual(protein_monomers.mask.sum(), 301)

    def test_filter_function(self):
        """reconstruct np.array from function and compare it to mask"""
        polymer_monomers = AtomFilter(self.run_.job.trajectory_order, sim=self.run_.job,
                                      monomer_id=self.run_.job.particle_ids['polymer'], molecule='full')

        bool_list = []
        for idx in enumerate(self.run_.job.trajectory_order):
            bool_list.append(polymer_monomers.filter_function())

        self.assertTrue(np.alltrue(polymer_monomers.mask == bool_list))


if __name__ == '__main__':
    ut.main(verbosity=2)
