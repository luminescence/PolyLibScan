import PolyLibScan.Analysis as lmp_lys
import pathlib2 as pl

import numpy as np

import unittest as ut

local_path = pl.Path(__file__).absolute().parent

class TestEpitopsy(ut.TestCase):

	def __init__(self, *args, **kwargs):
		super(TestEpitopsy, self).__init__(*args, **kwargs)
		self.project = lmp_lys.Project(local_path.joinpath('data/').as_posix(), 
                experimental_data=local_path.joinpath('data/static', 'ic50.h5').as_posix(),
                parameters=local_path.joinpath('data/static', 'parameters_hp.yml').as_posix(),
                protein_path='asd.pdb')
		self.job = self.project.jobs[0]

	def test_calc_protein_box(self):
		# test_margin
		box1 = self.job._calc_protein_box(margin=0)
		box2 = self.job._calc_protein_box(margin=10)

		box_diff = box2-box1
		print box_diff
		self.assertTrue(np.all((box_diff[0] + 10 ) == np.array([.0,.0,.0])))
		self.assertTrue(np.all((box_diff[1] - 10 ) == np.array([.0,.0,.0])))
		

	def test_add_run_to_epitopsy_map(self):
		pass

	def test_add_run_to_epitopsy_map(self):
		pass

	def test_create_epitopsy_map(self):
		pass

if __name__ == '__main__':
    ut.main(verbosity=2)
