import PolyLibScan.Analysis.density_container as dc
import pathlib2 as pl
import mock
import numpy as np

import unittest as ut

local_path = pl.Path(__file__).absolute().parent

def calc_protein_box(margin):
	box = np.array([[-10,-10,-10],[10,10,10]])
	box[0] -= 10
	box[1] += 10
	return box

def create_sim():
	sim = mock.Mock()
	sim.particle_ids = {'polymer': [1,2,3]}
	sim._calc_protein_box = calc_protein_box
	sim._parse = mock.Mock()
	attrs = {'load_traj_type_order.return_value': [1,2,3,4,5,1]}
	sim._parse.configure_mock(**attrs)
	return sim

class TestEpitopsy(ut.TestCase):

	def __init__(self, *args, **kwargs):
		super(TestEpitopsy, self).__init__(*args, **kwargs)
		pass

	def test_equal(self):
		sim = create_sim()
		map1 = dc.DensityContainer(sim, monomer_id=1)
		map2 = dc.DensityContainer(sim, monomer_id=1)
		self.assertTrue(map1==map2)

	def test_unequal(self):
		sim = create_sim()
		map1 = dc.DensityContainer(sim, monomer_id=1)
		map2 = dc.DensityContainer(sim, monomer_id=2)
		self.assertTrue(map1!=map2)

	def test_monomer_id(self):
		sim = create_sim()
		map1 = dc.DensityContainer(sim, monomer_id=1)
		self.assertEqual(map1.monomer_id, [1])
		map1.monomer_id = 2
		self.assertEqual(map1.monomer_id, [2])
		map1.monomer_id = [1,2,3]
		self.assertEqual(map1.monomer_id, [1,2,3])
		map1.monomer_id = np.array([1,2,3])
		self.assertEqual(map1.monomer_id, [1,2,3])
		map1.monomer_id = 'all'
		self.assertEqual(map1.monomer_id, [1,2,3])
		
		with self.assertRaises(ValueError):
			map1.monomer_id = 'something else'
		with self.assertRaises(ValueError):
			map1.monomer_id = 4

	def test_resolution(self):
		sim = create_sim()
		map1 = dc.DensityContainer(sim, monomer_id=1, resolution=10.0)
		self.assertEqual(map1.resolution, 10.0)
		with self.assertRaises(UserWarning):
			map1.resolution = 333

	def test_margin(self):
		sim = create_sim()
		map1 = dc.DensityContainer(sim, monomer_id=1, margin=20.0)
		self.assertEqual(map1.margin, 20.0)
		with self.assertRaises(UserWarning):
			map1.margin = 12

	# def test_add(self):
	# 	sim = create_sim()
	# 	map1 = dc.DensityContainer(sim, monomer_id=1)
	# 	map1.map = np.ones((2,2))
	# 	map2 = dc.DensityContainer(sim, monomer_id=1)
	# 	map2.map = np.ones((2,2))
	# 	map3 = map1 + map2
	# 	self.assertTrue(np.array_equal(map3.map, 2*np.ones((2,2))))

	# def test_iadd(self):
	# 	sim = create_sim()
	# 	map1 = dc.DensityContainer(sim, monomer_id=1)
	# 	map1.map = np.ones((2,2))
	# 	map2 = dc.DensityContainer(sim, monomer_id=1)
	# 	map2.map = np.ones((2,2))
	# 	map1 += map2
	# 	self.assertTrue(np.array_equal(map1.map, 2*np.ones((2,2))))

	def test_check_for_difference(self):
		sim = create_sim()
		map1 = dc.DensityContainer(sim, monomer_id=1)
		map2 = dc.DensityContainer(sim, monomer_id=1)
		self.assertEqual(map1._check_for_difference(map2), 'No differences found.')
		map1._resolution = 100
		self.assertEqual(map1._check_for_difference(map2), 'Resolution in maps differ: 100.0 and 3.0')
		map1._resolution = 3.0
		map1._margin = 15.0
		self.assertEqual(map1._check_for_difference(map2), 'Margins in maps differ: 15.0 and 20.0')
		map1._margin = 20.0
		map1._monomer_id = [2]
		self.assertEqual(map1._check_for_difference(map2), 'Monomers in maps differ: [2] and [1]')
		map1._monomer_id = [1]
		previous_box = map1.box
		map1.box = np.array([[1,1,1], [2,2,2]])
		self.assertEqual(map1._check_for_difference(map2), 'Boxes in maps differ: [[1 1 1]\n [2 2 2]] and [[-20. -20. -20.]\n [ 20.  20.  20.]]')

if __name__ == '__main__':
    ut.main(verbosity=2)
