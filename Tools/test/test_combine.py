from PolyLibScan.Tools.polymer2lmp import PolymerCreator
from PolyLibScan.Tools.pdb2lmp import ProteinCreator
from PolyLibScan.Tools.lmp_combine import EnvManipulator
from PolyLibScan.Tools.environment import Environment
import pathlib2 as pl

import numpy as np

import unittest as ut

local_path = pl.Path(__file__).absolute().parent

class TestCombine_Maker(ut.TestCase):

	def __init__(self, *args, **kwargs):
		super(TestCombine_Maker, self).__init__(*args, **kwargs)
		env = Environment(local_path.joinpath('data', 'updated_parameters8.yml').as_posix())
		self.protein = ProteinCreator(env, local_path.joinpath('data/1LYA.pdb').as_posix())
		self.polymer = PolymerCreator(env, ['ED', 'BP', 'CBS'], length=20, mode='random')

	def test_random_start_position(self):
		env = Environment(local_path.joinpath('data', 'updated_parameters8.yml').as_posix())
		self.protein.env = env
		self.polymer.env = env
		self.protein.create()
		self.polymer.create()

		polymers = filter(lambda x:x.mol_type=='polymer', env.molecules.values())
		env_manager = EnvManipulator(env, randomize_positions=True)
		box_before = polymers[0].box.copy()
		#print box_before
		env_manager.create_random_start_positions()
		#env.molecules[keys[1]].move(np.array([1000.0, 0.0, 0.0]))
		box_after = polymers[0].box.copy()
		#print box_after
		self.assertFalse(np.array_equal(box_before, box_after))

	def test_overlap(self):
		# Setup system
		env = Environment(local_path.joinpath('data', 'updated_parameters8.yml').as_posix())
		self.protein.env = env
		self.polymer.env = env
		self.protein.create()
		self.polymer.create()
		env_manager = EnvManipulator(env)
		
		keys = env.molecules.keys()
		center1 = (env.molecules[keys[0]].box[:, 0] + env.molecules[keys[0]].box[:, 1])/2.0
		center2 = (env.molecules[keys[1]].box[:, 0] + env.molecules[keys[1]].box[:, 1])/2.0
		
		dist = center2 - center1
		print center1, center2
		# move first molecule to the second - creating an overlap
		env.molecules[keys[0]].move(dist)
		center1 = (env.molecules[keys[0]].box[:, 0] + env.molecules[keys[0]].box[:, 1])/2.0
		center2 = (env.molecules[keys[1]].box[:, 0] + env.molecules[keys[1]].box[:, 1])/2.0
		
		dist = center2 - center1
		print center1, center2
		self.assertTrue(env_manager.molecules_within_margin(env.molecules[keys[0]].data['particles'],
                                                            env.molecules[keys[1]].data['particles']))

		# move one molecule far, far away
		env.molecules[keys[0]].move(np.array([1000.0, 0.0, 0.0]))
		self.assertFalse(env_manager.molecules_within_margin(env.molecules[keys[0]].data['particles'],
                                                             env.molecules[keys[1]].data['particles']))

	def test_random_start_position_distance(self):
		env = Environment(local_path.joinpath('data', 'updated_parameters8.yml').as_posix())
		self.protein.env = env
		self.polymer.env = env
		self.protein.create()
		self.polymer.create()
		env_manager = EnvManipulator(env)
		env_manager.create_random_start_positions()

		polymer = filter(lambda x:x.mol_type=='polymer', env.molecules.values())[0]
		protein = filter(lambda x:x.mol_type=='protein', env.molecules.values())[0]
		polymer_radius = np.linalg.norm(polymer.box[:,0] - polymer.box[:,1])/2.0
		protein_radius = np.linalg.norm(protein.box[:,0] - protein.box[:,1])/2.0
		polymer_center = (polymer.box[:,0] + polymer.box[:,1])/2.0
		protein_center = (protein.box[:,0] + protein.box[:,1])/2.0
		distance = np.linalg.norm( polymer_center-protein_center )

		self.assertLess(distance, polymer_radius+protein_radius)

	# def test_affinities(self):
	# 	env = Environment(local_path.joinpath('data', 'updated_parameters8.yml').as_posix())
	# 	self.protein.env = env
	# 	self.polymer.env = env
	# 	self.protein.create()
	# 	self.polymer.create()
	# 	env_manager = EnvManipulator(env)

	# 	mols = env.molecules
	# 	affinities = env.ff['affinity']

	# 	number_of_atomsTypes_before = len(env.atom_type)
	# 	env_manager.add_affinities(affinity_file=local_path.joinpath('data/affinities_Originals2.h5').as_posix())
	# 	number_of_atomsTypes_after = len(env.atom_type)

	# 	self.assertNotEqual(number_of_atomsTypes_before, number_of_atomsTypes_after)

	# 	sorted_atypes = sorted(env.atom_type.values(), key=lambda x:x.Id, reverse=True)
	# 	for i, a_type1 in enumerate(sorted_atypes):
	# 		for a_type2 in sorted_atypes:
	# 			self.assertTrue((a_type1, a_type2) in affinities)

	def test_make_particle_unique(self):
		env = Environment(local_path.joinpath('data', 'updated_parameters8.yml').as_posix())
		self.protein.env = env
		self.polymer.env = env
		self.protein.create()
		env_manager = EnvManipulator(env)

		molecule = env.molecules.values()[0]
		particle = molecule.data['particles'][0]
		old_type = particle.type_
		new_type = env_manager._make_particle_unique(particle)
		new_type = particle.type_
		self.assertTrue(new_type.name in env.atom_type)
		self.assertEqual(new_type.Id, particle.type_.Id)
		self.assertNotEqual(old_type.Id, new_type.Id)
		self.assertEqual(old_type.mass, new_type.mass)
		self.assertEqual(old_type.radius, new_type.radius)

if __name__ == '__main__':
    ut.main(verbosity=2)
