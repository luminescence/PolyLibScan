from PolyLibScan.Tools.environment import Environment
import PolyLibScan.Tools.lmp_force_field as ff
import pathlib2 as pl

import mock
import unittest as ut

local_path = pl.Path(__file__).absolute().parent


class Atom_type(object):
	def __init__(self, id_, name):
		self.Id = id_
		self.name = name
		self.interacting = True

class TestForceField(ut.TestCase):

	def __init__(self, *args, **kwargs):
		super(TestForceField, self).__init__(*args, **kwargs)
		#Environment(local_path.join('data', 'newLammpsConfig.cfg').as_posix())
		#self.affinity = ff.Interaction(env, 'test')

	def create_new_env(self):
		env = mock.MagicMock()
		env.pair_type = {}
		env.pair_type['test'] = mock.MagicMock()
		env.pair_type['test'].parameters = {'coeffs': [1.2, 2.5]}
		env.pair_type['test'].kind = 'morse'
		env.atom_type = {'test_atom1': Atom_type(id_=1, name='test_atom1')}
		env.atom_type['test_atom2'] = Atom_type(id_=2, name='test_atom2')
		env.atom_type['test_atom3'] = Atom_type(id_=3, name='test_atom3')
		return env
	'''
	def test_initialization(self):
		env = self.create_new_env()
		affinity = ff.Interaction(env, env.pair_type['test'])
		
		self.assertTrue(set(affinity._matrix.keys()) == set([1,2,3]))
		self.assertTrue(len(affinity._matrix.values()) == 3)
		
		self.assertTrue(affinity[(env.atom_type['test_atom2'],env.atom_type['test_atom3'])].pair_type.kind == 'morse')
	'''
	def test_add(self):
		env = self.create_new_env()
		affinity = ff.Interaction(env, env.pair_type['test'])
		# without hook
		env.atom_type['test_atom4'] = Atom_type(id_=4, name='test_atom4')
		affinity.update(env.atom_type['test_atom4'])

		self.assertTrue(affinity[(env.atom_type['test_atom4'], env.atom_type['test_atom3'])].pair_type.kind == 'morse')


class TestPair(ut.TestCase):

	def __init__(self, *args, **kwargs):
		super(TestPair, self).__init__(*args, **kwargs)
		self.pair_type_lj = mock.Mock()
		self.pair_type_lj.kind = 'lj/cut'
		self.pair_type_lj.parameters = {'coeffs': [1,2,3]}
		self.pair_type_lj.cutoff = 10.0
		self.pair_type_soft = mock.Mock()
		self.pair_type_soft.kind = 'soft'
		self.pair_type_soft.parameters = {'coeffs': [4,5,6]}
		self.pair_type_soft.cutoff = 6.0
		self.pair_type_morse = mock.Mock()
		self.pair_type_morse.kind = 'morse'
		self.pair_type_morse.parameters = {'coeffs': [7,8,9]}
		self.pair_type_morse.cutoff = 15.0
		self.atom_type1 = mock.Mock()
		self.atom_type1.hydrophobicity = 0
		self.atom_type1.surface_energy = 0.0
		self.atom_type1.radius = 2.0
		self.atom_type2 = mock.Mock()
		self.atom_type2.hydrophobicity = 0
		self.atom_type2.surface_energy = 0.0
		self.atom_type2.radius = 2.0
		self.atom_type3 = mock.Mock()
		self.atom_type3.hydrophobicity = 1.0
		self.atom_type3.surface_energy = 2.0
		self.atom_type3.radius = 2.0
		self.atom_type4 = mock.Mock()
		self.atom_type4.hydrophobicity = 1.0
		self.atom_type4.surface_energy = 0.0
		self.atom_type4.radius = 2.0

	def test_lj_pair(self):
		test_pair = ff.Pair(self.pair_type_lj, self.atom_type1, self.atom_type2, epsilon=5)
		self.assertEqual(test_pair.epsilon, 5)
		test_pair = ff.Pair(self.pair_type_lj, self.atom_type1, self.atom_type2)
		self.assertEqual(test_pair.epsilon, 1)
		self.assertEqual(test_pair.alpha, None)
		self.assertAlmostEqual(test_pair.cutoff, 10.0)

	def test_morse_pair(self):
		test_pair = ff.Pair(self.pair_type_morse, self.atom_type1, self.atom_type2, epsilon=5, alpha=5.5)
		self.assertEqual(test_pair.epsilon, 5)
		self.assertEqual(test_pair.alpha, 5.5)
		self.assertEqual(test_pair.cutoff, 15)
		test_pair = ff.Pair(self.pair_type_morse, self.atom_type1, self.atom_type2)
		self.assertEqual(test_pair.epsilon, 7)
		self.assertEqual(test_pair.alpha, 8)
		self.assertAlmostEqual(test_pair.cutoff, 15)

	def test_soft_pair(self):
		test_pair = ff.Pair(self.pair_type_soft, self.atom_type1, self.atom_type2, epsilon=6, alpha=5.5)
		self.assertEqual(test_pair.epsilon, 6)
		self.assertEqual(test_pair.alpha, 5.5)
		self.assertAlmostEqual(test_pair.cutoff, 6)
		test_pair = ff.Pair(self.pair_type_soft, self.atom_type1, self.atom_type2)
		self.assertEqual(test_pair.epsilon, 0.0)
		self.assertEqual(test_pair.alpha, None)
		self.assertAlmostEqual(test_pair.cutoff, 6)
		test_pair = ff.Pair(self.pair_type_soft, self.atom_type3, self.atom_type4, epsilon=10.0)
		self.assertEqual(test_pair.epsilon, 10.0)
		self.assertEqual(test_pair.alpha, None)
		self.assertAlmostEqual(test_pair.cutoff, 6)
		test_pair = ff.Pair(self.pair_type_soft, self.atom_type3, self.atom_type4)
		self.assertAlmostEqual(test_pair.epsilon, 4.0)
		self.assertEqual(test_pair.alpha, None)
		self.assertAlmostEqual(test_pair.cutoff, 6.0)

if __name__ == '__main__':
    ut.main(verbosity=2)
