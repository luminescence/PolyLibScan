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

	def test_whatever(self):
		pass

if __name__ == '__main__':
    ut.main(verbosity=2)
