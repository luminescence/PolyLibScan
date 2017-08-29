from PolyLibScan.Tools.lmpWriter import LmpWriter
import pathlib2 as pl

import numpy as np

import mock
import unittest as ut

local_path = pl.Path(__file__).absolute().parent

class TestWriter(ut.TestCase):

	def __init__(self, *args, **kwargs):
		super(TestWriter, self).__init__(*args, **kwargs)
		self.env = mock.Mock()
		self.env._calc_box = mock.MagicMock(return_value=np.array([[1,2],[1,2],[1,2]]))
		self.env.globals = {'atom_style': 'bond'}

	def test_entry_is_needed(self):
		writer = LmpWriter(self.env)
		
		self.env.globals.update({'atom_style': 'bond'})
		result1 = writer.entry_is_needed('bond')
		self.assertTrue(result1)

		self.env.globals.update({'atom_style': 'bond'})
		result1 = writer.entry_is_needed('dihedral')
		self.assertFalse(result1)

		self.env.globals = {'atom_style': 'hybrid'}
		self.env.globals.update({'atom_substyle1': 'bond'})
		self.env.globals.update({'atom_substyle2': 'angle'})
		result2 = writer.entry_is_needed('bond')
		self.assertTrue(result2)

		result3 = writer.entry_is_needed('dihedral')
		self.assertFalse(result3)



if __name__ == '__main__':
    ut.main(verbosity=2)
