from PolyLibScan.Tools.environment import Environment
import pathlib2 as pl
import mock
import numpy as np

import unittest as ut

local_path = pl.Path(__file__).absolute().parent


class TestEnvironment(ut.TestCase):

    def __init__(self, *args, **kwargs):
        super(TestEnvironment, self).__init__(*args, **kwargs)
        self.cfg_path = local_path.joinpath('data', 'updated_parameters8.yml')
        self.env = Environment(self.cfg_path.as_posix())

    def test_load_atom_types(self):

        self.assertEqual(self.env.atom_type['BB_bb'].mass, 100.0)
        self.assertEqual(self.env.atom_type['BB_bb'].radius, 2.0)
        self.assertEqual(self.env.atom_type['BB_bb'].Id, 1)

    def test_load_bond_types(self):
        self.assertEqual(self.env.bond_type['bond1'].name, 'bond1')
        self.assertEqual(self.env.bond_type['bond1'].Id, 1)
        self.assertEqual(self.env.bond_type['bond1'].kind, 'harmonic')
        self.assertEqual(self.env.bond_type['bond1'].parameters['coeffs'][0], 120)
        self.assertEqual(self.env.bond_type['bond1'].parameters['coeffs'][1], 4.0)

    # def test_forcefield_hooks(self):
    #     env = Environment(self.cfg_path.as_posix())
    #     # adding atom type
    #     new_atom_type = mock.MagicMock(Id=1, name='test_atom1')
    #     env.atom_type['test_type'] = new_atom_type
    #     # the container will change the id..
    #     self.assertTrue((new_atom_type, env.atom_type['CA']) in env.ff['affinity'])

    def test_calc_box_empty(self):   
        env = Environment(self.cfg_path.as_posix())
        self.assertEqual(env.box[0,0], 0.0)
        env.calc_box(add_margin=True)
        self.assertEqual(env.box[0,0], -50.0)

    def test_calc_box(self):
        first_box = self.env.box.copy()
        self.env.calc_box(add_margin=True)
        second_box = self.env.box.copy()

        self.assertNotEqual(first_box[0,0], second_box[0,0])


if __name__ == '__main__':
    ut.main(verbosity=2)
