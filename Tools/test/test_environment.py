from PolyLibScan.Tools.environment import Environment
import pathlib2 as pl
import mock
import numpy as np

import unittest as ut

local_path = pl.Path(__file__).absolute().parent


class TestEnvironment(ut.TestCase):

    def __init__(self, *args, **kwargs):
        super(TestEnvironment, self).__init__(*args, **kwargs)
        self.cfg_data = {}
        self.cfg_data['Bonds'] = {'bond1': {'kind': 'harmonic',
                                            'coef1': 120,
                                            'coef2': 4.0}
                                            }
        self.cfg_data['Atoms'] = {'CA': {'mass': 10.0,
                                         'radius': 2.0}}
        self.cfg_data['Pairs'] = {  'lj': {'kind': 'morse',
                                            'repulsive_only': 0,
                                            'single_type_parametrised': True,
                                            'cutoff': 45,
                                            'coef1': 0.0,
                                            'coef2': 0.4},
                                    'affinity': {'kind': 'morse',
                                                'repulsive_only': 0,
                                                'single_type_parametrised': True,
                                                'cutoff': 45,
                                                'coef1': 0.0,
                                                'coef2': 0.4}
                                    }
        self.cfg_data['Angles'] = {}
        self.cfg_data['Dihedrals'] = {}
        self.cfg_data['globals'] = {'affinity_file': '../data/affinities_Originals2.h5',
                                    'atom_style': 'angle',
                                    'bond_style': 'harmonic',
                                    'angle_style': 'harmonic',
                                    'pair_style': 'hybrid',
                                    'box_margin': 30,
                                    'surface_file': False}

        Environment._get_config = mock.MagicMock(return_value=self.cfg_data)
        self.env = Environment('test.cfg')

    def test_load_atom_types(self):

        self.assertEqual(self.env.atom_type['CA'].mass, 10.0)
        self.assertEqual(self.env.atom_type['CA'].radius, 2.0)
        self.assertEqual(self.env.atom_type['CA'].Id, 1)

    def test_load_bond_types(self):

        self.assertEqual(self.env.bond_type['bond1'].name, 'bond1')
        self.assertEqual(self.env.bond_type['bond1'].Id, 1)
        self.assertEqual(self.env.bond_type['bond1'].kind, 'harmonic')
        self.assertEqual(self.env.bond_type['bond1'].parameters['coeffs'][0], 120)
        self.assertEqual(self.env.bond_type['bond1'].parameters['coeffs'][1], 4.0)

    def test_hooks(self):
        self.assertTrue(len(self.env.atom_type._registered)==2)

    def test_forcefield_hooks(self):
        env = Environment('test.cfg')
        # adding atom type
        new_atom_type = mock.MagicMock(Id=1, name='test_atom1')
        env.atom_type['test_type'] = new_atom_type
        # the container will change the id..
        self.assertTrue((new_atom_type, env.atom_type['CA']) in env.ff['affinity'])

    def test_calc_box_empty(self):   
        env = Environment('test.cfg')
        self.assertEqual(env.box[0,0], 0.0)
        env.calc_box(add_margin=True)
        self.assertEqual(env.box[0,0], -30.0)

    def test_calc_box(self):
        first_box = self.env.box.copy()
        self.env.calc_box(add_margin=True)
        second_box = self.env.box.copy()

        self.assertNotEqual(first_box[0,0], second_box[0,0])


if __name__ == '__main__':
    ut.main(verbosity=2)
