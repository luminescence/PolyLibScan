import PolyLibScan.Tools.polymer2lmp as poly
import PolyLibScan.Tools.lmp_helpers as helpers
from PolyLibScan.Tools.environment import Environment
import numpy as np
import pathlib2 as pl

import mock
import unittest as ut

local_path = pl.Path(__file__).absolute().parent

class TestLmpObject(ut.TestCase):

    def __init__(self, *args, **kwargs):
        super(TestLmpObject, self).__init__(*args, **kwargs)
        self.cfg_data = {}
        self.cfg_data['Bonds'] = {'polymer': {
                                'kind': 'harmonic',
                                'coef1': 120,
                                'coef2': 4.0}
                                }
        self.cfg_data['Angles'] = {'polymer': {
                                'kind': 'harmonic',
                                'coef1': 120,
                                'coef2': 4.0}
                                }
        self.cfg_data['Atoms'] = {'CBS': {
                                    'mass': 10.0,
                                    'radius': 2.0},
                                  'BA': {
                                    'mass': 10.0,
                                    'radius': 2.0},
                                  'BP': {
                                    'mass': 10.0,
                                    'radius': 2.0}
                                    }
        self.cfg_data['Pairs'] = {'lj': {
                                    'kind': 'morse',
                                    'repulsive_only': 0,
                                    'single_type_parametrised': True,
                                    'cutoff': 45,
                                    'coef1': 0.0,
                                    'coef2': 0.4},
                                  'affinity': {
                                    'kind': 'morse',
                                    'repulsive_only': 0,
                                    'single_type_parametrised': True,
                                    'cutoff': 45,
                                    'coef1': 0.0,
                                    'coef2': 0.4}
                                    }
        self.cfg_data['Dihedrals'] = {}
        self.cfg_data['globals'] = {'affinity_file': '../data/affinties.h5',
                                    'atom_style': 'angle',
                                    'bond_style': 'harmonic',
                                    'angle_style': 'harmonic',
                                    'pair_style': 'hybrid',
                                    'box_margin': 30}
        Environment._get_config = mock.MagicMock(return_value=self.cfg_data)
        self.env = Environment('dummy.cfg')

    def test_move(self):
        polymer = poly.PolymerCreator(self.env, ['BA', 'BP', 'CBS'], 
                                        length=20, mode='random').create()
        shift = np.random.uniform(-5, 5, 3)
        pos1 = polymer.box
        polymer.move(shift)
        pos2 = polymer.box
        difference = pos2-pos1
        self.assertAlmostEqual(shift[0], difference[0,0])
        self.assertAlmostEqual(shift[1], difference[1,0])
        self.assertAlmostEqual(shift[2], difference[2,0])
        self.assertAlmostEqual(shift[0], difference[0,1])
        self.assertAlmostEqual(shift[1], difference[1,1])
        self.assertAlmostEqual(shift[2], difference[2,1])

    def test_bonds(self):
        polymer = poly.PolymerCreator(self.env, ['BA', 'BP', 'CBS'], 
                                      length=20, mode='random').create()
        for bond in polymer.data['bonds']:
            self.assertEqual(len(bond.members), 2)
            self.assertEqual(bond.type_.name, 'polymer')
        self.assertEqual(len(polymer.data['bonds']), len(polymer.data['particles'])-1)        

    def test_rotation(self):
        polymer = poly.PolymerCreator(self.env, ['BA', 'BP', 'CBS'], 
                                      length=20, mode='random').create()
        rot_axis = np.concatenate((np.array([1.0]),
                                   np.random.uniform(0.0, 2 * np.pi, 1),
                                   np.random.uniform(0.0,     np.pi, 1)) )
        axis = helpers.Cartesian_np(rot_axis)[0]
        angle = np.random.uniform(0.0, 2 * np.pi)

        init_length = np.linalg.norm((polymer.box[:,1]-polymer.box[:,0]))

        polymer.rotate(axis, angle)

        post_length = np.linalg.norm((polymer.box[:,1]-polymer.box[:,0]))
        self.assertAlmostEqual(init_length, post_length)

    def test_length(self):
        # random length | random order
        p_length = np.random.randint(5,100)
        polymer = poly.PolymerCreator(self.env, ['BA', 'BP', 'CBS'], 
                                      length=p_length, mode='random').create()
        self.assertEqual(p_length, len(polymer.data['particles']))
        # random length | repeating order
        p_length = np.random.randint(5,100)
        polymer = poly.PolymerCreator(self.env, ['BA', 'BP', 'CBS'], 
                                      length=p_length).create()
        self.assertEqual(p_length, len(polymer.data['particles']))
        # fixed length (each occurring once) | random mode
        polymer = poly.PolymerCreator(self.env, ['BA', 'BP', 'CBS'], 
                                      mode='random').create()
        self.assertEqual(3, len(polymer.data['particles']))
        # fixed length (each occurring once) | ordered composition
        polymer = poly.PolymerCreator(self.env, ['BA', 'BP', 'CBS'], 
                                      mode='random').create()
        self.assertEqual(3, len(polymer.data['particles']))

if __name__ == '__main__':
    ut.main(verbosity=2)
