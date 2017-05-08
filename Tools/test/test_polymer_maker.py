from PolyLibScan.Tools.polymer2lmp import PolymerCreator
from PolyLibScan.Tools.environment import Environment
import pathlib2 as pl
import numpy as np

import mock
import unittest as ut

local_path = pl.Path(__file__).absolute().parent

class TestPolymer_Maker(ut.TestCase):

    def __init__(self, *args, **kwargs):
        super(TestPolymer_Maker, self).__init__(*args, **kwargs)
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
                                  'Ani': {
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
        self.cfg_data['globals'] = {'affinity_file': '../data/affinities_Originals2.h5',
                                    'atom_style': 'angle',
                                    'bond_style': 'harmonic',
                                    'angle_style': 'harmonic',
                                    'pair_style': 'hybrid',
                                    'box_margin': 30}
        Environment._get_config = mock.MagicMock(return_value=self.cfg_data)
        self.env = Environment('dummy.cfg')

    def isProtein(self, molecule):
        self.assertTrue(hasattr(molecule, 'data'))
        self.assertTrue(hasattr(molecule, 'mol_type'))
        self.assertTrue(hasattr(molecule, 'Id'))
        self.assertTrue(hasattr(molecule, 'box'))

    def test_random_length(self):
        # random length | random order
        p_length = np.random.randint(5,100)
        poly_creator = PolymerCreator(self.env, ['BA', 'Ani', 'CBS'], 
                                    length=p_length, mode='random')
        polymer = poly_creator.create()
        self.isProtein(polymer)
        self.assertEqual(p_length, len(polymer.data['particles']))
        
        # random length | repeating order
        p_length = np.random.randint(5,100)
        poly_creator.length = p_length
        polymer = poly_creator.create()
        self.isProtein(polymer)
        self.assertEqual(p_length, len(polymer.data['particles']))
        
    def test_fixed_length(self):
        # fixed length (each occurring once) | random mode
        poly_creator = PolymerCreator(self.env, ['BA', 'Ani', 'CBS'], 
                                mode='random')
        polymer = poly_creator.create()
        self.isProtein(polymer)
        self.assertEqual(3, len(polymer.data['particles']))
        
        # fixed length (each occurring once) | ordered composition
        poly_creator.mode = 'cycle'
        polymer = poly_creator.create()
        self.isProtein(polymer)
        self.assertEqual(3, len(polymer.data['particles']))

    def test_fixed_sequence(self):
        sequence = 11*['CBS', 'Ani', 'BA']
        poly_creator = PolymerCreator(self.env, sequence, 
                                mode='cycle')
        polymer = poly_creator.create()
        self.assertEqual(len(polymer.data['particles']), len(sequence))
        self.assertEqual(polymer.data['particles'][5].type_.name, 'BA')

    def test_weight(self):
        p_length = np.random.randint(5,100)
        polymer = PolymerCreator(self.env, ['BA', 'Ani', 'CBS'], 
                                    length=p_length, 
                                    weights=[4, 0, 0],
                                    mode='random').create()
        for particle in polymer.data['particles']:
            self.assertEqual('BA', particle.type_.name)

    def test_autorepulsion(self):
        poly_creator = PolymerCreator(self.env, ['BA', 'Ani', 'CBS'], 
                                length=9,
                                mode='cycle')
        polymer = poly_creator.create()
        type_ids = polymer.env.atom_type.values()
        self.assertEqual(polymer.env.ff['affinity'][(type_ids[0], type_ids[1])].epsilon, 0)
        poly_creator.add_auto_repulsion(polymer)
        self.assertNotEqual(polymer.env.ff['affinity'][(type_ids[0], type_ids[1])].epsilon, 0)


if __name__ == '__main__':
    ut.main(verbosity=2)
