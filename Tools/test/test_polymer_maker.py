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
                                'coef1': 120,
                                'coef2': 4.0,
                                'kind': 'harmonic'},
                                'CBS': {
                                'coef1': 310, 
                                'coef2': 8.896201747106357, 
                                'kind': 'harmonic'},
                                'BP': {
                                'coef1': 310, 
                                'coef2': 5.7202124353115185, 
                                'kind': 'harmonic'},
                                'C5-BP': {
                                'coef1': 310, 
                                'coef2': 10.987927308594289, 
                                'kind': 'harmonic'}
                                }
        self.cfg_data['Angles'] = {'polymer': {
                                'kind': 'harmonic',
                                'coef1': 120,
                                'coef2': 4.0}
                                }
        self.cfg_data['Atoms'] = {'CBS':{
                                'bb': {
                                'charge': 0,
                                'hydrophobicity': 1.2469999999999999,
                                'mass': 69.0819,
                                'position': [0.0, 0.0, 0.0],
                                'radius': 2.7354867877074724},
                                'sc': {
                                'charge': 1,
                                'hydrophobicity': 1.2469999999999999,
                                'mass': 238.2464599999999,
                                'position': [0.0, 8.896201747106357, 0.0],
                                'radius': 6.266729891400137}},
                                'BP': {
                                'bb': {
                                'charge': 0,
                                'hydrophobicity': 2.937800000000001,
                                'mass': 69.0819,
                                'position': [0.0, 0.0, 0.0],
                                'radius': 2.6353298825293203},
                                'sc': {
                                'charge': -2,
                                'hydrophobicity': 2.937800000000001,
                                'mass': 306.168782,
                                'position': [0.0, 5.7202124353115185, 0.0],
                                'radius': 8.026821000359863}},
                                'C5-BP': {
                                'bb': {'charge': 0,
                                'hydrophobicity': 3.6151999999999997,
                                'mass': 69.0819,
                                'position': [0.0, 0.0, 0.0],
                                'radius': 2.7651371166249987},
                                'sc': {
                                'charge': -2,
                                'hydrophobicity': 3.6151999999999997,
                                'mass': 405.2998419999999,
                                'position': [0.0, 10.987927308594289, 0.0],
                                'radius': 10.578341599506002}}
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
        self.cfg_data['Dihedrals'] = {'polymer': {
                                    'kind': 'harmonic',
                                    'coef1': 0.1,
                                    'coef2': 0.5,
                                    'coef3': 1.5,
                                    'coef4': 120}
                                    }
        self.cfg_data['Monomers'] = {'CBS': {
                                    'Angles': [],
                                    'Bonds': [['bb', 'sc']],
                                    'Dihedrals': [],
                                    'Part': ['bb', 'sc']},
                                    'BP': {
                                    'Angles': [],
                                    'Bonds': [['bb', 'sc']],
                                    'Dihedrals': [],
                                    'Part': ['bb', 'sc']},
                                    'C5-BP': {'Angles': [],
                                    'Bonds': [['bb', 'sc']],
                                    'Dihedrals': [],
                                    'Part': ['bb', 'sc']},
                                    }
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
        poly_creator = PolymerCreator(self.env, ['BP', 'C5-BP', 'CBS'], 
                                    length=p_length, mode='random')
        polymer = poly_creator.create()
        self.isProtein(polymer)
        self.assertEqual(p_length*2, len(polymer.data['particles']))
        
        # random length | repeating order
        p_length = np.random.randint(5,100)
        poly_creator.length = p_length
        polymer = poly_creator.create()
        self.isProtein(polymer)
        self.assertEqual(p_length*2, len(polymer.data['particles']))
        
    def test_fixed_length(self):
        # fixed length (each occurring once) | random mode
        poly_creator = PolymerCreator(self.env, ['BP', 'C5-BP', 'CBS'], 
                                mode='random')
        polymer = poly_creator.create()
        self.isProtein(polymer)
        self.assertEqual(6, len(polymer.data['particles']))
        
        # fixed length (each occurring once) | ordered composition
        poly_creator.mode = 'cycle'
        polymer = poly_creator.create()
        self.isProtein(polymer)
        self.assertEqual(6, len(polymer.data['particles']))

    def test_fixed_sequence(self):
        sequence = 11*['CBS', 'C5-BP', 'BP']
        poly_creator = PolymerCreator(self.env, sequence, 
                                mode='cycle')
        polymer = poly_creator.create()
        self.assertEqual(len(polymer.data['particles']), len(sequence)*2)
        self.assertIn('BP', polymer.data['particles'][5].type_.name)

    def test_weight(self):
        p_length = np.random.randint(5,100)
        polymer = PolymerCreator(self.env, ['BP', 'C5-BP', 'CBS'], 
                                    length=p_length, 
                                    weights=[4, 0, 0],
                                    mode='random').create()
        for particle in polymer.data['particles']:
            self.assertIn('BP', particle.type_.name)

    def test_autorepulsion(self):
        poly_creator = PolymerCreator(self.env, ['BP', 'C5-BP', 'CBS'], 
                                length=9,
                                mode='cycle')
        polymer = poly_creator.create()
        type_ids = polymer.env.atom_type.values()
        self.assertEqual(polymer.env.ff['affinity'][(type_ids[0], type_ids[1])].epsilon, 0)
        poly_creator.add_auto_repulsion(polymer)
        self.assertNotEqual(polymer.env.ff['affinity'][(type_ids[0], type_ids[1])].epsilon, 0)

    def test_create_polymer_bonds(self):
        p_length = np.random.randint(5, 100)
        poly_creator = PolymerCreator(self.env, ['BP', 'C5-BP', 'CBS'],
                                length=p_length,
                                mode='random')
        polymer = poly_creator.create()
        self.isProtein(polymer)
        self.assertEqual(len(polymer.data['bonds']), (p_length*2)-1)

    def test_create_polymer_angles(self):
        p_length = np.random.randint(5, 100)
        poly_creator = PolymerCreator(self.env, ['BP', 'C5-BP', 'CBS'],
                                length=p_length,
                                mode='random')
        polymer = poly_creator.create()
        self.isProtein(polymer)
        self.assertEqual(len(polymer.data['angles']), ((p_length-2)*3)+2)

    def test_create_polymer_dihedrals(self):
        p_length = np.random.randint(5, 100)
        poly_creator = PolymerCreator(self.env, ['BP', 'C5-BP', 'CBS'],
                                length=p_length,
                                mode='random')
        polymer = poly_creator.create()
        self.isProtein(polymer)
        self.assertEqual(len(polymer.data['dihedrals']), ((p_length-3)+(p_length-1)))

if __name__ == '__main__':
    ut.main(verbosity=2)
