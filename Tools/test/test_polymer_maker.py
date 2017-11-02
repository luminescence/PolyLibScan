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
        self.cfg_path = local_path.joinpath('data', 'updated_parameters8.yml')
        self.env = Environment(self.cfg_path.as_posix())

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
        for monomer in polymer.data['monomers']:
            self.assertIn('BP', monomer.name)

    # def test_autorepulsion(self):
    #     poly_creator = PolymerCreator(self.env, ['BP', 'C5-BP', 'CBS'], 
    #                             length=9,
    #                             mode='cycle')
    #     polymer = poly_creator.create()
    #     type_ids = polymer.env.atom_type.values()
    #     self.assertEqual(polymer.env.ff['affinity'][(type_ids[0], type_ids[1])].epsilon, 0)
    #     poly_creator.add_auto_repulsion(polymer)
    #     self.assertNotEqual(polymer.env.ff['affinity'][(type_ids[0], type_ids[1])].epsilon, 0)

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
