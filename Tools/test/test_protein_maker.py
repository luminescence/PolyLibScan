from PolyLibScan.Tools.pdb2lmp import ProteinCreator
from PolyLibScan.Tools.environment import Environment
import pathlib2 as pl

import mock
import unittest as ut

local_path = pl.Path(__file__).absolute().parent

class TestProtein_Maker(ut.TestCase):

    def __init__(self, *args, **kwargs):
        super(TestProtein_Maker, self).__init__(*args, **kwargs)
        self.cfg_data = {}
        self.cfg_data['Bonds'] = {'peptide': {'kind': 'harmonic',
                                                'coef1': 120,
                                                'coef2': 4.0},
                                    'ghost': {'kind': 'harmonic',
                                                'coef1': 120,
                                                'coef2': 4.0}}
        self.cfg_data['Angles'] = {'peptide': {'kind': 'harmonic',
                                                'coef1': 120,
                                                'coef2': 4.0}}
        self.cfg_data['Dihedrals'] = {'peptide': {'kind': 'harmonic',
                                                    'coef1': 120,
                                                    'coef2': 4.0}}
        self.cfg_data['Atoms'] = {'CA': {
                                    'mass': 10.0,
                                    'radius': 2.0},
                                  'CA_ghost': {
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

    def test_Create_default_Protein(self):
        protein = ProteinCreator(self.env, local_path.joinpath('data/1LYA.pdb').as_posix()).create()
        self.isProtein(protein)
        self.assertTrue(len(protein.data['particles'])== 2*(338))

    def test_ghost_option_non_valid(self):
        p_creator = ProteinCreator(self.env, local_path.joinpath('data/1LYA.pdb').as_posix(), 
            with_ghosts='non-valid string')
        self.assertRaises(ValueError)

    def test_ghost_option(self):
        protein = ProteinCreator(self.env, local_path.joinpath('data/1LYA.pdb').as_posix(), 
                                with_ghosts=True).create()
        self.isProtein(protein)     
        self.assertTrue(len(protein.data['particles'])== 2*(338))
        protein = ProteinCreator(self.env, local_path.joinpath('data/1LYA.pdb').as_posix(), 
                                    with_ghosts=False).create()
        self.isProtein(protein)
        self.assertTrue(len(protein.data['particles'])== 338)

    def test_move_with_ghosts(self):
        pass

    def test_cg_lvl_non_valid(self):
        p_creator = ProteinCreator(self.env, local_path.joinpath('data/1LYA.pdb').as_posix(), 
                                  cg_lvl='non-valid string')
        self.assertRaises(ValueError)

    def test_change_to_resBased(self):
        p_creator = ProteinCreator(self.env, local_path.joinpath('data/1LYA.pdb').as_posix(), 
                                  cg_lvl='geometric_center')
        protein = p_creator.create()
        p_creator.change_to_res_based(protein, local_path.joinpath('data/amino_acids.yml').as_posix())

        self.isProtein(protein)
        self.assertEqual(protein.data['particles'][0].residue[0], protein.data['particles'][0].type_.name)
        self.assertIn('CA', protein.env.atom_type)
        self.assertIn('CA_ghost', protein.env.atom_type)
        self.assertIn('BP', protein.env.atom_type)
        self.assertEqual(len(protein.env.atom_type._defined_types), 23)
        # BP was never used
        self.assertEqual(len(protein.env.atom_type), 22)

    def test_cg_lvl_backbone(self):
        protein = ProteinCreator(self.env, local_path.joinpath('data/1LYA.pdb').as_posix(), 
                                  cg_lvl='backbone').create()
        self.isProtein(protein)

    def test_cg_lvl_geo_center(self):
        protein = ProteinCreator(self.env, local_path.joinpath('data/1LYA.pdb').as_posix(), 
                                  cg_lvl='geometric_center').create()
        self.isProtein(protein)
        first_particle = protein.data['particles'][0]
        self.assertEqual(first_particle.residue, ('GLY', 'A', (' ', 1, ' ')))
        self.assertAlmostEqual(first_particle.position[0],  9.93185711)
        self.assertAlmostEqual(first_particle.position[1],  81.30843353)
        self.assertAlmostEqual(first_particle.position[2],  34.65929031)
    '''    
    def test_cg_lvl_bbsc(self):
        
        protein = ProteinCreator(self.env, local_path.joinpath('data/1LYA.pdb').as_posix(), 
                                  cg_lvl='bb+sc').create()
        self.isProtein(protein)

    def test_cg_lvl_full(self):
        protein = ProteinCreator(self.env, local_path.joinpath('data/1LYA.pdb').as_posix(), 
                                  cg_lvl='full').create()
        self.isProtein(protein)
    '''

    def test_With_add_protein_properties(self):
        protein = ProteinCreator(self.env, local_path.joinpath('data/1LYA.pdb').as_posix(), 
                                  add_protein_binding=True, with_ghosts=False).create()
        self.isProtein(protein)
        self.assertTrue(len(protein.data['bonds'])== 338-1)
        self.assertTrue(len(protein.data['angles'])== 338-2)
        self.assertTrue(len(protein.data['dihedrals'])== 338-3)

    def test_Without_add_protein_properties(self):
        protein = ProteinCreator(self.env, local_path.joinpath('data/1LYA.pdb').as_posix(), 
                                  add_protein_binding=False, with_ghosts=False).create()
        self.isProtein(protein)
        self.assertTrue(len(protein.data['bonds'])== 0)
        self.assertTrue(len(protein.data['angles'])== 0)
        self.assertTrue(len(protein.data['dihedrals'])== 0)

if __name__ == '__main__':
    ut.main(verbosity=2)
