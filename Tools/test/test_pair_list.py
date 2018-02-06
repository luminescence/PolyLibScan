from PolyLibScan.Tools.lmp_control import PairList
import pathlib2 as pl
import unittest as ut
import yaml

local_path = pl.Path(__file__).absolute().parent

class TestPairList(ut.TestCase):

    def __init__(self, *args, **kwargs):
        super(TestPairList, self).__init__(*args, **kwargs)
        self.global_par = {'affinity_file': '../data/affinties.h5',
                      'angle_style': 'harmonic',
                      'atom_style': 'hybrid',
                      'atom_substyle1': 'angle',
                      'atom_substyle2': 'charge',
                      'bond_style': 'harmonic',
                      'box_margin': 50,
                      'pair_style': 'hybrid/overlay',
                      'pair_substyle1': 'coulomb',
                      'pair_substyle2': 'lj',
                      'pair_substyle3': 'hydrophobic'}
        self.Pairs = {'coulomb': {'cutoff': 25,
                 'kappa': 0.55,
                 'kind': 'coul/debye',
                 'repulsive_only': 0,
                 'single_type_parametrised': True},
                 'hydrophobic': {'coef1': -0.5,
                 'cutoff': 6.0,
                 'kind': 'soft',
                 'repulsive_only': 0,
                 'single_type_parametrised': False},
                 'lj': {'coef1': 1.0,
                 'coef2': 4.0,
                 'cutoff': 9.0,
                 'kind': 'lj/cut',
                 'repulsive_only': 1,
                 'single_type_parametrised': False}}

        

    def test_with_substyles(self):
        List = PairList(self.global_par, self.Pairs)

        self.assertTrue(List.has_substyles)
        self.assertEqual(List.styles['main'], 'hybrid/overlay')
        self.assertEqual(List.styles['sub'], ['coulomb', 'lj', 'hydrophobic'])

    def test_without_substyles(self):
        local_globals = self.global_par.copy()
        local_globals['pair_style'] = 'coulomb'
        List = PairList(local_globals, self.Pairs)

        self.assertFalse(List.has_substyles)
        self.assertEqual(List.styles['main'], 'coulomb')
        self.assertEqual(List.styles['sub'], [])

    def test_pair_style_str(self):
        List = PairList(self.global_par, self.Pairs)
        
        self.assertEqual(List.pair_style_str, 'hybrid/overlay coul/debye 0.55 25 lj/cut 9.0 soft 6.0')

        local_globals = self.global_par.copy()
        local_globals['pair_style'] = 'coulomb'
        List2 = PairList(local_globals, self.Pairs)
        
        self.assertEqual(List2.pair_style_str, 'coul/debye 0.55 25')

    def test_single_parameterized(self):
        List = PairList(self.global_par, self.Pairs)

        self.assertEqual(List.single_parameterized(), ['coulomb'])

    def test_pair_coef_str(self):
        List = PairList(self.global_par, self.Pairs)
        self.assertEqual(List.pair_coef_str('coulomb'), 'pair_coeff * * coul/debye  25')
        self.assertEqual(List.pair_coef_str('coulomb', atom_type1=3), 'pair_coeff 3 * coul/debye  25')

if __name__ == '__main__':
    ut.main(verbosity=2)
