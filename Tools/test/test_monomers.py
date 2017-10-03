from PolyLibScan.Tools.monomers import Monomer
from PolyLibScan.Tools.environment import Environment

import pathlib2 as pl
import numpy as np

import mock
import unittest as ut

local_path = pl.Path(__file__).absolute().parent


class Test_Monomers(ut.TestCase):

    def __init__(self, *args, **kwargs):
        super(Test_Monomers, self).__init__(*args, **kwargs)
        self.env = Environment(local_path.joinpath('data', 
            'updated_parameters3.yaml').as_posix())

    def test_bind_to(self):
        mono1 = Monomer([0, 0, 0], 'CBS', self.env.monomer_type['CBS'], 
            self.env, [])
        mono2 = Monomer([4, 0, 0], 'CBS', self.env.monomer_type['CBS'], 
            self.env, [])

        mono1.bind_with(mono2)

        self.assertEqual(len(mono1.bonds) , 2)
        self.assertEqual(len(mono2.bonds), 2)
        self.assertIs(mono1.bonds[1].Id, mono2.bonds[1].Id)
        self.assertIs(mono1.bonds[1], mono2.bonds[1])

    def test_angle_with(self):
        mono1 = Monomer([0, 0, 0], 'CBS', self.env.monomer_type['CBS'], 
            self.env, [])
        mono2 = Monomer([4, 0, 0], 'CBS', self.env.monomer_type['CBS'], 
            self.env, [])

        mono1.angle_with(mono2)

        self.assertEqual(len(mono1.angles) , 1)
        self.assertEqual(len(mono2.angles), 1)
        self.assertIs(mono1.angles[0].Id, mono2.angles[0].Id)
        self.assertIs(mono1.angles[0], mono2.angles[0])
    
    def test_dihedral_with(self):
        mono1 = Monomer([0, 0, 0], 'CBS', self.env.monomer_type['CBS'], 
            self.env, [])
        mono2 = Monomer([4, 0, 0], 'CBS', self.env.monomer_type['CBS'], 
            self.env, [])

        mono1.dihedral_with(mono2)

        self.assertEqual(len(mono1.dihedrals) , 1)
        self.assertEqual(len(mono2.dihedrals), 1)
        self.assertIs(mono1.dihedrals[0].Id, mono2.dihedrals[0].Id)
        self.assertIs(mono1.dihedrals[0], mono2.dihedrals[0])

    def test_create_particles(self):
        mono1 = Monomer([0, 0, 0], 'CBS', self.env.monomer_type['CBS'], 
             self.env, [])
        mono2 = Monomer([4, 0, 0], 'CBS', self.env.monomer_type['CBS'], 
            self.env, [])

        self.assertIn('bb', mono1.particles)

    def test_Monomer(self):
        mono1 = Monomer([0, 0, 0], 'CBS', self.env.monomer_type['CBS'], 
            self.env, [])

        self.assertEqual('CBS', mono1.name)

if __name__ == '__main__':
    ut.main(verbosity=2)
