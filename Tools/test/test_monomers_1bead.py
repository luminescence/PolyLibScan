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
            'updated_parameters8.yml').as_posix())

    def test_create_particles(self):
        mono1 = Monomer([0, 0, 0], 'CBS', self.env.monomer_type['CBS'], 
             self.env, [])
        mono2 = Monomer([4, 0, 0], 'CBS', self.env.monomer_type['CBS'], 
            self.env, [])

        self.assertIn('mono_bb', mono1.particles)
        self.assertNotIn('CBS_bb', mono1.particles) 
        self.assertEqual (mono1.name, mono2.name)

    def test_Monomer(self):
        mono1 = Monomer([0, 0, 0], 'CBS', self.env.monomer_type['CBS'], 
            self.env, [])

        self.assertEqual('CBS', mono1.name)

if __name__ == '__main__':
    ut.main(verbosity=2)
