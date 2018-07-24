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
        self.cfg_path = local_path.joinpath('data', 'parameters.yml')
        self.env = Environment(self.cfg_path.as_posix())


    def test_move(self):
        polymer = poly.PolymerCreator(self.env, ['Glu', 'BP', 'CBS'], 
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
        polymer = poly.PolymerCreator(self.env, ['Glu', 'BP', 'CBS'], 
                                      length=20, mode='random').create()
        for bond in polymer.data['bonds']:
            self.assertEqual(len(bond.members), 2)
            self.assertIn(bond.type_.name, ['polymer', 'Glu', 'BP', 'CBS'])
        self.assertEqual(len(polymer.data['bonds']), len(polymer.data['particles'])-1)        

    def test_rotation(self):
        polymer = poly.PolymerCreator(self.env, ['Glu', 'BP', 'CBS'], 
                                      length=20, mode='random').create()
        rot_axis = np.concatenate((np.array([1.0]),
                                   np.random.uniform(0.0, 2 * np.pi, 1),
                                   np.random.uniform(0.0,     np.pi, 1)) )
        axis = helpers.Cartesian_np(rot_axis)[0]
        angle = np.random.uniform(0.0, 2 * np.pi)

        init_length = np.linalg.norm((polymer.data['particles'][0].position - 
                                      polymer.data['particles'][1].position))

        polymer.rotate(axis, angle)

        post_length = np.linalg.norm((polymer.data['particles'][0].position - 
                                      polymer.data['particles'][1].position))
        self.assertAlmostEqual(init_length, post_length)

    def test_length(self):
        # random length | random order
        p_length = np.random.randint(5,100)
        polymer = poly.PolymerCreator(self.env, ['Glu', 'BP', 'CBS'], 
                                      length=p_length, mode='random').create()
        self.assertEqual(p_length, len(polymer.data['particles']))
        # random length | repeating order
        p_length = np.random.randint(5,100)
        polymer = poly.PolymerCreator(self.env, ['Glu', 'BP', 'CBS'], 
                                      length=p_length).create()
        self.assertEqual(p_length, len(polymer.data['particles']))
        # fixed length (each occurring once) | random mode
        polymer = poly.PolymerCreator(self.env, ['Glu', 'BP', 'CBS'], 
                                      mode='random').create()
        self.assertEqual(3, len(polymer.data['particles']))
        # fixed length (each occurring once) | ordered composition
        polymer = poly.PolymerCreator(self.env, ['Glu', 'BP', 'CBS'], 
                                      mode='random').create()
        self.assertEqual(3, len(polymer.data['particles']))

if __name__ == '__main__':
    ut.main(verbosity=2)
