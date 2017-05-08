import PolyLibScan.Tools.pdb2lmp as p2l
import PolyLibScan.Tools.lmp_types as lmp_t
import PolyLibScan.Tools.lmp_particlesAndInteractions as lmp_pi

import mock
import unittest

class TestSubClasses(unittest.TestCase):

    def test_atom_class(self):
        p = {'mass' : 1.0,
             'radius': 4.0,
             'interacting': 'True',
             'charge': -5}
        test_atom = lmp_t.AtomType('ca', p)
        self.assertEqual(test_atom.name, 'ca')
        self.assertEqual(test_atom.Id, None)
        self.assertEqual(test_atom.mass, 1.0)
        self.assertEqual(test_atom.charge, -5.0)
        self.assertEqual(test_atom.interacting, True)
        self.assertEqual(test_atom.radius, 4.0)

    def test_particle_class(self):
        molecule = mock.MagicMock()
        molecule.Id = 137
        molecule.name = 'test_molecule'
        p = {'mass' : 1.0,
             'radius': 4.0,
             'interacting': 'True',
             'charge': -5}
    	test_atom = lmp_t.AtomType('ca', p)
    	test_particle = lmp_pi.Particle(molecule, 1, test_atom, [1,2,3])
    	self.assertEqual(test_particle.Id, 1)
    	self.assertEqual(test_particle.position, [1,2,3])
        self.assertEqual(test_particle.mol_name, molecule.name)
        self.assertEqual(test_particle.mol_id, molecule.Id)

if __name__ == '__main__':
    unittest.main()
