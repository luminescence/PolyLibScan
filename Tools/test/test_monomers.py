from PolyLibScan.Tools import PolymerCreator
from PolyLibScan.Tools.monomers import Monomer
from PolyLibScan.Tools.environment import Environment
from operator import add

import pathlib2 as pl
import unittest as ut

local_path = pl.Path(__file__).absolute().parent


class BaseTestClass:
    """The base class is wrapped in another class so it itself won't be executed by the testing facility."""

    class TestMonomers(ut.TestCase):

        def generate_N_monomers(self):
            monomers = []
            monomer_position = [0, 0, 0]
            self.N_monomers = 4  # 4 allows to test dihedrals in the single bead model
            for x in range(self.N_monomers):
                monomers.append(Monomer(monomer_position, 'CBS', self.env.monomer_type['CBS'], self.env, []))
                monomer_position = map(add, monomer_position, [4, 0, 0])

            for idx, x in enumerate(monomers):
                if idx < self.N_monomers-1:
                    monomers[idx].bind_with(monomers[idx+1])
                    if idx > 0:
                        monomers[idx].bb_angle_with(monomers[idx-1], monomers[idx+1])
                        if idx < self.N_monomers - 2:
                            monomers[idx].bb_dihedral_with(monomers[idx-1], monomers[idx+1], monomers[idx+2])

            return monomers

        @staticmethod
        def number_of_unique_DOF(monomers, dof):
            return len(PolymerCreator.unique_ordered_DOF(monomers, dof))

        def test_Monomer(self):
            mono1, mono2, mono3, mono4 = self.generate_N_monomers()

            self.assertEqual('CBS', mono1.name)

        def test_bind_to(self):
            monomers = self.generate_N_monomers()

            self.assertEqual(len(monomers[0].bonds) + 1, len(monomers[1].bonds))
            self.assertEqual(len(monomers[0].bonds), len(monomers[-1].bonds))
            self.assertEqual(len(monomers[1].bonds), len(monomers[-2].bonds))
            self.assertIs(monomers[0].bonds[-1], monomers[1].bonds[-2])
            self.assertIs(monomers[1].bonds[-1], monomers[2].bonds[-2])
            self.assertIs(monomers[2].bonds[-1], monomers[3].bonds[-1])

        def test_number_of_bonds(self):
            monomers = self.generate_N_monomers()
            unique_bonds = self.number_of_unique_DOF(monomers, 'bonds')

            self.assertEqual(unique_bonds, self.beads_per_monomer * self.N_monomers - 1)

        def test_bb_angle_with(self):
            monomers = self.generate_N_monomers()

            unique_bb_angles = self.number_of_unique_DOF(monomers, 'angles')
            self.assertEqual(unique_bb_angles, self.N_monomers-2)

            self.assertEqual(len(monomers[0].angles), 1)
            self.assertEqual(len(monomers[0].angles) + 1, len(monomers[1].angles))
            self.assertEqual(len(monomers[0].angles), len(monomers[-1].angles))
            self.assertEqual(len(monomers[1].angles), len(monomers[-2].angles))
            self.assertIs(monomers[0].angles[0], monomers[1].angles[0])

        def test_bb_dihedral_with(self):
            monomers = self.generate_N_monomers()

            unique_bb_angles = self.number_of_unique_DOF(monomers, 'dihedrals')
            self.assertEqual(unique_bb_angles, self.N_monomers-3)

            self.assertEqual(len(monomers[0].dihedrals), 1)
            self.assertEqual(len(monomers[0].dihedrals), len(monomers[-1].dihedrals))
            self.assertEqual(len(monomers[1].dihedrals), len(monomers[-2].dihedrals))


class TestMonomers1bead(BaseTestClass.TestMonomers):

    def __init__(self, *args, **kwargs):
        super(BaseTestClass.TestMonomers, self).__init__(*args, **kwargs)
        self.env = Environment(local_path.joinpath('data',
                                                   'parameters_hp2.yml').as_posix())
        self.beads_per_monomer = 1

    def test_create_particles(self):
        mono1, mono2, mono3, mono4 = self.generate_N_monomers()

        self.assertIn('CBS_bb', mono1.particles)
        self.assertEqual(mono1.name, mono2.name)


class TestMonomers2bead(BaseTestClass.TestMonomers):

    def __init__(self, *args, **kwargs):
        super(BaseTestClass.TestMonomers, self).__init__(*args, **kwargs)
        self.env = Environment(local_path.joinpath('data',
                                                   'updated_parameters8.yml').as_posix())
        self.beads_per_monomer = 2

    def test_create_particles(self):
        mono1, mono2, mono3, mono4 = self.generate_N_monomers()

        self.assertIn('mono_bb', mono1.particles)
        self.assertIn('CBS_sc', mono1.particles)

    def test_angle_with(self):
        monomers = self.generate_N_monomers()

        for idx, x in enumerate(monomers):
            if idx < self.N_monomers - 1:
                monomers[idx].angle_with(monomers[idx + 1])
            else:
                monomers[idx].angle_with(monomers[idx - 1])

        self.assertEqual(len(monomers[0].angles), 2)
        self.assertEqual(len(monomers[-1].angles), 3)
        self.assertEqual(len(monomers[1].angles), 4)
        self.assertEqual(len(monomers[-2].angles), 5)

        unique_angles = self.number_of_unique_DOF(monomers, 'angles')
        self.assertEqual(unique_angles, self.beads_per_monomer*self.N_monomers - 2)

    def test_dihedral_with(self):
        """Currently, there is one dihedral missing. One would need to implement the edge-case
        of the 0th monomer binding sc-bb-bb-bb."""
        monomers = self.generate_N_monomers()

        for idx, x in enumerate(monomers):
            if idx < self.N_monomers - 1:
                monomers[idx].dihedral_with(monomers[idx + 1])

        self.assertEqual(len(monomers[0].dihedrals), 2)
        self.assertEqual(len(monomers[1].dihedrals), 3)

        unique_angles = self.number_of_unique_DOF(monomers, 'dihedrals')
        self.assertEqual(unique_angles, self.beads_per_monomer*self.N_monomers - 4)


if __name__ == '__main__':
    ut.main(verbosity=2)
