import PolyLibScan.Analysis as lmp_lys
import pathlib2 as pl
import mock
import numpy as np

import unittest as ut

local_path = pl.Path(__file__).absolute().parent

class TestPymolVis(ut.TestCase):

    def __init__(self, *args, **kwargs):
        super(TestPymolVis, self).__init__(*args, **kwargs)
        box = np.array([[-10, -10, -10], [10, 10, 10]])
        self.sim = mock.MagicMock()
        self.sim.db_path = local_path.joinpath('data')
        self.sim._calc_protein_box = mock.Mock(return_value=box)
        self.sim._particle_ids = [1,2,3]

    def test_monomer_id_prop(self):
        density = lmp_lys.DensityContainer(self.sim, 1)
        self.assertEqual(density.monomer_id[0], 1)

        density2 = lmp_lys.DensityContainer(self.sim, [1,2])
        self.assertEqual(density.monomer_id[0], 1)
        self.assertEqual(density.monomer_id[1], 2)

        density3 = lmp_lys.DensityContainer(self.sim, np.array([1]))
        self.assertEqual(density.monomer_id[0], 1)

        with self.assertException as e:
            density4 = lmp_lys.DensityContainer(self.sim, np.array([5]))

        with self.assertException as e:
            density5 = lmp_lys.DensityContainer(self.sim, '5')

        density6 = lmp_lys.DensityContainer(self.sim, 'all')
        self.assertEqual(density.monomer_id, [1,2,3])

    def test_create_empty_map(self):
        pass

    def test_map_path(self):
        pass

    def test_create_atom_type_filter(self):
        pass

    def test_add_run_to_epitopsy_map(self):
        pass

    def test_create_epitopsy_map(self):
        pass

    def test_save(self):
        pass

if __name__ == '__main__':
    ut.main(verbosity=2)
