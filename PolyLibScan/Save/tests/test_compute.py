import unittest as ut
import PolyLibScan.Save.compute as compute
import mock
import numpy as np
import pathlib2 as pl
import os
local_path = pl.Path(os.path.dirname(os.path.realpath(__file__)))

class Test_compute(ut.TestCase):
    '''Test Cases for compute.
    '''

    def __init__(self, *args, **kwargs):
        super(Test_compute, self).__init__(*args, **kwargs)

    def test_dist(self):
        dtype = {'names':['atom_type', 'x', 'y', 'z'],
                 'formats': [np.int] + 3*[np.float]}
        vec1 = np.array([(1, 3.0, 5.0, 0.0)], dtype=dtype)
        vec2 = np.array([(2, 0.0, 5.0, 4.0)], dtype=dtype)
        results = compute._dist(vec1, vec2)

        self.assertAlmostEqual(results, 5.0)

    def test_distance_to_active_site(self):
        dtype = {'names':['atom_type', 'x', 'y', 'z'],
                 'formats': [np.int] + 3*[np.float]}
        coordinates = np.array([(1, 3.0, 4.0, 0.0),
                                (1, 3.0, 5.0, 10.0),
                                (1, 3.0, 5.0, 15.0),
                                (2, 0.0, -5.0, 0.0),
                                (2, 0.0, -10.0, 0.0),
                                (3, 0.0, 0.0, 0.0),], dtype=dtype)
        poly_ids = np.array([2,3])
        active_site = np.array([1,2,3])

        class dummy_db(object):

            def __init__(self):
                self.sequence = [2, 2, 3]
                self.traj_type_order = coordinates['atom_type']

        results = compute.distance_to_active_site(coordinates,
                                                  db=dummy_db(),
                                                  polymer_ids=poly_ids,
                                                  active_site_no=active_site)

        self.assertAlmostEqual(results, 5.0)


if __name__ == '__main__':
    ut.main(verbosity=2)
