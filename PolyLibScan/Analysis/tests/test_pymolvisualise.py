import PolyLibScan.Analysis.pymol_visualisation as pymol
import pathlib2 as pl
import mock
import numpy as np

import unittest as ut

local_path = pl.Path(__file__).absolute().parent

class TestPymolVis(ut.TestCase):

    def __init__(self, *args, **kwargs):
        super(TestPymolVis, self).__init__(*args, **kwargs)
        job = mock.Mock()
        self.pdb_file = local_path.joinpath('data', 'static', '4cha.pdb')
        self.py_vis = pymol.PymolVisualisation(self.pdb_file)

    def test_protein_path(self):
        self.assertEqual(str(self.py_vis.protein_path), str(self.pdb_file))

    def test_init_protein_path(self):
        results1 = self.py_vis._init_protein_path(self.pdb_file)
        self.assertEqual(results1, self.pdb_file)

        results2 = self.py_vis._init_protein_path(None, search_path=self.pdb_file.parent)
        self.assertEqual(results2, self.pdb_file)


    
    # def test_poly_poses(self):
    #     self.py_vis._poly_poses()

if __name__ == '__main__':
    ut.main(verbosity=2)
