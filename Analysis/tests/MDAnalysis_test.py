import PolyLibScan.Analysis as lmp_lys
from PolyLibScan.Analysis.MDAnalysis_interface import MdaRun
import pathlib2 as pl

import unittest as ut

local_path = pl.Path(__file__).absolute().parent


class TestMdaInterface(ut.TestCase):

    def __init__(self, *args, **kwargs):
        super(TestMdaInterface, self).__init__(*args, **kwargs)
        self.project = lmp_lys.Project(local_path.joinpath('data_TRP-TRP_BP_monomer/').as_posix(),
                                       experimental_data=None)
        self.job = self.project.jobs[0]
        first_run = self.job[0]
        self.mda_run = MdaRun(first_run)
        self.mda_job = MdaJob(self.job)
        self.mda_project = MdaProject(self.project)
        self.trypsin_bead_length = 224

    def test_end_distance_BP_trypsin_A(self):
        # test_margin
        trypsin_A_selection = 'bynum 1-%s' % self.trypsin_bead_length
        BP_monomer_selection = 'bynum %s' % (2 * self.trypsin_bead_length + 1)
        distance_BP_trypsin_A = self.mda_run.comp_min_distance_between_selections(trypsin_A_selection,
                                                                                  BP_monomer_selection)
        self.assertAlmostEqual(distance_BP_trypsin_A[2], 21.3, places=1)    # comparison: pymol
        self.assertEqual(len(distance_BP_trypsin_A), 3)  # three snapshots: 0, 2000 and 4000 time steps

    def test_project_job_run(self):
        sel1 = 'bynum 1'
        sel2 = 'bynum 2'

        run_distance = self.mda_run.comp_min_distance_between_selections(sel1, sel2)
        job_distance = self.mda_job.comp_min_distance_between_selections(sel1, sel2)
        project_distance = self.mda_project.comp_min_distance_between_selections(sel1, sel2)

        self.assertTrue((run_distance == job_distance[0]).all())
        self.assertTrue((job_distance == project_distance[0].to_pandas().transpose()).all()[0])


if __name__ == '__main__':
    ut.main(verbosity=2)
