import pathlib2 as pl
import unittest as ut

import xarray as xr

import PolyLibScan.Analysis as lmp_lys
from PolyLibScan.Analysis.MDAnalysis_interface import MdaRun, MdaJob, MdaProject

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

        self.sel1 = 'bynum 1'
        self.sel2 = 'bynum 2'

    def test_end_distance_BP_trypsin_A(self):
        # test_margin
        trypsin_A_selection = 'bynum 1-%s' % self.trypsin_bead_length
        BP_monomer_selection = 'bynum %s' % (2 * self.trypsin_bead_length + 1)
        distance_BP_trypsin_A = self.mda_run.comp_min_distance_between_selections(trypsin_A_selection,
                                                                                  BP_monomer_selection)
        self.assertAlmostEqual(distance_BP_trypsin_A[2], 21.3, places=1)    # comparison: pymol

    def test_range(self):
        distance_BP_trypsin_A = self.mda_run.comp_min_distance_between_selections(self.sel1,
                                                                                  self.sel2)
        self.assertEqual(len(distance_BP_trypsin_A), 3)  # three snapshots: 0, 2000 and 4000 time steps

    def test_indexation(self):
        """check that the index of the pd.Series is correct"""
        distance_BP_trypsin_A_end = self.mda_run.comp_min_distance_between_selections(self.sel1,
                                                                                      self.sel2,
                                                                                      snapshots=-1)

        self.assertTrue((distance_BP_trypsin_A_end.index == 2).all())

    def test_project_job_run(self):
        run_distance = self.mda_run.comp_min_distance_between_selections(self.sel1, self.sel2)
        job_distance = self.mda_job.comp_min_distance_between_selections(self.sel1, self.sel2)
        project_distance = self.mda_project.comp_min_distance_between_selections(self.sel1, self.sel2)

        self.assertTrue((run_distance == job_distance[0]).all())
        self.assertTrue((job_distance == project_distance[0].to_pandas().transpose()).all()[0])

    def test_timestep_selection(self):
        for hierarchy_lvl, mda_instance in enumerate([self.mda_run, self.mda_job, self.mda_project]):
            run_distance = mda_instance.comp_min_distance_between_selections(self.sel1, self.sel2)
            start_distance = mda_instance.comp_min_distance_between_selections(self.sel1, self.sel2, snapshots=0)
            # check if last one can be accessed through -1
            end_distance = mda_instance.comp_min_distance_between_selections(self.sel1, self.sel2, snapshots=-1)

            if hierarchy_lvl == 0:
                self.assertEqual(float(run_distance[2]), float(end_distance))
                self.assertEqual(float(run_distance[0]), float(start_distance))
            if hierarchy_lvl == 1:
                self.assertEqual(float(run_distance[0][2]), float(end_distance[0]))
                self.assertEqual(float(run_distance[0][0]), float(start_distance[0]))
            if hierarchy_lvl == 2:
                self.assertEqual(float(run_distance[0][0][2]), float(end_distance[0][0]))
                self.assertEqual(float(run_distance[0][0][0]), float(start_distance[0][0]))

    def test_naming(self):
        dataarray = self.mda_project.comp_min_distance_between_selections(self.sel1, self.sel2)
        self.assertEqual(dataarray.dims, ('jobs', 'runs', 'snapshots'))
        self.assertEqual(len(dataarray.jobs), 1)
        self.assertEqual(len(dataarray.runs), 1)
        self.assertEqual(len(dataarray.snapshots), 3)

    def test_structure(self):
        dataarray = self.mda_project.comp_min_distance_between_selections(self.sel1, self.sel2)
        self.assertEqual(type(dataarray), xr.DataArray)

    def test_molecule_selection(self):
        selection_string = self.mda_run.molecule_selection_string(selection='protein')
        id_strings = selection_string.split(' or ')

        self.assertEqual(len(id_strings), self.trypsin_bead_length * 2)


if __name__ == '__main__':
    ut.main(verbosity=2)
