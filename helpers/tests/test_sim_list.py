import unittest as ut
import PolyLibScan.helpers.sim_list as h_sim_list
import pathlib2 as pl
import os


local_path = pl.Path(os.path.dirname(os.path.realpath(__file__)))
data_folder = local_path.joinpath('data')

class Test_msms(ut.TestCase):
    '''Test Cases for compute.
    '''

    def __init__(self, *args, **kwargs):
        super(Test_msms, self).__init__(*args, **kwargs)
        self.sim_list_file1 = data_folder.joinpath('sim.list1')
        if self.sim_list_file1.exists():
        	self.sim_list_file1.unlink()
        self.sim_list_file2 = data_folder.joinpath('sim.list2')
        
    def test_create_empty(self):
    	sim_list = h_sim_list.SimList(self.sim_list_file1.as_posix(), '/data/test', 10)

    	self.assertTrue(self.sim_list_file1.exists())

    def test_mark_complete(self):
    	sim_list = h_sim_list.SimList(self.sim_list_file1.as_posix(), '/data/test', 10)
    	sim_list.mark_complete(0)

    	# check if the sim was saved interally
    	self.assertTrue(len(sim_list.sims) == 1)
    	self.assertEqual(sim_list.sims[0].id, 0)

    	# check if the sim was saved in the file
    	sim_list2 = h_sim_list.SimList(self.sim_list_file1.as_posix(), '/data/test', 10)
    	self.assertTrue(len(sim_list2.sims) == 1)
    	self.assertEqual(sim_list2.sims[0].id, 0)

    def test_iter(self):
    	sim_list = h_sim_list.SimList(self.sim_list_file1.as_posix(), '/data/test', 3)
    	self.assertEqual(list(sim_list), [0,1,2])
    	sim_list.mark_complete(1)
    	self.assertEqual(list(sim_list), [0,2])


    def tearDown(self):
    	if self.sim_list_file1.exists():
        	self.sim_list_file1.unlink()


if __name__ == '__main__':
    ut.main(verbosity=2)
