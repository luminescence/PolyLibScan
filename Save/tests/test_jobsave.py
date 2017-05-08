import unittest as ut
import PolyLibScan.Save.JobSave as js
import mock 
import pathlib2 as pl
import os
local_path = pl.Path(os.path.dirname(os.path.realpath(__file__)))

class Test_JobSave(ut.TestCase):
	'''Test Cases for JobSave.
	'''

	def __init__(self, *args, **kwargs):
		super(Test_JobSave, self).__init__(*args, **kwargs)
		self.path = local_path.joinpath('data')

	def test_set_paths(self):
		obj = js.JobSave(self.path.as_posix(), db_name='test2.h5')
		results = obj._set_paths(self.path)

		self.assertIn('meta', results)
		self.assertIn('input', results)
		self.assertIn('logs', results)
		self.assertIn('output', results)
		self.assertTrue(results['meta'].exists())
		self.assertTrue(results['input'].exists())
		self.assertTrue(results['logs'].exists())
		self.assertTrue(results['output'].exists())
		obj._db.close()
		obj.db_path.unlink()

	def test_has_error(self):
		obj = js.JobSave(self.path.as_posix(), db_name='test2.h5')
		self.assertTrue(obj.has_error())

		obj._db.close()
		obj.db_path.unlink()

	def test_overwrite(self):
		with self.assertRaises(Exception):
			js_obj = js.JobSave(self.path.as_posix(), db_name='test.h5')

	
if __name__ == '__main__':
    ut.main(verbosity=2)
