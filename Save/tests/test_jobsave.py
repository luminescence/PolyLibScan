import unittest as ut
import PolyLibScan.Save as js
import mock 
import pathlib2 as pl
import os

local_path = pl.Path(os.path.dirname(os.path.realpath(__file__)))

class Test_JobSave(ut.TestCase):
	'''Test Cases for JobSave.
	'''

	def __init__(self, *args, **kwargs):
		super(Test_JobSave, self).__init__(*args, **kwargs)
		self.path = {'input': local_path.joinpath('data', 'jobs', 'abcd', 'input').as_posix(),
		'output': local_path.joinpath('data', 'jobs', 'abcd', 'output').as_posix(),
		'root': local_path.joinpath('data', 'jobs', 'abcd').as_posix(),
		'script': local_path.joinpath('data', 'static', 'cluster.in').as_posix(),
		'logs': local_path.joinpath('data', 'jobs', 'abcd', 'logs').as_posix(),
		'fifo': local_path.joinpath('data', 'jobs', 'abcd', 'fifo').as_posix(),
		'local_root': ''}
		self.db_file1 = pl.Path(self.path['root']).joinpath('test.h5')
		self.db_file2 = pl.Path(self.path['root']).joinpath('test2.h5')

	def test_set_paths(self):
		if self.db_file2.exists():
			self.db_file2.unlink()
		obj = js.JobSave(self.path, db_name='test2.h5')
		results = obj._set_paths(self.path)
		self.assertIn('meta', results)
		self.assertIn('input', results)
		self.assertIn('logs', results)
		self.assertIn('output', results)
		self.assertIn('fifo', results)
		self.assertTrue(results['meta'].exists())
		self.assertTrue(results['input'].exists())
		self.assertTrue(results['logs'].exists())
		self.assertTrue(results['output'].exists())
		self.assertTrue(results['fifo'].exists())
		obj.db.close()
		obj.db_path.unlink()


	def test_run(self):
		if self.db_file2.exists():
			self.db_file2.unlink()
		js_obj = js.JobSave(self.path, db_name='test2.h5')
		self.assertEqual(len(js_obj.runs), 1)
		js_obj.db.close()
		self.db_file2.unlink()

	def test_get_version(self):
		if self.db_file2.exists():
			self.db_file2.unlink()
		js_obj = js.JobSave(self.path, db_name='test2.h5')
		self.assertTrue(len(js_obj.__git_hash__), 20)

	# def test_has_error(self):
	# 	#if self.db_file2.exists():
	# 	#	self.db_file2.unlink()
	# 	obj = js.JobSave(self.path, db_name='test2.h5')
	# 	self.assertFalse(obj.has_error())

	# 	obj.db.close()
	# 	obj.db_path.unlink()

	# def test_overwrite(self):
	# 	with self.assertRaises(Exception):
	# 		js_obj = js.JobSave(self.path, db_name='test.h5')

	def test_git_hashes(self):
		'''since the hashes and versions will change in the future,
		we will just check for failure.
		'''
		if self.db_file2.exists():
			self.db_file2.unlink()
		js_obj = js.JobSave(self.path, db_name='test2.h5')
		self.assertNotEqual(js_obj.__git_hash__, ' ')
		js_obj.db.close()
		js_obj.db_path.unlink()

if __name__ == '__main__':
    ut.main(verbosity=2)
