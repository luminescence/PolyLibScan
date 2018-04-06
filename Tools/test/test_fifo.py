import gzip
import os
import pathlib2 as pl
import time

import PolyLibScan.Tools.fifo as fifo
from PolyLibScan.Tools.fifo_scripts import fifo_traj_compression, fifo_processor

import unittest as ut

local_path = pl.Path(__file__).absolute().parent

class TestFiFo(ut.TestCase):

	def __init__(self, *args, **kwargs):
		super(TestFiFo, self).__init__(*args, **kwargs)
		self.test_fifo =  fifo.BaseFiFo({}, 
								   local_path.joinpath('data', 'fifo_file.fifo').as_posix(),
								   'path.to.nowhere',
								   local_path.joinpath('fifo.out').as_posix(),
								   2000)
		self.test_fifo.terminate()

	def test_additional_arguments(self):
		with self.assertRaises(NotImplementedError):
			self.test_fifo.additional_arguments()		

	def test_create(self):
		self.test_fifo.create()
		self.assertTrue(os.path.exists(self.test_fifo.fifo_path))
		with self.assertRaises(OSError):
			self.test_fifo.create()
		os.remove(self.test_fifo.fifo_path)

	def test_terminate(self):
		with self.assertRaises(OSError):
			self.test_fifo.terminate()
		self.test_fifo.create()
		self.test_fifo.terminate()
		self.assertFalse(os.path.exists(self.test_fifo.fifo_path))

	def test_fifo_exists(self):
		self.assertFalse(self.test_fifo._fifo_exists())
		self.test_fifo.create()
		self.assertTrue(self.test_fifo._fifo_exists())
		self.test_fifo.terminate()
		self.assertFalse(self.test_fifo._fifo_exists())

	def test_path(self):
		new_path = local_path.joinpath('data', 'fifo_2file.fifo').as_posix()
		new_path2 = local_path.joinpath('data', 'fifo_3file.fifo').as_posix()
		
		# creating fifo file via setting fifo_path
		self.test_fifo.fifo_path = new_path
		self.assertTrue(self.test_fifo._fifo_exists())
		self.assertTrue(os.path.exists(new_path))

		# switching path
		self.test_fifo.fifo_path = new_path2
		self.assertFalse(os.path.exists(new_path))
		self.assertTrue(os.path.exists(new_path2))

		# deleting path
		del self.test_fifo.fifo_path
		self.assertFalse(self.test_fifo._fifo_exists())

		# create file with having no path previously
		self.test_fifo.fifo_path = new_path
		self.assertTrue(self.test_fifo._fifo_exists())


class TestCompressionFiFo(ut.TestCase):

	def __init__(self, *args, **kwargs):
		super(TestCompressionFiFo, self).__init__(*args, **kwargs)
		self.compression_fifo = fifo.TrajCompressionFifo({},
														 local_path.joinpath('data', 'example_snapshot.xyz').as_posix(),
														 fifo_traj_compression.__file__,
														 local_path.joinpath('fifo.out').as_posix(),
														 2000)

	def test(self):
		output = 'snapshot.xyz.gz'

		# note that the subprocess takes a moment
		self.compression_fifo._start(output)

		# therefore, wait until output exists
		while not pl.Path(output).exists():
			time.sleep(0.1)

		compressed_content = []
		with gzip.open(output, 'r') as gz_file:
			for line in gz_file:
				compressed_content.append(line)

		# remove output file
		os.remove(output)

		original_content = []
		with open(self.compression_fifo.fifo_path, 'r') as org_file:
			for line in org_file:
				original_content.append(line)

		self.assertEqual(compressed_content, original_content)

class TestDistanceFiFo(ut.TestCase):

	def __init__(self, *args, **kwargs):
		super(TestDistanceFiFo, self).__init__(*args, **kwargs)

		polymer_sequence = [23, 23, 24, 24, 24,
							23, 24, 24, 24, 24,
							24, 23, 24, 24, 24,
							24, 24, 24, 25, 24]
		monomer_ids = list(set(polymer_sequence))

		# dummy needed because distance_fifo.parent.lmp_settings will be checked
		class DummyParent(object):
			def __init__(self):
				self.lmp_settings = {'monomer_ids': monomer_ids,
									 'poly_sequence': polymer_sequence}

		self.distance_fifo = fifo.DistanceFifo(DummyParent(),
											   local_path.joinpath('data','example_snapshot.xyz').as_posix(),
											   fifo_processor.__file__,
											   local_path.joinpath('fifo.out').as_posix(),
											   2000)

	def test(self):
		output = 'output.txt'
		# note that the subprocess takes a moment
		self.distance_fifo._start(output)

		# therefore, wait until output exists
		def subprocess_finished():
			try:
				filesize = os.path.getsize(output)
				if filesize > 0:
					return True
				else: 
					return False
			except OSError:
				return False

		while not subprocess_finished():
			time.sleep(0.1)

		original_content = []
		with open(output, 'r') as org_file:
			for line in org_file:
				original_content.append(line)

		# remove output file
		os.remove(output)

		closest_contact = float(original_content[0].split()[-1])

		self.assertAlmostEqual(closest_contact, 30.4, places=1) # can be confirmed with pymol


if __name__ == '__main__':
    ut.main(verbosity=2)
