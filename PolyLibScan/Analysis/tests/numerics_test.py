import pathlib2 as pl
import PolyLibScan.Analysis.numerics as numerics
import numpy as np

import unittest as ut

local_path = pl.Path(__file__).absolute().parent

class TestNumerics(ut.TestCase):

	def __init__(self, *args, **kwargs):
		super(TestNumerics, self).__init__(*args, **kwargs)


	def test_calc_box(self):
		test_coordinates1 = np.array([[1.0, 1.0, 1.0],[-1.0, -1.0, -1.0],[0.0, 0.0, 0.0]])
		results1 = np.array([[-1.0, -1.0, -1.0], [1.0, 1.0, 1.0]])
		results2 = np.array([[-11.0, -11.0, -11.0], [11.0, 11.0, 11.0]])
		self.assertTrue(np.all(numerics.calc_box(test_coordinates1, margin=0) == results1))
		self.assertTrue(np.all(numerics.calc_box(test_coordinates1, margin=10) == results2))

	def test_in_box(self):
		test_coordinates_outside = np.array([[-11.0, -0.0, -0.0],[0.0, 10.0, 0.0],[0.0, 0.0, 10.0]])
		test_coordinates_inside = np.array([[0.9, 0.9, 0.9],[0.5, 0.5, 0.5]])
		box = np.array([[-1.0, -1.0, -1.0], [1.0, 1.0, 1.0]])

		self.assertTrue(numerics.in_box(test_coordinates_inside[0], box))
		self.assertTrue(numerics.in_box(test_coordinates_inside[1], box))
		self.assertFalse(numerics.in_box(test_coordinates_outside[0], box))
		self.assertFalse(numerics.in_box(test_coordinates_outside[1], box))
		self.assertFalse(numerics.in_box(test_coordinates_outside[2], box))


if __name__ == '__main__':
    ut.main(verbosity=2)
