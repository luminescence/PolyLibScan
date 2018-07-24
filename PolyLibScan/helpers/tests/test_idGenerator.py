from PolyLibScan.helpers.idGenerator import IdGen
import pathlib2 as pl
import mock
import numpy as np

import unittest as ut

local_path = pl.Path(__file__).absolute().parent


class TestIdGen(ut.TestCase):
    def __init__(self, *args, **kwargs):
        super(TestIdGen, self).__init__(*args, **kwargs)
        self.counter = IdGen()

    def test_missing_counter(self):
        self.assertTrue('missingKeyword' not in self.counter.keys())
        self.assertEqual(self.counter['missingKeyword'], 1)
        self.assertTrue('missingKeyword' in self.counter.keys())

    def test_get_item(self):
        self.assertEqual(self.counter['test_get_item1'], 1)
        self.assertEqual(self.counter['test_get_item1'], 2)
        self.assertEqual(self.counter['test_get_item1'], 3)

    def test_contains(self):
        self.assertEqual(self.counter['contains_test'], 1)
        self.assertTrue('contains_test' in self.counter)


if __name__ == '__main__':
    ut.main(verbosity=2)
