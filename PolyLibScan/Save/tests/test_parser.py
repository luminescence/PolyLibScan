import unittest as ut
import PolyLibScan.Save as js
import mock
import pathlib2 as pl
import os

local_path = pl.Path(os.path.dirname(os.path.realpath(__file__)))


class Test_Parser(ut.TestCase):
    """Test Cases for JobSave."""

    def __init__(self, *args, **kwargs):
        super(Test_Parser, self).__init__(*args, **kwargs)
        self.parser = js.parser.Parser()

    def test_version(self):
        formatted_array = self.parser.version({'AA': 'a' * 12,
                                               'B': 'b' * 24})

        self.assertEqual(formatted_array.dtype[0], 'S2')
        self.assertEqual(formatted_array.dtype[1], 'S24')


if __name__ == '__main__':
    ut.main(verbosity=2)
