import unittest as ut
import PolyLibScan.Database.db
import pathlib2 as pl
import os
import numpy as np

local_path = pl.Path(os.path.dirname(os.path.realpath(__file__)))

class Test_msms(ut.TestCase):
    '''Test Cases for compute.
    '''

    def __init__(self, *args, **kwargs):
        super(Test_msms, self).__init__(*args, **kwargs)
        pass