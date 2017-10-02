from PolyLibScan.Tools.container import Container
import pathlib2 as pl
import mock
import numpy as np

import unittest as ut

local_path = pl.Path(__file__).absolute().parent

class TestContainer(ut.TestCase):
    def __init__(self, *args, **kwargs):
        super(TestContainer, self).__init__(*args, **kwargs)
        self.c = Container('test')
        obj = mock.MagicMock()
        self.c['item1'] = obj

    def test_Container_name(self):
        self.assertEqual(self.c.name, 'test')

    def test_generator(self):
        self.assertEqual(self.c['item1'].Id, 1)

    def test_contains(self):
        obj = mock.MagicMock()
        obj.name = 'item12'
        self.c['item12'] = obj
        self.assertTrue(obj.name in self.c)

    def test_iter(self):
        self.assertTrue(len([item for item in self.c])>0)

    def test_unique_key(self):
        obj1 = mock.MagicMock()
        obj2 = mock.MagicMock()
        self.c['test'] = obj1
        self.c['test'] = obj2
        # calculating the auto-generated key name
        new_key = 'test'+ str(self.c['test'].Id+1)
        # looking if new key exists and if the Id is as expected
        self.assertEqual(self.c[new_key].Id, self.c['test'].Id + 1)

    def test_defined_types(self):
        obj1 = mock.MagicMock()
        obj1.name = 'defined_type'
        self.c.define_type('test_define', obj1)
        self.assertEqual(id(self.c['test_define']), id(obj1))
