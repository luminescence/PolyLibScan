import unittest as ut
import PolyLibScan.Database.db as db
import pathlib2 as pl
import os
import numpy as np

local_path = pl.Path(os.path.dirname(os.path.realpath(__file__)))

class Test_database(ut.TestCase):
    '''Test Cases for compute.
    '''

    def __init__(self, *args, **kwargs):
        super(Test_database, self).__init__(*args, **kwargs)
        self.data_folder = local_path.parent
        self.db_file_name = self.data_folder.joinpath('db_test.h5')

    def test_open(self):
        db_h5 = db.Database(self.db_file_name)
        self.assertTrue(db_h5._handle.isopen)
        db_h5.close()
        db_h5 = db.Database(self.db_file_name.as_posix())
        self.assertTrue(db_h5._handle.isopen)
        db_h5.close()

    def test_repr__(self):
        db_h5 = db.Database(self.db_file_name)
        self.assertIsInstance(repr(db_h5), basestring)
        db_h5.close()

    def test_close(self):
        db_h5 = db.Database(self.db_file_name)
        db_h5.close()
        self.assertFalse(db_h5._handle.isopen)

    def test_save_table(self):
        db_h5 = db.Database(self.db_file_name, 'w')
        types = [('a', '<f8'), ('b', '|S10')]
        data = np.array([(3.14, 'huhu'), (2.72, 'asdf')], dtype=types)
        db_h5._save_table(data, '/', 'huhu')
        self.assertTrue('huhu' in db_h5._handle.root)
        db_Data = db_h5._handle.root._f_get_child('huhu')[:]
        self.assertAlmostEqual(db_Data['a'][0], 3.14)
        db_h5.close()

    def test_save_array(self):
        db_h5 = db.Database(self.db_file_name, 'w')
        data = np.ones([10,10])
        db_h5._save_array(data, '/', 'huhu2')
        self.assertTrue('huhu2' in db_h5._handle.root)
        db_Data = db_h5._handle.root._f_get_child('huhu2')[:]
        self.assertAlmostEqual(db_Data[3,8], 1.0)
        db_h5.close()

    def test_create_group(self):
        db_h5 = db.Database(self.db_file_name, 'w')
        db_h5._create_group('test_group', '/')
        self.assertTrue('test_group' in db_h5._handle.root)
        db_h5.close()

    def test_load_table(self):
        db_h5 = db.Database(self.db_file_name, 'w')
        data = np.ones([10,10])
        db_h5._save_array(data, '/', 'huhu2')
        db_data = db_h5._load_table('/', 'huhu2')
        self.assertEqual(db_data[3,3], 1.0)
        db_h5.close()

if __name__ == '__main__':
    ut.main(verbosity=2)
