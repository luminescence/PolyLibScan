from PolyLibScan.Tools.lmp_control import LmpController
import pathlib2 as pl
import unittest as ut
import yaml

local_path = pl.Path(__file__).absolute().parent

class TestControl(ut.TestCase):

    def __init__(self, *args, **kwargs):
        super(TestControl, self).__init__(*args, **kwargs)
        self.Id = 0
        self.parameter = {'stoichiometry': [1,1]}

        self.path = {'input': local_path.joinpath('data', 'jobs', 'abcd', 'input').as_posix(),
                     'output': local_path.joinpath('data', 'jobs', 'abcd', 'output').as_posix(),
                     'root': local_path.joinpath('data', 'jobs', 'abcd').as_posix(),
                     'local_root': ''}

        with open(local_path.joinpath('data', 'parameters_hp2.yml').as_posix()) as f:
            self.parametrisation = yaml.load(f)

    def setUp(self):
        self.controller = LmpController(self.Id, self.parameter, self.path, self.parametrisation, fifos={})

    def tearDown(self):
        self.controller.instance.close()
        del self.controller

    def test_variable(self):
        self.controller.variable('c', 'string', '2')
        self.assertEqual(self.controller.instance.variables['c'].value, 2.0)

    def test_dict_as_lammps_variables(self):
        test_dict = {'a': '2',
                     'b': '10000.0'}

        self.controller.set_dictionary_as_lammps_variables(test_dict, var_style='string')

        for var in test_dict:
            self.assertEqual(self.controller.instance.variables[var].value, float(test_dict[var]))

    def test_list_conversion(self):
        test_list = ['This', 'is', 'my', 'list', 'to', 'test!']
        self.assertEqual(LmpController.convert_python_list_to_lammps_list(test_list), '"This is my list to test!"')
        self.assertEqual(LmpController.convert_python_list_to_lammps_list(test_list, with_quotes=False),
                         'This is my list to test!')


if __name__ == '__main__':
    ut.main(verbosity=2)
