import os
import PolyLibScan.helpers.config

class JobConfig(PolyLibScan.helpers.config.Config2):
    '''Child of misc config base class.
    '''

    def __init__(self, config_path):
        self.config_path = config_path
        self.read(self.config_path)

    def _copy_data(self, data):
        self.lmp_parameter = data['lammps_parameter']
        self.lmp_path = data['lammps_path']
        self.fifo = data['fifo']

        self.sim_parameter = data['sim_parameter']
        self.sim_path = data['sim_path']

    def save(self, path=None):
        if not path:
            root, ext = os.path.splitext(self.config_path)
            path = root + '_with_setup.yml'
        data = {}
        data['lammps_parameter'] = self.lmp_parameter
        data['lammps_path'] = self.lmp_path
        data['fifo'] = self.fifo
        data['sim_parameter'] = self.sim_parameter
        data['sim_path'] = self.sim_path
        with open(path, 'w') as f:
            PolyLibScan.helpers.config.yaml.dump(data, f)
        return path

    def copy(self):
         new_conf = JobConfig(self.config_path)
         new_conf.lmp_parameter = self.lmp_parameter.copy()
         new_conf.lmp_path = self.lmp_path.copy()
         new_conf.fifo = self.fifo.copy()
         new_conf.sim_parameter = self.sim_parameter.copy()
         new_conf.sim_path = self.sim_path.copy()
         return new_conf
