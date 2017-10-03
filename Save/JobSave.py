import pathlib2 as pl
import shutil
import numpy as np
import sys
import tarfile
import PolyLibScan.Tools.config as cfg
import PolyLibScan.Database.db as DB
import parser
import compute

class JobSave(object):

    __git_hash__ = sys.modules['PolyLibScan'].__git_hash__

    def __init__(self, root_path, db_name='jobdata.h5', overwrite=False):
        self.root = pl.Path(root_path)
        setup = self.root.joinpath('config_with_setup.yml')
        self.config = cfg.JobConfig(setup.as_posix())
        self.path = self._set_paths(self.root)
        self.parse = parser.Parser()
        self.db_path = pl.Path(root_path).joinpath(db_name)
        if overwrite:
            self.db = DB.JobDataBase(self.db_path, mode='w')
        else:
            self.db = DB.JobDataBase(self.db_path, mode='a')
        self.db.create_groups()
        self.saved = False
        self.runs = self.read_runs()

    def _set_paths(self, root_path):
        path = {dir_name: root_path.joinpath(dir_name)
                        for dir_name in ['input', 'output', 'logs', 'fifo']}
        path['meta'] = root_path.joinpath('config_with_setup.yml')
        path['p_list'] = root_path.joinpath('particle_list.npy')
        return path

    def save(self):
        self.save_meta_data(self.config)
        self.save_trajectory_meta(self.runs[0].path['full_traj'].as_posix())
        self.save_particle_list(self.path['p_list'].as_posix())
        self.save_runs(self.runs)
        self.save_endstates(self.runs)
        self.saved = True

    def save_trajectory_meta(self, path):
        type_list, step_size = self.parse.trajectory_meta(path)
        traj_info = self.parse.traj_info(path)
        self.db.traj_type_order = type_list
        self.db.traj_info = traj_info

    def save_meta_data(self, config):
        meta_data = self.parse.meta(config)
        self.db.sequence = meta_data['sequence']
        self.db.weights = meta_data['weights']
        self.db.misc = meta_data['misc']
        self.db.active_site = meta_data['active_site']
        self.db.parameter = meta_data['parameter']

    def save_particle_list(self, path):
        particle_data = self.parse.particle_list(path)
        self.db.particle_list = particle_data

    def save_runs(self, runs):
        polymer_ids = np.array(np.unique(self.db.sequence['ID']))
        active_site_pos = self.db.active_site['xyz']
        for run in runs:
            self.db.start_trajectories_save(run.start_traj, run.Id)
            self.db.end_trajectories_save(run.end_traj, run.Id)
            if 'full_traj' in run.path:
                data = self.parse.trajectory(run.path['full_traj'].as_posix())
                self.db.trajectory_save(data['xyz'], run.Id)
            # Save Energy timeseries
            self.db.energie_ts_save(run.energy, run.Id)
            # Save Distance Timeseries
            if 'distance' in run.path:
                self.db.distance_ts_save(run.distance, run.Id)

    def save_endstates(self, runs):
        polymer_ids = np.array(np.unique(self.db.sequence['ID']))
        active_site_pos = self.db.active_site['xyz']
        endstates = np.zeros(len(runs), dtype=[('ID', '>i2'), 
                                               ('Energy', '>f4'), 
                                               ('TotalEnergy', '>f4'), 
                                               ('Distance', '>f4')])
        for run in runs:
            distance = compute.distance_to_active_site(run.end_traj, polymer_ids, 
                                                       active_site_pos)
            endstates[run.Id] = (run.Id, run.energy[-1,0], run.energy[-1,1], distance)
        self.db.end_states = endstates
    
    def read_runs(self):
        runs = []
        input_files = self.path['input'].glob('*')
        return [Run.from_id(self, self.path, input_file.name) 
                for input_file in input_files]

    def get_Error(self):
        if self.has_error():
            with open(self.path['logs'].joinpath('slurm.err').as_posix()) as f:
                return f.read()
        else:
            return 'No error occurred.' 

    def has_error(self):
        '''Check if the job ran through correctly by looking at 
        last line.
        '''
        slurm_error_file = self.path['logs'].joinpath('slurm.err')
        #Check last line
        if slurm_error_file.exists():
            with open(slurm_error_file.as_posix()) as f:
                last_line = f.readlines()[-1].strip()
            if last_line.endswith('Finished simulations. Starting to clean up.Not'):
                return True
        return False

    def save_versions(self, lmp_version=None):
        versions = {'PolyLibScan': str(self.__git_hash__)}
        if lmp_version:
            versions['LAMMPS'] = str(lmp_version)
        else:
            versions['LAMMPS'] = ''
        data = self.parse.version(versions)
        self.db.versions = data

    def remove_local(self):
        shutil.rmtree(self.config.lmp_path['local_root'])

    def clean_up(self, force=False):
        if self.saved or force:
            for run in self.runs:
                run.cleanup()

        if not self.has_error():
            for file_ in ['log.lammps', 'particle_list.npy']:
                if self.root.joinpath(file_).exists():
                    self.root.joinpath(file_).unlink()
            for file_ in ['slurm.out','slurm.err']:
                if self.path['logs'].joinpath(file_).exists():
                    self.path['logs'].joinpath(file_).unlink()    
            for dir_ in ['input', 'output', 'logs', 'fifo']:
                if self.path[dir_].exists():
                    self.path[dir_].rmdir()
            if self.config.sim_parameter['local'] == 1:
                self.remove_local()

class Run(object):

    @classmethod
    def from_id(cls, job, path, str_ID):
        run_id = int(str_ID)
        out_path = path['output']
        in_path = path['input']
        start_file = in_path.joinpath(str_ID)
        energy = out_path.joinpath('Energy'+str_ID)
        start_traj = out_path.joinpath('trajectoryS%s.xyz' % str_ID)
        end_traj = out_path.joinpath('trajectoryE%s.xyz' % str_ID)
        distance_file = out_path.joinpath('distance_as_polymer%s' % str_ID)
        if not distance_file.exists():
            distance_file = None
        traj_file = out_path.joinpath('full_trajectory%s.xyz.gz' % str_ID)
        if not traj_file.exists():
            traj_file = None
        return cls(job, run_id, start_file, energy, start_traj, end_traj, distance_file, traj_file)

    def __init__(self, job, ID, start_file, energy, start_traj, end_traj, distance=None, full_traj=None):
        self.parent = job
        self.Id = ID
        self.path = {}
        self.path['start_file'] = start_file
        self.path['energy'] =  energy
        self.path['start_traj'] = start_traj
        self.path['end_traj'] = end_traj
        if distance:
            self.path['distance'] = distance
            self.distance = self.parent.parse.distance(self.path['distance'].as_posix())
        if full_traj:
            self.path['full_traj'] = full_traj

        self.start_traj = self.parent.parse.xyz(self.path['start_traj'].as_posix())
        self.end_traj = self.parent.parse.xyz(self.path['end_traj'].as_posix())
        self.energy = self.parent.parse.energy(self.path['energy'].as_posix())

    def cleanup(self):
        for path in self.path.values():
            if path.exists():
                path.unlink()

    def __repr__(self):
        return 'LammpsRun - id %d' % self.Id
