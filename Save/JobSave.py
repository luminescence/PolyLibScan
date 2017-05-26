import pathlib2 as pl
import shutil
import database
import sys
import tarfile
import PolyLibScan.Tools.config as cfg

class JobSave(object):

    __git_hash__ = sys.modules['PolyLibScan'].__git_hash__

    def __init__(self, root_path, db_name='jobdata.h5', overwrite=False):
        self.root = pl.Path(root_path)
        setup = self.root.joinpath('config_with_setup.yml')
        if not setup.exists():
            raise FileNotFoundError('%s does not exist. You have to run sims first' % setup.as_posix())
        self.config = cfg.JobConfig(self.root.joinpath('config_with_setup.yml').as_posix())
        self.path = self._set_paths(self.root)
        self.db_path = pl.Path(root_path).joinpath(db_name)
        # initialize database
        self._db = database.Database(self.db_path, mode='a')
        self._db.create_groups()
        self.saved = False
        self.runs = self._read_runs()

    def _set_paths(self, root_path):
        path = {}
        path['meta'] = root_path.joinpath('config_with_setup.yml')
        path['input'] = pl.Path(self.config.lmp_path['input'])
        path['logs'] = root_path.joinpath('logs/')
        path['output'] = pl.Path(self.config.lmp_path['output'])
        path['fifo'] = pl.Path(self.config.lmp_path['fifo'])
        return path

    def save(self, path=None, overwrite=False):
        self._db.save_meta_data(self.config)
        self._db.save_trajectory_meta(self.runs[0].full_traj.as_posix())
        self._db.save_runs(self.runs)
        self.saved = True

    def _read_runs(self):
        runs = []
        # get inputs
        input_files = self.path['input'].glob('*')
        for input_file in input_files:
            run_id = int(input_file.name)
            energy = self.path['output'].joinpath('Energy'+input_file.name)
            if self.path['output'].joinpath('trajectoryS%s.xyz' % input_file.name).exists():
                start_traj = self.path['output'].joinpath('trajectoryS%s.xyz' % input_file.name)
            else:
                start_traj = None

            if self.path['output'].joinpath('trajectoryE%s.xyz' % input_file.name).exists():
                end_traj = self.path['output'].joinpath('trajectoryE%s.xyz' % input_file.name)
            elif self.path['output'].joinpath('trajectory%s.xyz' % input_file.name).exists():
                end_traj = self.path['output'].joinpath('trajectory%s.xyz' % input_file.name)
            else:
                end_traj = None
            distance_file = self.path['output'].joinpath('distance_as_polymer%s' % input_file.name)
            if not distance_file.exists():
                distance_file = None
            traj_file = self.path['output'].joinpath('full_trajectory%s.xyz.gz' % input_file.name)
            if not distance_file.exists():
                traj_file = None
            if energy.exists() and end_traj != None:
                runs.append(Run(run_id, input_file, energy, 
                                start_traj, end_traj, distance_file, traj_file))
        return runs

    def info(self):
        print('%d Simulations will be saved.' % len(self.runs))

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
        self._db.save_versions(versions)

    def remove_local(self):
        shutil.rmtree(self.config.lmp_path['local_root'])

    def clean_up(self, force=False):
        if self.saved or force:
            for run in self.runs:
                run.energy.unlink()
                run.end_traj.unlink()
                run.input_path.unlink()
                run.start_traj.unlink()
                if run.distance:
                    run.distance.unlink()
                if run.full_traj:
                    run.full_traj.unlink()

        if not self.has_error():
            if self.path['logs'].joinpath('slurm.out').exists():
                self.path['logs'].joinpath('slurm.out').unlink()
            if self.path['logs'].joinpath('slurm.err').exists():
                self.path['logs'].joinpath('slurm.err').unlink()
            if self.root.joinpath('log.lammps').exists():
                self.root.joinpath('log.lammps').unlink()
            if self.root.joinpath('meta.dat').exists():
                self.root.joinpath('meta.dat').unlink()
            if self.path['input'].exists():
                self.path['input'].rmdir()
            if self.path['output'].exists():
                self.path['output'].rmdir() 
            if self.path['logs'].exists():
                self.path['logs'].rmdir()
            if self.config.sim_parameter['local'] == 1:
                self.remove_local()

class Run(object):

    def __init__(self, Id, input_path, energy, start_traj, end_traj, distance=None, full_traj=None):
        self.Id = Id
        self.energy = energy
        self.start_traj = start_traj
        self.end_traj = end_traj
        self.input_path = input_path
        if distance:
            self.distance = distance
        if full_traj:
            self.full_traj = full_traj

    def __repr__(self):
        return 'LammpsRun - id %d' % self.Id
