import tables as tb
import numpy as np
import PolyLibScan.helpers.db as DB
import parser
import compute

class Database(DB.Database):
    """docstring for Database"""
    def __init__(self, path, mode='w'):
        super(Database, self).__init__(path, mode)
        self.parse = parser.Parser()
        
    def create_groups(self):
        '''Creates the default groups in the database.
        '''
        self._create_group('energies', self._handle.root)
        self._create_group('distances', self._handle.root)
        self._create_group('start_trajectories', self._handle.root)
        self._create_group('end_trajectories', self._handle.root)
        self._create_group('meta', self._handle.root)
        self._create_group('trajectory', '/meta')
        self._create_group('trajectories', self._handle.root)

    def get_active_site_data(self):
        return self._load_table('/meta', 'active_site')

    def get_sequence_Data(self):
        return self._load_table('/meta', 'sequence')

    def save_meta_data(self, config):
        meta_data = self.parse.meta(config)
        self._save_table(meta_data['sequence'], '/meta', 'sequence')
        self._save_table(meta_data['weights'], '/meta', 'weights')
        self._save_table(meta_data['misc'], '/meta', 'misc')
        self._save_table(meta_data['active_site'], '/meta', 'active_site')
        self._save_table(meta_data['parameter'], '/meta', 'parameter')
        
    def save_distance_ts(self, path, ID):
        distance_ts = self.parse.distance(path)
        self._save_table(distance_ts, '/distances', 'd%d' % ID)
        return distance_ts

    def save_energy_ts(self, path, ID):
        energy = self.parse.energy(path)
        self._save_array(energy, '/energies', 'e%d' % ID)
        return energy

    def save_trajectory(self, path, kind, ID):
        '''Save trajectory data to one of two tables.
        tables specified in 'kind'.
        
        kind: ['start_trajectories' | 'end_trajectories']
        '''
        traj = self.parse.xyz(path)
        self._save_table(traj, '/%s' % kind, 't%d' % ID)
        return traj

    def save_runs(self, runs):
        endstates = np.zeros(len(runs), dtype=[('ID', '>i2'), 
                                               ('Energy', '>f4'), 
                                               ('TotalEnergy', '>f4'), 
                                               ('Distance', '>f4')])
        polymer_ids = np.array(np.unique(self.get_sequence_Data()['ID']))
        active_site_pos = self.get_active_site_data()['xyz']
        for run in runs:
            # Save First trajectory
            B_traj = self.save_trajectory(run.start_traj.as_posix(), 
                                'start_trajectories', run.Id)
            # Save Last trajectory
            E_traj = self.save_trajectory(run.end_traj.as_posix(), 
                                'end_trajectories', run.Id)
            self.save_complete_trajectory(run.full_traj.as_posix(), run.Id)
            # Save Energy timeseries
            energy = self.save_energy_ts(run.energy.as_posix(), run.Id)
            # Save Distance Timeseries
            if run.distance:
                distance_ts = self.save_distance_ts(run.distance.as_posix(), run.Id)
            # endstate calculations
            distance = compute.distance_to_active_site(E_traj, polymer_ids, 
                                                       active_site_pos)
            endstates[run.Id] = (run.Id, energy[-1,0], energy[-1,1], distance)
        # Save end-states
        self._save_table(endstates, '/', 'end_state')

    def save_trajectory_meta(self, path):
        type_list, step_size = self.parse.trajectory_meta(path)

        self._save_traj_type_order(type_list)
        self._save_traj_info(len(type_list), step_size)

    def _save_traj_type_order(self, data):
        self._save_array(data,'/meta/trajectory/', 'type_order')

    def _save_traj_info(self, particle_number, step_size):
        data = np.array([('particle_number', particle_number),
                         ('step_size', step_size)], dtype=[('key', '|S15'),
                                                           ('value', np.int)])
        self._save_table(data, '/meta/trajectory', 'info')

    def save_complete_trajectory(self, path, id_):
        data = self.parse.trajectory(path)
        self._save_array(data, '/trajectories', 'traj_%d' % id_, compress=True)

    def save_versions(self, version_info):
        np_format = [('program', '|S8'), ('version', '|S40')]
        data = np.array(version_info.items(), dtype=np_format)
        
        self._save_table(data, '/meta', 'versions')
