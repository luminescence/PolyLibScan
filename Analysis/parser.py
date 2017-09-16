import json
import numpy as np
import PolyLibScan.Database.db as DB

class Parser(DB.Database):

    def __init__(self, path, mode='r'):
        super(Parser, self).__init__(path, mode)

    def end_state(self):
        if not self.is_open():
            self.open()
        if 'end_state' in self._handle.root:
            data = self._load_table('/', 'end_state')
        else:
            data = self._load_table('/', 'end_energies')
        return data

    def energy_series(self, Id, column='binding'):
        if column == 'binding':
            try:
                return self._load_ctable('/energies', 'e%d' % Id, col=[0,5])
            except IndexError:
                return self._load_ctable('/energies', 'e%d' % Id, col=0)
        elif column == 'total':
            return self._load_ctable('/energies', 'e%d' % Id, col=[1,5])
        elif column == 'potential':
            return self._load_ctable('/energies', 'e%d' % Id, col=[2,5])
        elif column == 'kinetic':
            return self._load_ctable('/energies', 'e%d' % Id, col=[3,5])
        elif column == 'temperature':
            return self._load_ctable('/energies', 'e%d' % Id, col=[4,5])

    def distance_series(self, Id):
        return self._load_ctable('/distances', 'd%d' % Id)

    def meta(self, table_name='misc'):
        data = self._load_table('/meta', table_name)
        meta_data = {key: value for key,value in data}
        if table_name == 'misc':
            if 'polymer_name' in meta_data:
                meta_data['poly_name'] = meta_data['polymer_name'].replace('[x]', '')
        return meta_data

    def active_site(self):
        return self._load_table('/meta', 'active_site')

    def weights(self):
        data= self._load_table('/meta', 'weights')
        return {key: value for key,value in data}

    def sequence(self):
        return self._load_table('/meta', 'sequence')

    def xyz(self, Id, start=False):
        if start:
            group = self._handle.root.start_trajectories
        else:
            group = self._handle.root.end_trajectories
        return self._load_table(group, 't%d' % Id)

    def traj_meta(self):
        data = self._load_table('/meta/trajectory', 'info')
        meta_data = {key: value for key,value in data}
        return meta_data

    def load_traj_type_order(self):
        if not self.is_open():
            self.open()
        return self._load_ctable('/meta/trajectory/', 'type_order')

    def trajectory(self, Id):
        if not self.is_open():
            self.open()
        return self._load_ctable('/trajectories', 'traj_%d' % Id, as_iterator=True)

    def save_hist_data(self, distance, energy):
        '''Save the observables of distance and energy to the database.

        Since this involves writing to the database, it opens the 
        database in write mode and writes to the /histogram table.
        The timestep information is in the distance argument.
        '''
        if self.is_open():
            self.close()
        self.open(mode='a')
        results_dtype = [('distance', np.float), ('frequency', np.int),
                         ('energy', np.float)]
        data = np.empty(distance.shape[0], results_dtype)
        data['distance'] = distance['distance']
        data['frequency'] = distance['frequency']
        data['energy'] = energy['energy']
        self._save_table(data, '/', 'histogramm')
        self.close()

    def hist_data(self):
        if not self.is_open():
            self.open()
        hist_data = self._load_table('/', 'histogramm')
        self.close()
        return hist_data

    def lmp_parameters(self):
        if not self.is_open():
            self.open()
        try:
            # this table does not exist in legacy sims
            data = {key:val for key,val in self._load_table('/meta', 'parameter')}
        except DB.tb.NoSuchNodeError:
            # take default. This is sensible since old simulations had the same settings.
            data = {'dielectric_par': 78, 'timestep': 1}
        return data

