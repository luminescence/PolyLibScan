import numpy as np
import plotting
import pymol_visualisation as pym
import itertools as it

class Run(plotting.Run):

    def __init__(self, parent, Id, energy, distance=None, with_pymol=True):
        self.job = parent
        self._parse = self.job._parse
        self.Id = Id
        if with_pymol:
            self.pymol = pym.PymolVisRun(self)
        self._distance = -1.0
        if not distance:
            self._distance = distance
        self.energy = energy

    def distance():
        doc = "The distance property."
        def fget(self):
            if self._distance < 0.0:
                self._distance = self._distance_to_active_site()
            return self._distance
        def fset(self, value):
            self._distance = value
        return locals()
    distance = property(**distance())

    def show_trajectory(self):
        self.job.project._visualize.trajectory(self)

    def meta(self):
        return self.job.meta

    def weights(self):
        return self.job.weights

    def sequence(self):
        return self.job.sequence

    def coordinates(self):
        if not self._parse.db.is_open():
            self._parse.open()
        data = {'start': self._parse.start_trajectories_load(self.Id),
                'end': self._parse.end_trajectories_load(self.Id)}
        return data

    def _create_atom_type_filter(self, particle_order, monomer_id=None):
        if monomer_id == None:
            monomer_id = self.monomer_id
        mask = np.in1d(particle_order, monomer_id)
        iterator = it.cycle(mask)
        def atom_type_filter(whatever):
            return iterator.next()
        return atom_type_filter

    def trajectory(self, molecule='full'):
        particle_order = self.job.trajectory_order
        ids = self.job.particle_ids
        if molecule == 'protein':
            type_filter = self._create_atom_type_filter(particle_order, 
                            monomer_id=ids['protein'])
        if molecule == 'polymer':
            type_filter = self._create_atom_type_filter(particle_order, 
                            monomer_id=ids['polymer'])
        if molecule == 'full':
            alll_ids = np.concatenate(ids['protein'], ids['polymer'])
            type_filter = self._create_atom_type_filter(particle_order, 
                            monomer_id=ids)
        else:
            raise AttributeError("molecule must be 'protein', 'polymer' or 'full'.")
        traj_iterator = self.job._parse.trajectory_load(self.Id)

        monomer_coords = it.ifilter(type_filter, traj_iterator)
        time_step_coords = np.zeros(len(self.sequence()), dtype=[('type_', np.int), ('xyz', np.float,3)])
        full_time_steps = self.job.lmp_parameters['time_steps']/self.job.trajectory_meta['step_size']+1
        for time_step in xrange(full_time_steps):
            for i, coord in it.izip(xrange(len(self.sequence())), monomer_coords):
                time_step_coords[i] = coord
            yield time_step_coords

    def polymer_trajectory(self):
        particle_order = self.job.trajectory_order
        type_filter = self._create_atom_type_filter(particle_order, 
                            monomer_id=self.job.particle_ids['polymer'])
        traj_iterator = self.job._parse.trajectory_load(self.Id)

        monomer_coords = it.ifilter(type_filter, traj_iterator)
        xyz = np.zeros(len(self.sequence()), dtype=[('xyz', np.float,3)])
        full_time_steps = self.job.lmp_parameters['time_steps']/self.job.trajectory_meta['step_size']+1
        for time_step in xrange(full_time_steps):
            for i, coord in it.izip(xrange(len(self.sequence())), monomer_coords):
                xyz[i] = coord
            yield xyz

    def _distance_to_active_site(self):
        xyz = self.coordinates()['end']
        mask = np.in1d(xyz['atom_type'], self.job.particle_ids['polymer'])
        poly_coords = xyz[mask]
        # subtracting the active site id by one since lammps starts 
        # counting from one, while Numpy does not.
        active_site = xyz[self.job.active_site['xyz']-1]
        min_dist = 1000
        for resi in active_site:
            for mono in poly_coords:
                md = self._dist(resi, mono)
                if md < min_dist:
                    min_dist = md
        return min_dist

    def binding_energy(self):
        if not self._parse.db.is_open():
            self._parse.open()
        return self._parse.energie_ts_load(self.Id, col=[0,5])

    def total_energy(self):
        if not self._parse.db.is_open():
            self._parse.open()
        return self.job._parse.energie_ts_load(self.Id, col=[1,5])
        
    def potential_energy(self):
        if not self._parse.db.is_open():
            self._parse.open()
        return self.job._parse.energie_ts_load(self.Id, col=[2,5])

    def kinetic_energy(self):
        if not self._parse.db.is_open():
            self._parse.open()
        return self.job._parse.energie_ts_load(self.Id, col=[3,5])

    def temperature(self):
        if not self._parse.db.is_open():
            self._parse.open()
        return self.job._parse.energie_ts_load(self.Id, col=[4,5])

    def distance_time_series(self):
        if not self._parse.db.is_open():
            self._parse.open()
        return self._parse.distance_ts_load(self.Id)

    @staticmethod
    def _dist(a, b):
        an = np.array([a['x'], a['y'], a['z']])
        bn = np.array([b['x'], b['y'], b['z']])
        c = an - bn
        return np.sqrt(np.sum(c**2))

    def __repr__(self):
        return 'LammpsRun - id %d' % self.Id
