import numpy as np
import plotting

class Run(plotting.Run):

    def __init__(self, parent, Id, energy, distance=None):
        self.job = parent
        self._parse = self.job._parse
        self.Id = Id
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
        if not self._parse.is_open():
            self._parse.open()
        data = {'start': self._parse.xyz(self.Id, start=True),
                'end': self._parse.xyz(self.Id)}
        return data

    def _distance_to_active_site(self):
        xyz = self.coordinates()['end']
        mask = np.in1d(xyz['atom_type'], self.job._particle_ids['polymer'])
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
        if not self._parse.is_open():
            self._parse.open()
        return self._parse.energy_series(self.Id, column='binding')

    def total_energy(self):
        if not self._parse.is_open():
            self._parse.open()
        return self.job._parse.energy_series(self.Id, column='total')
        
    def potential_energy(self):
        if not self._parse.is_open():
            self._parse.open()
        return self.job._parse.energy_series(self.Id, column='potential')

    def kinetic_energy(self):
        if not self._parse.is_open():
            self._parse.open()
        return self.job._parse.energy_series(self.Id, column='kinetic')

    def temperature(self):
        if not self._parse.is_open():
            self._parse.open()
        return self.job._parse.energy_series(self.Id, column='temperature')

    def distance_time_series(self):
        if not self._parse.is_open():
            self._parse.open()
        return self._parse.distance_series(self.Id)

    def _dist(self, a, b):
        an = np.array([a['x'], a['y'], a['z']])
        bn = np.array([b['x'], b['y'], b['z']])
        c = an - bn
        return np.sqrt(np.sum(c**2))

    def __repr__(self):
        return 'LammpsRun - id %d' % self.Id