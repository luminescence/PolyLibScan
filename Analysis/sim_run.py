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

    def trajectory(self, molecule='full'):
        particle_order = self.job.trajectory_order
        ids = self.job.particle_ids
        if molecule == 'protein':
            type_filter = AtomFilter(particle_order, sim=self.job,
                                                  monomer_id=ids['protein'], molecule=molecule)
        elif molecule == 'polymer':
            type_filter = AtomFilter(particle_order, sim=self.job,
                                                  monomer_id=ids['polymer'], molecule=molecule)
        elif molecule == 'full':
            all_ids = np.concatenate((ids['protein'], ids['polymer']))
            type_filter = AtomFilter(particle_order, sim=self.job,
                                                  monomer_id=all_ids, molecule=molecule)
        else:
            raise AttributeError("molecule must be 'protein', 'polymer' or 'full'.")
        traj_iterator = self.job._parse.trajectory_load(self.Id)

        monomer_coords = it.ifilter(type_filter.filter_function, traj_iterator)
        seq_length = len(self.sequence())
        time_step_coords = np.zeros(seq_length, dtype=[('type_', np.int), ('xyz', np.float,3)])
        full_time_steps = self.job.lmp_parameters['time_steps']/self.job.trajectory_meta['step_size']+1
        for time_step in xrange(full_time_steps):
            for i, atom_id, coord in it.izip(xrange(seq_length), self.sequence()['ID'], monomer_coords):
                time_step_coords[i][0] = atom_id
                time_step_coords[i][1] = coord
            yield time_step_coords

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


class AtomFilter(object):

    def __init__(self, particle_order, sim, monomer_id, molecule='full', filter_specification='type'):
        self.particle_order = particle_order
        self.sequence = sim.sequence
        self.monomer_id = monomer_id
        self.molecule = molecule
        self.filter_specification = filter_specification
        self.mask = self.generate_mask()
        self.filter_function = self.create_atom_type_filter()

    def generate_mask(self):
        """return a 1d np.array (lenghth = number of particles in MD) filled with True or False"""

        polymer_length = len(self.sequence)
        protein_length = len(self.particle_order) - polymer_length

        if self.filter_specification == 'type':
            mask = np.in1d(self.particle_order, self.monomer_id)
        elif self.filter_specification == 'id':
            if self.molecule != 'polymer':
                raise NotImplementedError('Filtering by id is only implemented for the polymer so far! ')
            mask = np.array([False] * (polymer_length + protein_length))
            absolute_id = [x+protein_length for x in self.monomer_id]
            mask[absolute_id] = True
        else:
            raise NotImplementedError("Filter is unknown. Use 'type' or 'id'!")

        # if molecule == 'full', nothing needs to be done
        if self.molecule == 'polymer':
            mask[:protein_length] = [False] * protein_length
        elif self.molecule == 'protein':
            mask[protein_length:] = [False] * polymer_length

        return mask

    def create_atom_type_filter(self):
        iterator = it.cycle(self.mask)

        def atom_type_filter(dummy_variable=True):
            # one variable is required but there is nothing to be passed
            return iterator.next()

        return atom_type_filter
