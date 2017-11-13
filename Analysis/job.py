#import parser
import numpy as np
import pathlib2 as pl
import pandas as pd
import warnings
import plotting
import bayesModels as bayes
import itertools as it
import pymol_visualisation
import sim_run
import numerics as num_
import PolyLibScan.Database.db as DB

warnings.filterwarnings("ignore")

class Job(bayes.Job):

    def __init__(self, project, db_path, with_pymol=True):
        self.project = project
        self.poly_type = None
        self.db_path = pl.Path(db_path)
        if with_pymol:
            self.pymol = pymol_visualisation.PymolVisJob(self)
        self._parse = DB.JobDataBase(self.db_path, 'r')
        self.Id = None
        self.meta = self._parse.misc
        self.lmp_parameters = {key:val for key,val in self._parse.parameter}
        self.trajectory_meta = self._parse.traj_info
        self.trajectory_order = self._parse.traj_type_order
        self.sequence = self._parse.sequence
        self.weights = dict(self._parse.weights)
        self.active_site = self._parse.active_site
        self._runs = self._read_runs(with_pymol=with_pymol)
        self.particle_ids = self._get_particle_ids()
        self._charge = None
        self._parse.close()
        self._distance_frequency = None
        self._energy_distance_distribution = None
        super(Job, self).__init__()

    def __len__(self):
        return len(self._runs)

    def __getitem__(self, key):
        return self._runs[key]

    def __iter__(self):
        return iter(self._runs)

    def __repr__(self):
        info =  (self.meta['poly_name'], self.meta['protein'], len(self))
        return 'Lammps Run: %s-%s | Runs: %d' % info

    def __str__(self):
        info =  (self.meta['poly_name'], self.meta['protein'], self.Id, len(self))
        weights = zip(*self.weights.items())    
        out = []
        out += ['Lammps Run: %s-%s | ID: %d | Runs: %d' % info]
        out += ['Monomers: %s' % '-'.join(weights[0])]
        out += ['Weights:  %s' % '-'.join(map(str, weights[1]))]
        return '\n'.join(out)

    def _get_particle_ids(self):
        # in case there are no sim runs
        if len(self) == 0:
            raise Exception('No sim data available.')
            
        p_type_ids = {}
        p_type_ids['polymer'] = np.unique(self.sequence['ID'])
        all_particles_set = set(self[0].coordinates()['end']['atom_type'])
        protein_set = all_particles_set - set(p_type_ids['polymer'])
        p_type_ids['protein'] = np.array(list(protein_set))
        return p_type_ids

    def charge():
        doc = "The charge property."
        def fget(self):
            if self.project.parameters:
                if not self._charge:
                    self._charge = self._calculate_polymer_charge()
                return self._charge
            else:
                raise AttributeError('Set parameters by path or dict.')
        def fset(self, value):
            self._charge = value
        def fdel(self):
            del self._charge
        return locals()
    charge = property(**charge())

    def _calculate_polymer_charge(self):
        total_charge = 0
        for monomer in self.sequence['monomer']:
            p_name, sub_name = monomer.split('_')
            total_charge += self.project.parameters['Atoms'][p_name][sub_name]['charge']
        return total_charge

    def _read_runs(self, with_pymol=True):
        '''Reads in all runs from the database and 
        creates a run object. All runs are sorted by id.
        '''
        runs = []
        # get inputs
        for data in self._parse.end_states:
            if 'Distance' in data.dtype.names:
                runs.append(sim_run.Run(self, data['ID'], data['Energy'], data["Distance"], with_pymol=True))
            else:
                runs.append(sim_run.Run(self, data['ID'], data['Energy'], with_pymol=True))
        return sorted(runs, key=lambda x:x.Id)

    def distance_frequency():
        doc = "The distance_frequency property."
        def fget(self):
            if self._energy_distance_distribution == None:
                try:
                    results = self._parse.histogramm
                except DB.tb.NoSuchNodeError:
                    results = self._calc_distance_density(self)
                    hist_array = self._parse_hist_data(results[['distance', 'frequency']],
                                               results['energy'])
                    self._parse.histogram = hist_array
                self._distance_frequency = results[['distance', 'frequency']]
                self._energy_distance_distribution = results[['distance', 'energy']]
            return self._distance_frequency
        def fset(self, value):
            self._distance_frequency = value
        def fdel(self):
            self._distance_frequency = None
        return locals()
    distance_frequency = property(**distance_frequency())

    def energy_distance_distribution():
        doc = "The energy_distance_distribution property."
        def fget(self):
            if self._energy_distance_distribution == None:
                try:
                    results = self._parse.hist_data()
                except parser.DB.tb.NoSuchNodeError:
                    results = self._calc_distance_density(self)
                    hist_array = self._parse_hist_data(results[['distance', 'frequency']],
                                               results[ 'energy'])
                    self._parse.histogram = hist_array
                self._distance_frequency = results[['distance', 'frequency']]
                self._energy_distance_distribution = results[['distance', 'energy']]
            return self._energy_distance_distribution
        def fset(self, value):
            self._energy_distance_distribution = value
        def fdel(self):
            self._energy_distance_distribution = None
        return locals()
    energy_distance_distribution = property(**energy_distance_distribution())

    def _parse_hist_data(self, distance, energy):
        '''Save the observables of distance and energy to the database.

        Since this involves writing to the database, it opens the 
        database in write mode and writes to the /histogram table.
        The timestep information is in the distance argument.
        '''
        
        results_dtype = [('distance', np.float), ('frequency', np.int),
                         ('energy', np.float)]
        data = np.empty(distance.shape[0], results_dtype)
        data['distance'] = distance['distance']
        data['frequency'] = distance['frequency']
        data['energy'] = energy
        return data


    def _calc_distance_density(self, runs):
        dist_v = np.zeros(2000, dtype=[('count', np.int), ('energy', np.float)])
        for run in runs:
            energy = run.binding_energy()[:,0]
            distance = run.distance_time_series()['distance'][-len(energy):]
            # discretize energy and distance timeseries into bins ofr histogram
            distance_density = num_.binning(distance, energy, dist_v)
        reduced_hist = num_.discard_tailing_zeros(dist_v['count'])
        reduced_len = reduced_hist.shape[0]
        energy_mean = dist_v['energy'][:reduced_len] / reduced_hist

        results_dtype = [('distance', np.float), ('frequency', np.int),
                         ('energy', np.float)]
        results = np.empty(reduced_len, results_dtype)
        results['distance'] = np.arange(reduced_len, dtype=np.float)/10
        results['frequency'] = reduced_hist
        results['energy'] = energy_mean
        return results

    def to_dataFrame(self):
        if self._parse.db.is_open():
            self._parse.open()
        data = pd.DataFrame(self._parse.end_states)
        data = data.set_index('ID')
        data.set_index([[self.Id for i in xrange(data.shape[0])], data.index], inplace=True)
        data.index.names = ['PolymerId', 'RunId']
        self._parse.close()
        return data

    def _calc_protein_box(self, margin=20):
        """calculate the minimal box of the protein coordinates, add margin
        and round box_coordinates to the next integer.
        The coordinates are taken from the first timestep; this assumes that
        the protein stays fixes during the simulation.

        Arguments:
            margin: margin added to each side of the box [float/integer]

        Output:
            box -- 2x3 numpy array
        """
        protein_mask = np.in1d(self.trajectory_order, self.particle_ids['protein'])
        trajectory_iterator = self._parse.trajectory_load(0)
        protein_coords = np.array(list(it.islice(trajectory_iterator, len(protein_mask))))[protein_mask]
        box = num_.calc_box(protein_coords, margin=margin)
        discrete_box = np.zeros((2,3))
        discrete_box[0] = np.floor(box[0])
        discrete_box[1] = np.ceil(box[1])
        return discrete_box

    def energies(self):
        return self.project.endstate_matrix.loc[self.Id ,(self.meta['poly_name'], 'Energy')]

    def distances(self):
        return self.project.endstate_matrix.loc[self.Id ,(self.meta['poly_name'], 'Distance')]

    def mean_energy(self):
        return np.mean(self.energies())

