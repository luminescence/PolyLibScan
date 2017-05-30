import parser
import numpy as np
import tqdm
import pathlib2 as pl
import pandas as pd
import collections as col
import yaml
import warnings
import plotting
import bayesModels as bayes
import PolyLibScan.helpers.idGenerator as idGen
import itertools as it
import numerics
import visualize
import concurrent.futures as concurrent
import pymol_visualisation

warnings.filterwarnings("ignore")

class Project(plotting.Project, bayes.Project):
    def __init__(self, project_path, experimental_data=None, parameters=None):
        self.path = pl.Path(project_path)
        self._id_gen = idGen.IdGen()
        self._endstate_matrix = None
        self._parameters = None
        if parameters is not None:
            self.parameters = parameters
        self._visualize = visualize.Visualize(self)
        # stores jobs in list and dict form; the dict stores jobs by polymer type
        self.jobs, self.polymer_types = self.read_jobs(self.path)

        # Experimental Data needs information from job objects
        self._experimental_data = None
        if experimental_data is not None:
            self.experimental_data = experimental_data
            for p_type in self.polymer_types.values():
                p_type.ic50 = self.experimental_data[p_type.polymer_type]

    def read_jobs(self, path):
        '''Read in database files found in folders of project folder.
        The files are read in alphabetical order of their foldername.
        Each job receives an polymer type Id. The Id will stay the same 
        due to the alphabetical sorting.
        The created job objects are stored in a list and in a dict that
        aranges the jobs in lists of polymer types.
        '''
        # creating containers
        jobs = []
        sims_by_type = col.defaultdict(list)

        # sorting makes sure that the 'poly_id' stay the same during read-ins.
        job_path_list = sorted(path.glob('*'), key=lambda x:x.name)
        for folder in tqdm.tqdm(job_path_list, desc='Reading Jobs'):
            if folder.is_dir() and folder.joinpath('jobdata.h5').exists():
                # creating job and setting polymer type id
                job = Job(self, folder.joinpath('jobdata.h5').as_posix())
                poly_id = self._id_gen[job.meta['poly_name']]
                job.Id = poly_id
                # storing job in the two containers
                jobs.append(job)
                sims_by_type[job.meta['poly_name']].append(job)
        polymer_types = {name: PolymerTypeSims(self, sims) for name, sims in sims_by_type.items()}
        return jobs, polymer_types

    def __repr__(self):
        out  = []
        out += ['Project of LAMMPS simulation runs.']
        if 'protein_name' in self.jobs[0].meta:
            out += ['Protein Name: %s (PDB_id: %s).' % (self.jobs[0].meta['protein_name'], 
                                                        self.jobs[0].meta['protein'])]
        else:
            out += ['Protein PdbId: %s.' % self.jobs[0].meta['protein']]
        out += ['Containing %04d Jobs.' % len(self.jobs)]
        out += ['           %04d Polymer Types.' % len(self.polymer_types)]
        out += ['  - having %04d Runs each. (%d per Job)' % (self.endstate_matrix.shape[0],
                                                            len(self.jobs[0]))]
        return '\n'.join(out)

    def polymer_names(self):
        return self.polymer_types.keys()

    def experimental_data():
        doc = "Experimental_data can be set from a DataFrame a Series or a path to a h5 DB-file."
        def fget(self):
            return self._experimental_data
        def fset(self, value):
            if isinstance(value, pd.core.frame.DataFrame):
                try:
                    series = value.loc[:, self.jobs[0].meta['protein']]
                except IndexError:
                    raise IndexError("No Jobs were found in project-folder")
            elif isinstance(value, pd.core.series.Series):
                series = value
            elif isinstance(value, basestring):
                store = pd.HDFStore(value)
                df = store['ic50ext_pdb']
                store.close()
                try:
                    series = df.loc[:, self.jobs[0].meta['protein'].lower()]
                except IndexError:
                    raise IndexError("No Jobs were found in project-folder")
            common_polymers = list(set(self.endstate_matrix.columns.levels[0]) & set(series.index))
            self._experimental_data = series[common_polymers]
        def fdel(self):
            del self._experimental_data
        return locals()
    experimental_data = property(**experimental_data())

    def parameters():
        doc = "Parameters can be set via a path to a yaml file or directly via a dict."
        def fget(self):
            return self._parameters
        def fset(self, value):
            if isinstance(value, basestring):
                with open(value) as f:
                    data = yaml.load(f)
            elif isinstance(value, dict):
                data = value
            else:
                raise ValueError('parameters must be set either with a dict or a path.')
            self._parameters = data
        def fdel(self):
            del self._parameters
        return locals()
    parameters = property(**parameters())

    def endstate_matrix():
        doc = "The endstate_matrix property."
        def fget(self):
            if not isinstance(self._endstate_matrix, pd.core.frame.DataFrame):
                self._endstate_matrix = self.get_distance_matrix()
            return self._endstate_matrix
        def fset(self, value):
            self._endstate_matrix = value
        def fdel(self):
            del self._endstate_matrix
        return locals()
    endstate_matrix = property(**endstate_matrix())

    def get_distance_matrix(self):
        '''construct dataframe with energy and binding distance
        from all simulations.
        '''
        complete_df = pd.concat([p_type.data_frame() 
            for p_type in tqdm.tqdm(self.polymer_types.values(), 
                                    desc='Creating endstate Dataframe')], axis=1)
        return complete_df

class PolymerTypeSims(plotting.PolymerTypeSims, bayes.PolymerTypeSims):
    '''Stores all simulation runs of one polymer type.
    '''
    def __init__(self, project, polymer_sims, ic50=None):
        self.project = project
        self.sims = polymer_sims
        self.register(self.sims)
        self.ic50 = ic50
        self._distance_probability = None
        self._energy_distance_distribution = None
        self.polymer_type = self.sims[0].meta['poly_name']
        self.weights = self.sims[0].weights
        super(PolymerTypeSims, self).__init__()

    def distance_probability():
        doc = "The distance_probability property."
        def fget(self):
            if self._energy_distance_distribution == None:
                results = self._calc_distance_distribution(self.sims)
                self._distance_probability = results[['distance', 'density']]
                self._energy_distance_distribution = results[['distance', 'energy']]
            return self._distance_probability
        def fset(self, value):
            self._distance_probability = value
        def fdel(self):
            self._distance_probability = None
        return locals()
    distance_probability = property(**distance_probability())


    def energy_distance_distribution():
        doc = "The energy_distance_distribution property."
        def fget(self):
            if self._energy_distance_distribution == None:
                results = self._calc_distance_distribution(self.sims)
                self._distance_probability = results[['distance', 'density']]
                self._energy_distance_distribution = results[['distance', 'energy']]
            return self._energy_distance_distribution
        def fset(self, value):
            self._energy_distance_distribution = value
        def fdel(self):
            self._energy_distance_distribution = None
        return locals()
    energy_distance_distribution = property(**energy_distance_distribution())

    def _calc_distance_distribution(self, sims, thread_num=2):
        dist_v = np.zeros(2000, dtype=[('count', np.int), ('energy', np.float)])
        with concurrent.ThreadPoolExecutor(max_workers=thread_num) as executor:
            for job in tqdm.tqdm(sims, desc='Calculating distance distribution of Simulation'):
                job_count = job.distance_frequency['frequency']
                dist_v['count'][:len(job_count)] += job_count
                dist_v['energy'][:len(job_count)] += job.energy_distance_distribution['energy']
        # Shortening the array by trailing zeros
        reduced_hist = numerics.discard_tailing_zeros(dist_v['count'])
        reduced_len = reduced_hist.shape[0]
        
        # Mean Energy
        energy_mean = dist_v['energy'][:reduced_len] / len(self.sims)
        # set elements to zero, who are devided by zero
        energy_mean[np.isnan(energy_mean)] = 0.0

        # Calculating norm
        n = np.sum(reduced_hist).astype(np.float)
        hist_normed = reduced_hist/n
        # results
        results_dtype = [('distance', np.float), ('density', np.float),
                         ('energy', np.float)]
        results = np.empty(reduced_len, results_dtype)
        results['distance'] = np.arange(reduced_len, dtype=np.float)/10
        results['density'] = hist_normed
        results['energy'] = energy_mean
        return results

    def cumulative_binding_probability(self, distance):
        return numerics.cdf(self.distance_probability['density'], distance)

    def cumulative_mean_energy(self, distance):
        mean_energy = (numerics.cdf(self.distance_probability['density'] * self.energy_distance_distribution['energy'], distance)
                                   / self.cumulative_binding_probability(distance))
        return mean_energy

    def distance_cumulative_probability(self):
        '''Cumulative distance function: 
        Return a 2xN Array that denotes the cumulative density within a 
        given radius around the active site.
        '''
        cum_distance = np.empty_like(self.distance_probability)
        cum_distance['distance'] = self.distance_probability['distance'] + 0.1
        cum_distance['density'] = numerics.cumulative_bins(self.distance_probability['density'])
        return cum_distance

    def register(self, children):
        for child in children:
            child.poly_type = self

    def data_frame(self):
        sub_frames = np.empty(len(self.sims), dtype=object)
        for i, job in enumerate(self.sims):
            sub_frames[i] = job.to_dataFrame()
            job._parse.close()
        data_frame = pd.concat(sub_frames)
        col_index = pd.MultiIndex.from_product([[self.polymer_type], data_frame.columns], 
                                                names=['PolyType', 'Results'])
        data_frame.set_axis(1, col_index)
        return data_frame

    def __repr__(self):
        out  = []
        out += ['Polymer Type %s.' % self.polymer_type]
        out += ['Containing %d Jobs.' % len(self.sims)]
        return '\n'.join(out)


class Job(bayes.Job):

    def __init__(self, project, db_path):
        self.project = project
        self.poly_type = None
        self.db_path = pl.Path(db_path)
        self.pymol = pymol_visualisation.PymolVisualisation(self)
        self._parse = parser.Parser(self.db_path, 'r')
        self.Id = None
        self.meta = self._parse.meta('misc')
        self.lmp_parameters = self._parse.lmp_parameters()
        # self.lmp_parameter = self._parse.meta('parameter')
        self.sequence = self._parse.sequence()
        self.weights = self._parse.weights()
        self.active_site = self._parse.active_site()
        self._runs = self._read_runs()
        self._particle_ids = self._get_particle_ids()
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
            total_charge += self.project.parameters['Atoms'][monomer]['charge']
        return total_charge

    def _read_runs(self):
        '''Reads in all runs from the database and 
        creates a run object. All runs are sorted by id.
        '''
        runs = []
        # get inputs
        for data in self._parse.end_state():
            if 'Distance' in data.dtype.names:
                runs.append(Run(self, data['ID'], data['Energy'], data["Distance"]))
            else:
                runs.append(Run(self, data['ID'], data['Energy']))
        return sorted(runs, key=lambda x:x.Id)

    def distance_frequency():
        doc = "The distance_frequency property."
        def fget(self):
            if self._energy_distance_distribution == None:
                try:
                    results = self._parse.hist_data()
                except parser.DB.tb.NoSuchNodeError:
                    results = self._calc_distance_density(self)
                    self._parse.save_hist_data(results,
                                               results)
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
                    self._parse.save_hist_data(results[['distance', 'frequency']],
                                               results[ 'energy'])
                self._distance_frequency = results[['distance', 'frequency']]
                self._energy_distance_distribution = results[['distance', 'energy']]
            return self._energy_distance_distribution
        def fset(self, value):
            self._energy_distance_distribution = value
        def fdel(self):
            self._energy_distance_distribution = None
        return locals()
    energy_distance_distribution = property(**energy_distance_distribution())

    def _calc_distance_density(self, runs):
        dist_v = np.zeros(2000, dtype=[('count', np.int), ('energy', np.float)])
        for run in runs:
            energy = run.binding_energy()
            distance = run.distance_time_series()['distance'][-len(energy):]
            # discretize energy and distance timeseries into bins ofr histogram
            distance_density = numerics.binning(distance, energy, dist_v)
        reduced_hist = numerics.discard_tailing_zeros(dist_v['count'])
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
        if self._parse.is_open():
            self._parse.open()
        data = pd.DataFrame(self._parse.end_state())
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
        protein_mask = np.in1d(self._parse.load_traj_type_order(), self._particle_ids['protein'])
        trajectory_iterator = self._parse.trajectory(0)
        protein_coords = np.array(list(it.islice(trajectory_iterator, len(protein_mask))))[protein_mask]
        box = numerics.calc_box(protein_coords, margin=margin)
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
