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
import pymol_visualisation
import sim_run
import poly_type
import numerics as num_

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
        # stores jobs in list (jobs)
        # and dict form (polymer_types); polymer_types stores jobs sorted by polymer type
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
        polymer_types = {name: poly_type.PolymerTypeSims(self, sims) 
                              for name, sims in sims_by_type.items()}
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

    def _scatter_data(self, with_errors=False, with_labels=False, with_crossvalidation=False, 
                           confidence_interval=0.95, min_dist_to_ac=10):
        if self.experimental_data is not None:
            polymer_list = list(set(self.endstate_matrix.columns.levels[0]) & set(self.experimental_data.index))
            experimental = self.experimental_data[polymer_list]
            experimental.sort_values(inplace=True)
        else:
            polymer_list = self.endstate_matrix.columns.levels[0]
            
        full_distance_matrix = self.endstate_matrix.swaplevel(0,1, axis=1)['Distance'][polymer_list]
        full_energy_matrix = self.endstate_matrix.swaplevel(0,1, axis=1)['Energy'][polymer_list]
        
        distance_matrix = full_distance_matrix[full_distance_matrix<min_dist_to_ac]
        energy_matrix = full_energy_matrix[full_distance_matrix<min_dist_to_ac]

        if with_errors:
            df1 = num_.distance_with_error(distance_matrix)
            df2 = num_.energy_with_error(energy_matrix, 0.95)
            results = pd.concat([df1, df2], axis=1)
        else:
            results = pd.DataFrame(index=distance_matrix.columns)
            results['energy_mean'] = energy_matrix.mean()
            results['dist_mean'] = distance_matrix.count()/float(full_distance_matrix.shape[0])
        # drop polymer types that are not present in the min-radius
        results.dropna(inplace=True)
        if self.experimental_data is not None:
            results['color_by_inhibition'] = experimental[results.index].apply(lambda x:'b' if x>0 else 'r')

        return results


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
                runs.append(sim_run.Run(self, data['ID'], data['Energy'], data["Distance"]))
            else:
                runs.append(sim_run.Run(self, data['ID'], data['Energy']))
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
            energy = run.binding_energy()[:,0]
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
