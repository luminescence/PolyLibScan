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
import PolyLibScan.Tools.config as cfg

warnings.filterwarnings("ignore")

class Job(bayes.Job):

    def __init__(self, project, db_path, with_pymol=True):
        self.project = project
        self.poly_type = None
        self.db_path = pl.Path(db_path)
        self._parse = DB.JobDataBase(self.db_path, 'r')
        if with_pymol:
            self.pymol = pymol_visualisation.PymolVisJob(self, protein_path=self.project.protein_path)
        self.Id = None
        self.meta = self._parse.misc
        self.lmp_parameters = {key:val for key,val in self._parse.parameter}
        self.trajectory_meta = self._parse.traj_info
        self.trajectory_order = self._parse.traj_type_order
        self.particle_list = self._parse.particle_list
        self.sequence = self._parse.sequence
        self.weights = dict(self._parse.weights)

        # set self.active_site if there's a protein present
        self.config = cfg.JobConfig(self.db_path.parent.joinpath('config_with_setup.yml').as_posix())
        if 'stoichiometry' in self.config.sim_parameter:
            if self.config.sim_parameter['stoichiometry'][0] > 0:
                self.active_site = self._parse.active_site
        else:
            self.active_site = self._parse.active_site

        self._runs = self._read_runs(with_pymol=with_pymol)
        self.particle_ids = self._get_particle_ids()
        self._charge = None
        self._hydrophobicity = None
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
        out  = []
        out += ['Lammps Run: %s-%s | ID: %d | Runs: %d' % info]
        out += ['Monomers: %s' % '-'.join(weights[0])]
        out += ['Weights:  %s' % '-'.join(map(str, weights[1]))]
        return '\n'.join(out)

    def _get_particle_ids(self):
        # in case there are no sim runs
        if len(self) == 0:
            raise Exception('No sim data available.')

        # this approach works if:
        # 1. there is only one protein and one polymer molecule
        # 2. the protein was created first, the polymer second
        p_type_ids = {}
        p_type_ids['polymer'] = np.unique(self.sequence['ID'])
        polymer_length = len(self.sequence)
        p_type_ids['protein'] = np.unique(self.trajectory_order[:-polymer_length])

        return p_type_ids

    def hydrophobicity():
        doc = "The total hydrophobicity."
        def fget(self):
            if self.project.parameters:
                if self._hydrophobicity is None:
                    self._hydrophobicity = self._calculate_polymer_property('hydrophobicity')
                return self._hydrophobicity
            else:
                raise AttributeError('Set parameters by path or dict.')
        def fset(self, value):
            self._hydrophobicity = value
        def fdel(self):
            del self._hydrophobicity
        return locals()
    hydrophobicity = property(**hydrophobicity())

    def charge():
        doc = "The total charge."
        def fget(self):
            if self.project.parameters:
                if self._charge is None:
                    self._charge = self._calculate_polymer_property('charge')
                return self._charge
            else:
                raise AttributeError('Set parameters by path or dict.')
        def fset(self, value):
            self._charge = value
        def fdel(self):
            del self._charge
        return locals()
    charge = property(**charge())

    def _calculate_polymer_property(self, property_):
        total = []
        for monomer in self.sequence['monomer']:
            p_name, sub_name = monomer.split('_')
            total.append(self.project.parameters['Atoms'][p_name][sub_name][property_])
        if property_ == 'charge':
            return sum(total)
        elif property_ == 'hydrophobicity':
            return sum(filter(lambda x:x>0, total))

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

    def distance_histogram(sel_attr):
        def calc_both(self):
            try:
                results = self._parse.histogramm
            except DB.tb.NoSuchNodeError:
                results = self._calc_distance_density(self)
                hist_array = self._parse_hist_data(results[['distance', 'frequency']],
                                                   results['energy'])
                self._parse.histogram = hist_array
            self._distance_frequency = results[['distance', 'frequency']]
            self._energy_distance_distribution = results[['distance', 'energy']]

        def fget(self, selected_attribute=sel_attr):
            if getattr(self, selected_attribute) is None:
                calc_both(self)
            return getattr(self, selected_attribute)

        def fset(self, value, selected_attribute=sel_attr):
            setattr(self, selected_attribute, value)

        def fdel(self, selected_attribute=sel_attr):
            setattr(self, selected_attribute, None)

        export_dict = locals()
        del export_dict['calc_both']
        del export_dict['sel_attr']

        return export_dict

    energy_distance_distribution = property(**distance_histogram(sel_attr='_energy_distance_distribution'))
    distance_frequency = property(**distance_histogram(sel_attr='_distance_frequency'))

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
        # 2000 is an arbitrary choice for the number of bins. 
        # Because the length of each bin is 0.1 Angstrom, 
        # It assumes that the distance never exceeds 200 Angstrom.
        bins = 2000
        bin_length = 0.1
        dist_v = np.zeros(bins, dtype=[('count', np.int), ('energy', np.float)])
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
        results['distance'] = np.arange(reduced_len, dtype=np.float) * bin_length
        results['frequency'] = reduced_hist
        results['energy'] = energy_mean
        return results

    def to_dataFrame(self, timestep=None):
        if self._parse.db.is_open():
            self._parse.open()
        if timestep is None:
            data = pd.DataFrame(self._parse.end_states)
            data = data.set_index('ID')
            self._parse.close()
        else:
            data = self.create_state_df_on(timestep)
        data.set_index([[self.Id for i in xrange(data.shape[0])], data.index], inplace=True)
        data.index.names = ['PolymerId', 'RunId']
        return data

    def create_state_df_on(self, timestep):
        values = np.zeros(len(self),
                          dtype=[('Energy', np.float), 
                                 ('TotalEnergy', np.float),
                                 ('Distance', np.float)])
        # binding energy
        for i,run in enumerate(self):
            b_ene = run.binding_energy()
            values['Energy'][i] = b_ene[[b_ene[:,1]== timestep]][0,0]
        # total_energy
        for i,run in enumerate(self):
            total_ene = run.total_energy()
            values['TotalEnergy'][i] = total_ene[[total_ene[:,1]== timestep]][0,0]
        # distance
        for i,run in enumerate(self):
            dist = run.distance_time_series()
            values['Distance'][i] = dist[[dist['time_step']== timestep]]['distance']
        data = pd.DataFrame(values)
        return data

    def gen_intermediate_state(self):
        pass

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

    def get_particle_parameters(self, parameter):
        """return a list of the requested parameter for all particles in the system"""
        particle_types = self.particle_list['name']
        atom_possibilities = self.project.parameters['Atoms']
        monomer_types = []  # amino acids, BP, BA, ...
        bead_types = []  # bb or sc
        for x in particle_types:
            # some monomers have NAME_bb, refering to the back bone bead
            # amino acids don't have this appendix so we need try/except
            try:
                monomer, bead = x.split('_')
            except ValueError:
                monomer = x
                bead = 'bb'
            # ghost atoms are not in the trajectories so they need to be neglected
            if monomer != 'ghost':
                monomer_types.append(monomer)
                bead_types.append(bead)
        particle_parameter = [atom_possibilities[mono][pos][parameter] for mono, pos in zip(monomer_types, bead_types)]
        return particle_parameter
