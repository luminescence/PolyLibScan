import MDAnalysis as mda
import os
import pandas as pd
import tqdm


class MdaRun(object):

    def __init__(self, run):
        self.run = run
        self.job = run.job
        self.initial_xyz_file = 'starting_point.xyz'

    def temporarily_provide_xyz(self, input_generator):
        """generate required xyz file, pass generator, remove xyz file"""

        # xyz file of the starting positions is required for mda
        starting_coordinates = self.job._parse.db._load_table('/start_trajectories', 't%s' % self.run.Id)
        self.generate_xyz_file_from_array(starting_coordinates, filename=self.initial_xyz_file)

        # receiving a generator of a generator, iterating to give back generators
        for contained_generator in input_generator:
            yield contained_generator

        os.remove(self.initial_xyz_file)

    @staticmethod
    def generate_xyz_file_from_array(snapshot, filename):
        with open(filename, 'w') as file_handle:
            no_atoms = len(snapshot)
            snapshot.dtype = [('id', int), ('x', float), ('y', float), ('z', float)]

            snapshot_df = pd.DataFrame(snapshot[['x', 'y', 'z']], index=snapshot['id'])

            file_handle.write(str(no_atoms) + '\n')
            file_handle.write('generic line \n')
            snapshot_df.to_csv(file_handle, header=False, sep=' ')

    def assign_masses_in_universe(self, mda_universe):
        weights = self.get_particle_weights()

        atom_group = mda_universe.select_atoms('all')
        atom_group.masses = weights

    def get_particle_weights(self):
        particle_types = self.job._parse.db._load_table('/meta', 'particle_list')['name']
        atom_possibilities = self.job.project.parameters['Atoms']
        monomer_types = []  # amino acids, BP, BA, ...
        bead_types = []  # bb or sc
        for x in particle_types:
            monomer, bead = x.split('_')
            monomer_types.append(monomer)
            bead_types.append(bead)
        weights = [atom_possibilities[mono][pos]['mass'] for mono, pos in zip(monomer_types, bead_types)]
        return weights

    def universe_creator(self, snapshot):
        no_particles = len(self.job.trajectory_order)
        mda_universe = mda.Universe(self.initial_xyz_file, snapshot['xyz'].reshape(1, no_particles, 3))
        self.assign_masses_in_universe(mda_universe)
        return mda_universe

    def stream_trajectory_iterator(self, func):
        for x in tqdm.tqdm_notebook(self.run.trajectory()):
            mda_universe = self.universe_creator(x)
            yield func(mda_universe)

    @staticmethod
    def trajectory_capturer(func):
        """store everything from trajectory in pd.Series"""
        output_list = []
        for x in func:
            output_list.append(x)

        return pd.Series(output_list)

    def compute_radius_of_gyration(self):
        """wrap function for universe; can't be done otherwise because self. is not defined in class"""
        @self.trajectory_capturer
        @self.temporarily_provide_xyz
        @self.stream_trajectory_iterator
        def rg_from_universe(mda_universe):
            return mda_universe.select_atoms('all').radius_of_gyration()

        return rg_from_universe


class MdaJob(object):

    def __init__(self, job):
        self.job = job
        self.MdaRuns = [MdaRun(x) for x in self.job]

    def func_on_all_runs(self, func):
        all_results = []

        for sel_MdaRun in tqdm.tqdm_notebook(self.MdaRuns):
            all_results.append(getattr(sel_MdaRun, func)())

        return pd.DataFrame(all_results).transpose()
