from functools import partial
import inspect
import numpy as np
import os

import MDAnalysis as mda
import pandas as pd
from scipy.spatial.distance import cdist
import tqdm
import xarray as xr


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
        weights = [atom_possibilities[mono][pos]['mass'] for mono, pos in zip(monomer_types, bead_types)]
        return weights

    def universe_creator(self, snapshot):
        no_particles = len(self.job.trajectory_order)
        mda_universe = mda.Universe(self.initial_xyz_file, snapshot['xyz'].reshape(1, no_particles, 3))
        self.assign_masses_in_universe(mda_universe)
        return mda_universe

    def stream_trajectory_iterator(self, func, snapshots='all'):
        snapshots_range = self.determine_snapshots_range(snapshots)

        for ts, x in tqdm.tqdm_notebook(enumerate(self.run.trajectory())):
            if ts in snapshots_range:
                mda_universe = self.universe_creator(x)
                yield func(mda_universe)

    def determine_snapshots_range(self, snapshots):
        """:return range of snapshots; convert multiple formats into range.
        Note that the order of the input is lost completely!"""
        # get trajectory info
        timesteps_per_snapshot = dict(self.job._parse.db._load_table('/meta/trajectory', 'info'))['step_size']
        total_no_timesteps = dict(self.job._parse.db._load_table('/meta', 'parameter'))['time_steps']
        no_snapshots = (total_no_timesteps / timesteps_per_snapshot) + 1  # include 0th time step

        # convert to range
        if snapshots == 'all':
            snapshots_range = range(0, no_snapshots)
        elif type(snapshots) == int:
            snapshots_range = [snapshots]
        elif type(snapshots) == list or type(snapshots) == tuple:
            snapshots_range = snapshots
        else:
            raise ValueError('Could not understand snapshots given!!!')

        # account for counting backwards, i.e. last timestep == -1
        snapshots_range = [no_snapshots + x if x < 0 else x for x in snapshots_range]
        return snapshots_range

    def generator_to_series(self, input_generator, snapshots):
        """store everything from generator (in our case trajectory) in pd.Series"""
        output_list = []
        for x in input_generator:
            output_list.append(x)

        series = pd.Series(output_list)
        snapshots_as_indeces = self.determine_snapshots_range(snapshots)
        snapshots_as_indeces = sorted(snapshots_as_indeces)
        series.index = snapshots_as_indeces

        series.index.rename('snapshots', inplace=True)

        return series

    # scientifically meaningful methods
    # start with comp_

    def comp_radius_of_gyration(self, snapshots='all'):
        """wrap function for universe"""
        def rg_from_universe(mda_universe):
            return mda_universe.select_atoms('all').radius_of_gyration()

        return self.generator_to_series(self.temporarily_provide_xyz(
            self.stream_trajectory_iterator(rg_from_universe, snapshots=snapshots)), snapshots=snapshots)

    def comp_min_distance_between_selections(self, sel1, sel2, snapshots='all'):
        """wrap function for universe"""

        def min_dist_from_universe(mda_universe):
            group_a = mda_universe.select_atoms(sel1)
            group_b = mda_universe.select_atoms(sel2)

            dist_a_b = cdist(group_a.positions, group_b.positions)

            min_dist = np.min(dist_a_b)
            return min_dist

        return self.generator_to_series(self.temporarily_provide_xyz(
            self.stream_trajectory_iterator(min_dist_from_universe, snapshots=snapshots)), snapshots=snapshots)


class MdaJob(object):

    def __init__(self, job):
        self.job = job
        self.MdaRuns = [MdaRun(x) for x in self.job]

        # add relevant methods to this instance
        run_methods = self.available_methods()
        for method in run_methods:
            method_name = method[0]
            if method_name[:5] == 'comp_':    # filter for relevant methods
                setattr(self, method_name, partial(self.func_on_all_runs, method_name))

    def func_on_all_runs(self, func, *args, **kwargs):
        """:return pandas.DataFrame with results of func for all runs"""
        all_results = []

        for sel_MdaRun in tqdm.tqdm_notebook(self.MdaRuns):
            all_results.append(getattr(sel_MdaRun, func)(*args, **kwargs))

        dataframe = pd.DataFrame(all_results).transpose()
        dataframe.columns.rename('runs', inplace=True)
        return dataframe

    def available_methods(self):
        """:return a list of all methods for the first MdaRun instance. These can be used with 'func_on_all_runs'"""
        return inspect.getmembers(self.MdaRuns[0], predicate=inspect.ismethod)
    
    
class MdaProject(object):
    
    def __init__(self, project):
        self.project = project
        self.MdaJobs = [MdaJob(x) for x in self.project.jobs]

        # add relevant methods to this instance
        run_methods = self.available_methods()
        for method in run_methods:
            method_name = method[0]
            if method_name[:5] == 'comp_':    # filter for relevant methods
                setattr(self, method_name, partial(self.func_on_all_jobs, method_name))

    def func_on_all_jobs(self, func, *args, **kwargs):
        """:return xarray.DataArray with results of func for all jobs"""
        all_results = []

        for sel_MdaJob in tqdm.tqdm_notebook(self.MdaJobs):
            all_results.append(getattr(sel_MdaJob, func)(*args, **kwargs))

        dataarray_list = [df.to_xarray().to_array() for df in all_results]
        dataarray = xr.concat(dataarray_list, dim="jobs")
        # Unfortunately, 'runs' are called 'variable' Although columns are named properly in DataFrame. Rename it!
        dataarray = dataarray.rename({'variable': 'runs'})
        # Add coordinates for readability
        dataarray.coords["jobs"] = range(len(self.MdaJobs))
        return dataarray

    def available_methods(self):
        """:return a list of all partials for the first MdaJob instance. Note: they were imported with partial in jobs,
         therefore we cannot look for methods at the project level. The partials can be used with 'func_on_all_runs'"""
        return inspect.getmembers(self.MdaJobs[0], predicate=self.ispartial)

    @staticmethod
    def ispartial(attr):
        return type(attr) == partial
