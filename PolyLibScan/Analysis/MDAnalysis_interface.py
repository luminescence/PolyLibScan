from functools import partial
import inspect
import numpy as np
import os

import MDAnalysis as mda
import pandas as pd
from scipy.spatial.distance import cdist
import xarray as xr

from PolyLibScan.helpers.jupyter_compatibility import agnostic_tqdm
from PolyLibScan.Analysis.sim_run import AtomFilter


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
        masses = self.job.get_particle_parameters(parameter='mass')

        atom_group = mda_universe.select_atoms('all')
        atom_group.masses = masses

    def universe_creator(self, snapshot):
        no_particles = len(self.job.trajectory_order)
        mda_universe = mda.Universe(self.initial_xyz_file, snapshot['xyz'].reshape(1, no_particles, 3))
        self.assign_masses_in_universe(mda_universe)
        return mda_universe

    def stream_trajectory_iterator(self, func):
        for ts, x in agnostic_tqdm(enumerate(self.run.trajectory())):
            if ts in self.snapshots_range:
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

    def generator_to_series(self, input_generator):
        """store everything from generator (in our case trajectory) in pd.Series"""
        output_list = []
        for x in input_generator:
            output_list.append(x)

        series = pd.Series(output_list)
        snapshots_as_indeces = sorted(self.snapshots_range)
        series.index = snapshots_as_indeces

        series.index.rename('snapshots', inplace=True)

        return series

    def return_universe(self, snapshots='all'):
        """return universe for development purposes"""

        self.snapshots_range = self.determine_snapshots_range(snapshots)

        def universe(mda_universe):
            return mda_universe

        return self.generator_to_series(self.temporarily_provide_xyz(
            self.stream_trajectory_iterator(universe)))

    def molecule_selection_string(self, selection='all'):
        """return selection string for MDAnalysis"""

        if selection in ('all', 'full'):
            return 'all'
        elif selection in ('polymer', 'protein'):
            molecule = selection
            type_filter = AtomFilter(self.job.trajectory_order, self.job.sequence, self.job.particle_ids[molecule], molecule=molecule)

            particle_ids = np.where(type_filter.mask)
            particle_ids = particle_ids[0] + 1 #    MDAnalysis counts from 1 not 0

            id_strings = ['bynum ' + str(x) for x in particle_ids]

            selection_string = ' or '.join(id_strings)

            return selection_string
        else:
            raise ValueError('Cannot interpret selection %s' % selection)

    # scientifically meaningful methods
    # start with comp_

    def comp_radius_of_gyration(self, selection='all', snapshots='all'):
        """wrap function for universe"""

        self.snapshots_range = self.determine_snapshots_range(snapshots)

        def rg_from_universe(mda_universe):
            return mda_universe.select_atoms(selection).radius_of_gyration()

        return self.generator_to_series(self.temporarily_provide_xyz(
            self.stream_trajectory_iterator(rg_from_universe)))

    def comp_min_distance_between_selections(self, sel1, sel2, snapshots='all'):
        """wrap function for universe"""

        self.snapshots_range = self.determine_snapshots_range(snapshots)

        def min_dist_from_universe(mda_universe):
            group_a = mda_universe.select_atoms(sel1)
            group_b = mda_universe.select_atoms(sel2)

            dist_a_b = cdist(group_a.positions, group_b.positions)

            min_dist = np.min(dist_a_b)
            return min_dist

        return self.generator_to_series(self.temporarily_provide_xyz(
            self.stream_trajectory_iterator(min_dist_from_universe)))


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

        for sel_MdaRun in agnostic_tqdm(self.MdaRuns):
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

        for sel_MdaJob in agnostic_tqdm(self.MdaJobs):
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
