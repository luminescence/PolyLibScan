import tqdm
import pathlib2 as pl
import pandas as pd
import collections as col
import yaml
import warnings
import plotting
import bayesModels as bayes
import PolyLibScan.helpers.idGenerator as idGen
import pymol_visualisation as pym
import visualize
import poly_type
import numerics as num_
import job as job_class

warnings.filterwarnings("ignore")

class Project(plotting.Project, bayes.Project):
    '''

    with_pymol: Loads pymol visualisation option

    Pymol:
    you have to start pymol with the '-R' flag for this to function.
    '''
    def __init__(self, 
                 project_path, 
                 experimental_data=None, 
                 parameters=None, 
                 protein_path=None, 
                 with_pymol=True):
        self.path = pl.Path(project_path)
        self._id_gen = idGen.IdGen()
        self._endstate_matrix = None
        self.parameters = self._init_parameters(parameters)
        self._visualize = visualize.Visualize(self)
        # stores jobs in list (jobs)
        # and dict form (polymer_types); polymer_types stores jobs sorted by polymer type
        self.jobs, self.polymer_types = self.read_jobs(self.path.joinpath('jobs'), with_pymol=with_pymol)
        self.protein_path = self._init_protein(protein_path)
        if with_pymol:
            self.pymol = pym.PymolVisProject(self, self.protein_path)

        # Experimental Data needs information from job objects
        self.experimental_data = self._init_experimental_data(experimental_data)
        if self.experimental_data is not None:
            for p_type in self.polymer_types.values():
                if p_type.name in self.experimental_data:
                    p_type.ic50 = self.experimental_data[p_type.name]
                else:
                    p_type.ic50 = -1

    def search_static(self, file_name):
        file_list = list(self.path.joinpath('static').glob(file_name))
        if len(file_list) == 1:
            return file_list[0]
        else:
            raise IOError('No file found with name: %s' % file_name)

    def _init_parameters(self, file_path):
        if file_path:
            return file_path
        else:
            try:
                default_file_name = 'parameters.yml'
                return self.search_static(default_file_name).as_posix()
            except IOError:
                print 'Parameter file not specified and no file "%s" found in static folder.' % default_file_name
                return None

    def _init_protein(self, file_path):
        if file_path:
            return file_path
        else:
            try:
                default_file_name = '*.pdb'
                return self.search_static(default_file_name).as_posix()
            except IOError:
                print 'pdb file not specified and no file "%s" found in static folder.' % default_file_name
                return None
    
    def _init_experimental_data(self, file_path):
        if file_path:
            return file_path
        else:
            try:
                default_file_name = 'ic50.h5'
                return self.search_static(default_file_name).as_posix()
            except IOError:
                print 'Inhibition file not specified and no file "%s" found in static folder.' % default_file_name 
                return None

    def read_jobs(self, path, with_pymol=True):
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
                job = job_class.Job(self, folder.joinpath('jobdata.h5').as_posix(), with_pymol=with_pymol)
                poly_id = self._id_gen[job.meta['poly_name']]
                job.Id = poly_id
                # storing job in the two containers
                jobs.append(job)
                sims_by_type[job.meta['poly_name']].append(job)
        polymer_types = {name: poly_type.PolymerTypeSims(self, sims, with_pymol=with_pymol) 
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
            else:
                self._experimental_data = None
                return    
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

    def _scatter_data(self, subset=None, with_errors=False, with_labels=False, with_crossvalidation=False, 
                            property_=None, confidence_interval=0.95, min_dist_to_ac=10, 
                            ignore_experiment=False):
        if subset:
            polymer_list = subset
            if self.experimental_data is not None:
                experimental = self.experimental_data[polymer_list]    
        elif (self.experimental_data is not None) and (not ignore_experiment):
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
            df2 = num_.energy_with_error(energy_matrix, confidence_interval)
            results = pd.concat([df1, df2], axis=1)
        else:
            results = pd.DataFrame(index=distance_matrix.columns)
            results['energy_mean'] = energy_matrix.mean()
            results['dist_mean'] = distance_matrix.count()/float(full_distance_matrix.shape[0])
        # drop polymer types that are not present in the min-radius
        results.dropna(inplace=True)
        if (self.experimental_data is not None) and (not ignore_experiment):
            results['color_by_inhibition'] = experimental[results.index].apply(lambda x:'b' if x>0 else 'r')
        else:
            results['color_by_inhibition'] = 'k'

        if with_crossvalidation:
            input_data = results.loc[:, ('energy_mean', 'dist_mean')]
            classification = results['color_by_inhibition'].apply(lambda x:x=='b')
            true_predictions, model_predictions, probabilities = self._cross_validate(input_data, classification)
            results['true_predictions'] = true_predictions
            results['model_predictions'] = model_predictions
            results['probabilities'] = probabilities
        if not property_:
            results['property'] = 'k'
        elif property_ == 'charge':
            results['property'] = pd.Series(data={name: pt.charge_average().sequence 
                                                for name,pt in self.polymer_types.iteritems()})
        elif property_ == 'hydrophobicity':
            results['property'] = pd.Series(data={name: pt.hydrophobic_average().sequence 
                                                for name,pt in self.polymer_types.iteritems()})
        else:
            raise AttributeError("only 'charge' and 'hydrophobicity' supported as properties.")

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