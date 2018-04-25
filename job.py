import getpass
import itertools as it
import os
import collections as col
import numpy as np
import warnings
import yaml
from lammps import lammps

import Save
import Tools
import helpers.git as _git
import helpers.sim_list
from Tools import lmp_control

__git_hash__ = _git.get_git_hash(__file__)

class Job(object):

    def __init__(self, config_path):
        super(Job, self).__init__()
        self.config = self.read_config(config_path)
        self.config_with_setup = os.path.join(os.path.dirname(config_path), 'config_with_setup.yml')
        if os.path.exists(self.config_with_setup):
            with open(self.config_with_setup) as f:
                setup_config = yaml.load(f)
            self.config.sim_parameter['poly_sequence'] = setup_config['sim_parameter']['poly_sequence']
        self.username = getpass.getuser()
        if self.config.sim_parameter['local'] == 1:
            self._switch_to_local()
        else:
            new_paths = self.create_folders(self.config.lmp_path['root'])
            self.config.lmp_path.update(new_paths)
        # create_sim_list depends on the results of the previous condition
        self.sim_list = self.create_sim_list()
        with open(self.config.sim_path['config']) as f:
            self.physical_model = yaml.load(f)

    def read_config(self, path):
        return Tools.config.JobConfig(path)

    def _create_fifos(self):
        fifo = {}
        for name, data in self.config.fifo.items():
            data['path'] = os.path.join(self.config.lmp_path['fifo'], data['path'])
            fifo[name] = Tools.FiFo.create_from_dict(name, self, data)
        return fifo

    def set_active_site(self):
        as_data = self.sim.activeSiteParticles(self.protein, self.config.sim_path['active_site'])
        # the mapping ensures that vanilla python objects are used as the initial data is numpy-based.
        return  {'xyz': map(int, as_data['xyz']), 
                'chain': map(str, as_data['chain']), 
                'pdb_id': map(int, as_data['pdb_id']), 
                'iCode': map(str, as_data['iCode'])}

    def setup_env(self):
        self.env = Tools.Environment(self.config.sim_path['config'])
        # Protein
        self.protein_creator = Tools.ProteinCreator(
            self.env,
            self.config.sim_path['protein'],
            with_ions=True,
            surface_file=self.config.sim_path['surface_db'],
            protonation_file=self.config.sim_path['protonation_db'],
            ph=self.config.sim_parameter['pH']
        )

        number_of_proteins = self.config.sim_parameter['stoichiometry'][0]
        number_of_polymers = self.config.sim_parameter['stoichiometry'][1]

        # protein
        if number_of_proteins == 0:
            pass
        elif number_of_proteins == 1:
            self.protein = self.protein_creator.create()
        else:
            raise NotImplementedError(
                'Multiple protein molecules are currently not supported! \n You will need to edit: Tools.lmp_combine.create_random_start_positions')

        # polymer
        if 'poly_sequence' in self.config.sim_parameter:
            self.polymer_creator = Tools.PolymerCreator(self.env, self.config.sim_parameter['poly_sequence'],
                                                        mode='cycle')
        else:
            self.polymer_creator = Tools.PolymerCreator(self.env,
                                                        self.config.sim_parameter['monomers'],
                                                        weights=self.config.sim_parameter['weights'],
                                                        length=self.config.sim_parameter['poly_length'])
        if number_of_polymers == 1:
            self.poly = self.polymer_creator.create()

        self.sim = Tools.EnvManipulator(self.env, auto_repulsion=False)
        self.sim.create_random_start_positions()
        self.setup_writer = Tools.LmpWriter(self.env)

        # save particle list
        p_list_path = os.path.join(self.config.sim_path['root'], 'particle_list.npy')
        self.save_particle_list(p_list_path)

        self.lmp_settings = {}
        self.lmp_settings['stoichiometry'] = self.config.sim_parameter['stoichiometry']
        if number_of_polymers == 1:
            self.config.sim_parameter['poly_sequence'] = [monomers.type_.name for monomers in self.poly.data['monomers']]
            self.config.sim_parameter['named_sequence'] = [particle.type_.name for particle in self.poly.data['particles']]
            self.config.sim_parameter['id_sequence'] = [particle.type_.Id for particle in self.poly.data['particles']]
            self.lmp_settings['monomer_ids'] = self._get_monomer_ids()
            self.lmp_settings['poly_sequence'] = self.config.sim_parameter['id_sequence']
        elif number_of_polymers > 1:
            raise NotImplementedError('Multiple polymer molecules are currently not supported!')
        if number_of_proteins == 1:
            self.config.sim_parameter['active_site'] = self.set_active_site()
            self.lmp_settings['active_site_ids'] = self.config.sim_parameter['active_site']['xyz']
            self.lmp_settings['bb_id'] = self.env.atom_type['BB_bb'].Id
            self.lmp_settings['ghost_id'] = self.env.atom_type['BB_ghost_bb'].Id
        elif number_of_proteins > 1:
            raise NotImplementedError('Multiple protein molecules are currently not supported!')
        self.lmp_settings['particle_ids'] = self.particle_list['p_id'].tolist()

        if bool(self.config.lmp_parameter):
            warnings.warn('The following parameters in the physical model will be overwritten: %s' %
                          self.config.lmp_parameter)
            self.physical_model['MD_parameters'].update(self.config.lmp_parameter)

        # the config with added information is always saved to the root directory
        self.config.save(self.config_with_setup)
        # fifos can only be created if the monomer ids are known
        self.fifo = self._create_fifos()

    def save_particle_list(self, path):
        dtype = [('xyz', np.int), ('p_id', np.int), ('name', 'S10'), 
                 ('chain', 'S1'), ('atom', 'S6'), ('res_id', np.int), 
                 ('iCode', 'S1'), ('charge', np.float)]

        molecules = sorted(self.env.molecules.values(), key=lambda x: x.Id)
        particle_gen = it.chain(*[molecule.data['particles'] for molecule in molecules])
        particle_count = sum([len(mol.data['particles']) for mol in self.env.molecules.values()])
        particle_dat = np.empty(particle_count, dtype)
        for i,p in enumerate(particle_gen):
            particle_dat[i] = (p.Id, p.type_.Id, p.residue.name , p.residue.chain, 
                               p.residue.id[0] , p.residue.id[1], p.residue.id[2],
                               p.charge)
        self.particle_list = particle_dat
        np.save(path, particle_dat)

    def terminate_fifos(self):
        for name, fifo in self.fifo.items():
            fifo.terminate()

    def setup_job_save(self):
        self.compactor = Save.JobSave(self.config.lmp_path)

    def save(self):
        '''Save the data of the completed simulations to HDF5 database.
        '''
        if not hasattr(self, 'compactor'):
            self.setup_job_save()
        try:
            self.compactor.save()
            self.compactor.save_versions(lmp_version=lammps().version())
        finally:
            self.compactor.db.close()

    def clean_up(self):
        try:
            self.terminate_fifos()
        except OSError:
            pass
        self.compactor.clean_up()

    def generate_new_sim(self, index):
        # Create new Start Conditions
        self.sim.create_random_start_positions()
        # Create New Setup File
        self.setup_writer.write('%s/%05d' % (self.config.lmp_path['input'], index))

    def run(self):
        for i in self.sim_list:
            self.generate_new_sim(i)
            # start next LAMMPS run
            lmp_controller = lmp_control.LmpController(i, self.lmp_settings, self.config.lmp_path, self.physical_model, self.fifo)
            lmp_controller.lmps_run()
            # report completed simulation so restarting jobs will know
            # also, it notes the machine and folder, so scattered info can be retrieved
            self.sim_list.mark_complete(i)

    def create_local_env(self, local_dir='/data/'):
        '''Create local job-folder and needed subfolders.
        The job-folder is placed in the folder under the local_dir.
        The local dir is either composed via the 'local_root' path in 
        config.sim_path if it exists or the functions 'local_dir' argument.
        Note that the config takes precedent!
        '''

        # get the folder of the the project and job name
        name_comp = self.config.sim_path['root'].split('/')[-3:]
        # depending on the existence of the 'jobs' folder 
        # choose project-root and job-root
        if name_comp[1] == 'jobs':
            project_root = name_comp[0]
        else:
            project_root = name_comp[1]
        job_root = name_comp[2]
        folder_name = '%s-%s' % (project_root, job_root)
        if 'local_root' in self.config.sim_path and self.config.sim_path['local_root'] != '':
            local_userdir = self.config.sim_path['local_root']
        else:
            local_userdir = os.path.join(local_dir, self.username)
        local_folder = os.path.join(local_userdir, folder_name)
        if not os.path.exists(local_folder):
            os.mkdir(local_folder)
        new_paths = self.create_folders(local_folder)
        return new_paths

    def create_folders(self, root_folder):
        '''Create the three needed folder in the root_folder
        in case they do not already exist.
        '''
        new_paths = {'local_root': root_folder}
        for sub_folder in ['input', 'output', 'fifo']:
            folder = os.path.join(root_folder, sub_folder)
            if not os.path.exists(folder):
                os.mkdir(folder)
            new_paths[sub_folder] = folder
        return new_paths

    def _switch_to_local(self):
        '''if the data of the simulations 
        is to be stored locally, the lmp paths 
        are changed to local folders.
        '''
        new_paths = self.create_local_env()
        self.config.lmp_path.update(new_paths)

    def create_sim_list(self):
        list_path = os.path.join(self.config.sim_path['root'], 'sim.list')
        if self.config.sim_parameter['local'] == 1: 
            data_folder = self.config.lmp_path['local_root'] 
        else: 
            data_folder = self.config.lmp_path['root']
        return helpers.sim_list.SimList(list_path,
                                        data_folder,
                                        self.config.sim_parameter['sampling_rate'])

    def _get_monomer_ids(self):
        monomer_ids = set([p.type_.Id for p in self.poly.data['particles']])
        return sorted(monomer_ids)      
