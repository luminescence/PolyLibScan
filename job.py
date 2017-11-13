import getpass
import itertools as it
import os

import numpy as np
import yaml
from lammps import lammps

import Save
import Tools
import helpers.git as _git
import helpers.time as tm
from Tools.lmp_control import lmp_controler

__git_hash__ = _git.get_git_hash(__file__)

class Job(lmp_controler):

    def __init__(self, config_path):
        super(Job, self).__init__()
        self.config = self.read_config(config_path)
        self.config_with_setup = os.path.join(os.path.dirname(config_path), 'config_with_setup.yml')
        if os.path.exists(self.config_with_setup):
            with open(self.config_with_setup) as f:
                setup_config = yaml.load(f)
            self.config.sim_parameter['poly_sequence'] = setup_config['sim_parameter']['poly_sequence']

        self.username = getpass.getuser()

        self.config.sim_path['sim_list'] = os.path.join(self.config.sim_path['root'], 'sim.list')
        open(self.config.sim_path['sim_list'], 'a').close()
        if self.config.sim_parameter['local'] == 1:
            self._switch_to_local()
        else:
            new_paths = self.create_folders(self.config.lmp_path['root'])
            self.config.lmp_path.update(new_paths)

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

        self.protein = self.protein_creator.create()
        # Polymer
        if 'poly_sequence' in self.config.sim_parameter:
            self.polymer_creator = Tools.PolymerCreator(self.env, self.config.sim_parameter['poly_sequence'], mode='cycle')
        else:
            self.polymer_creator = Tools.PolymerCreator(self.env, 
                self.config.sim_parameter['monomers'], weights=self.config.sim_parameter['weights'], 
                length=self.config.sim_parameter['poly_length'])
        self.poly = self.polymer_creator.create()

        self.sim = Tools.EnvManipulator(self.env, auto_repulsion=False)
        self.sim.create_random_start_positions()
        self.setup_writer = Tools.LmpWriter(self.env)

        # Update Lammps Parameters
        self.config.sim_parameter['poly_sequence'] = [monomers.type_.name for monomers in self.poly.data['monomers']]
        self.config.sim_parameter['named_sequence'] = [particle.type_.name for particle in self.poly.data['particles']]
        self.config.sim_parameter['id_sequence'] = [particle.type_.Id for particle in self.poly.data['particles']]
        self.config.sim_parameter['active_site'] = self.set_active_site()
        self.config.lmp_parameter['active_site_ids'] = self.config.sim_parameter['active_site']['xyz']
        self.config.lmp_parameter['monomer_ids'] = self._get_monomer_ids()
        self.config.lmp_parameter['bb_id'] = self.env.atom_type['BB_bb'].Id
        self.config.lmp_parameter['ghost_id'] = self.env.atom_type['BB_ghost_bb'].Id
        # the config with added information is always saved to the root directory
        self.config.save(self.config_with_setup)
        # save particle list
        p_list_path = os.path.join(self.config.sim_path['root'], 'particle_list.npy')
        self.save_particle_list(p_list_path)
        # fifos can only be created if the monomer ids are known
        self.fifo = self._create_fifos()

    def save_particle_list(self, path):
        dtype = [('xyz', np.int), ('p_id', np.int), ('name', 'S6'), 
                 ('chain', 'S1'), ('atom', 'S6'), ('res_id', np.int), ('iCode', 'S1')]

        molecules = sorted(self.env.molecules.values(), key=lambda x: x.Id)
        particle_gen = it.chain(*[molecule.data['particles'] for molecule in molecules])
        particle_count = sum([len(mol.data['particles']) for mol in self.env.molecules.values()])
        particle_dat = np.empty(particle_count, dtype)
        for i,p in enumerate(particle_gen):
            particle_dat[i] = (p.Id, p.type_.Id, p.residue.name, p.residue.chain, 
                                p.residue.id[0], p.residue.id[1], p.residue.id[2]) 
        particle_dat.tofile(path)

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
        start_idx = self._get_last_uncompleted_index()
        end_idx = self.config.sim_parameter['sampling_rate']
        # check if simulations is already completed
        if start_idx == -1:
            return 
        for i in xrange(start_idx,  end_idx):
            self.generate_new_sim(i)
            # start next LAMMPS run
            self.lmps_run(i, self.config.lmp_parameter, self.config.lmp_path, fifos=self.fifo)
            # report completed simulation so restarting jobs will know
            # also, it notes the machine and folder, so scattered info can be retrieved
            self._mark_complete(i)
        if start_idx != -1:
            self._mark_complete(-1)

    def create_local_env(self, local_dir='/data/'):
        '''Create unique local job-folder and create the 
        '''
        name_comp = self.config.sim_path['root'].split('/')[-3:]
        if name_comp[1] == 'jobs':
            del name_comp[1]
        else:
            del name_comp[0]
        folder_name = '-'.join(name_comp)
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

    def _get_last_uncompleted_index(self):
        '''get the last line of the sim_list file 
        and return the index.
        '''
        with open(self.config.sim_path['sim_list']) as f:
            completed_sims = f.read().split('\n')
        if len(completed_sims) > 1:
            last_line = completed_sims[-2]
        else:
            return 0
        return int(last_line[:5])

    def _switch_to_local(self):
        '''if the data of the simulations 
        is to be stored locally, the lmp paths 
        are changed to local folders.
        '''
        new_paths = self.create_local_env()
        self.config.lmp_path.update(new_paths)

    def _mark_complete(self, index):
        if self.config.sim_parameter['local'] == 1:
            path = self.config.lmp_path['local_root']
        else:
            path = self.config.lmp_path['root']
        with open(self.config.sim_path['sim_list'], 'a') as f:
            info = '%05d;%s;%s;%s\n' % (index, localhost(), path, tm.time_string())
            f.write(info)

    def _get_monomer_ids(self):
        monomer_ids = set([p.type_.Id for p in self.poly.data['particles']])
        return sorted(monomer_ids)      


def localhost():
    """return the nodename of the computer.
    """
    return os.uname()[1]
