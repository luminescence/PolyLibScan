import os as _os
import numpy as np
import jinja2 as ji
import PolyLibScan.helpers.numpy as np_help
import itertools as it

class PymolPose(object):
    '''Parent class of pymol Visualisation
    '''
    atom_names = ['C', 'N', 'O', 'S', 'H']
    def __init__(self, pymol):
        self.pymol = pymol
        self._module_folder = _os.path.dirname(__file__)
        self.pymol_handle = self.pymol.pymol_handle


    def _poly_poses(self, sim, state='end'):
        '''Create generator of data for polymer-pdb via the pymol_polymer.tpl
        Each yield (run) the coordinates of the polymer are updated, while the constant
        information is left unchanged
        '''
        if state not in ['start', 'end']:
            raise AttributeError("state must have value of 'start' or 'end'.")
        mask = np.in1d(sim._parse.load_traj_type_order(), sim.particle_ids['polymer'])
        pose_data = self._create_pose_array(sim, mask)
        for run in sim:
            np_help.copy_fields(pose_data, self.traj_data(run, state, mask), ['x','y', 'z'])
            yield pose_data

    def _protein_poses(self, sim, state='end'):
        '''Create generator of data for polymer-pdb via the pymol_polymer.tpl
        Each yield (run) the coordinates of the polymer are updated, while the constant
        information is left unchanged
        '''
        if state not in ['start', 'end']:
            raise AttributeError("state must have value of 'start' or 'end'.")
        mask = np.in1d(sim._parse.load_traj_type_order(), sim.particle_ids['protein'])
        pose_data = self._create_pose_array(sim, mask)
        for run in sim:
            np_help.copy_fields(pose_data, self.traj_data(run, state, mask), ['x','y', 'z'])
            yield pose_data


    def _create_pose_array(self, sim, mask):
        '''
        '''
        pose_list_type = [('id', '>i2'), ('ele', '|S1'), ('res_name', '|S3'), 
                          ('x', np.float), ('y', np.float), ('z', np.float)]
        elements = dict(zip(sim.particle_ids['polymer'], self.atom_names))
        pose_data = np.zeros(mask.sum(), dtype=pose_list_type)
        pose_data['id'] = np.arange(1, mask.sum()+1)
        pose_data['ele'] = map(lambda x:elements[x], sim[0].coordinates()['end'][mask]['atom_type'])
        pose_data['res_name'] = sim.sequence['monomer']
        return pose_data

    def traj_data(self, run, state, mask):
        return run.coordinates()[state][mask]

    def add_polymers(self, state='end'):
        polymer_pdb_path, name = self._create_polymer_pdb(state=state)
        self.pymol_handle.load(polymer_pdb_path)
        self.pymol_handle.show_as('sphere', name)
        self.pymol_handle.do('cmd.set("sphere_scale", 2.0, "all")')

    def add_protein(self, state='end'):
        protein_pdb_path, name = self._create_protein_pdb(state=state)
        self.pymol_handle.load(protein_pdb_path)
        self.pymol_handle.show_as('sphere', name)
        self.pymol_handle.do('cmd.set("sphere_scale", 2.0, "all")')


class Project(PymolPose):
    """docstring for PymolVisPolyType"""
    def __init__(self, pymol):
        super(Project, self).__init__(pymol)


class Type(PymolPose):
    """docstring for PymolVisPolyType"""
    def __init__(self, pymol):
        super(Type, self).__init__(pymol)

    def _create_polymer_pdb(self, state='end'):
        if state in ['start', 'end']:
            gen = it.chain.from_iterable((self._poly_poses(sim, state=state) for sim in self.sims))
            file_name = '%s_%s_poses.pdb' % (self.poly_type.name, state)
        else:
            raise AttributeError("state supports only 'start' or 'end' as values.")
        output_path = self.type_folder.joinpath(file_name).as_posix()
        pymol_name = file_name[:-4] # remove '.pdb'
        with open('%s/pymol_polymer.tpl' % self._module_folder) as f:
            template = ji.Template(f.read())
        with open(output_path, 'w') as f2:
            f2.write(template.render(pymol=gen))
        return output_path, pymol_name

    def _create_protein_pdb(self, state='end'):
        if state in ['start', 'end']:
            gen = it.chain.from_iterable((self._poly_poses(sim, state=state) for sim in self.sims))
            file_name = '%s_%s_protein_poses.pdb' % (self.poly_type.name, state)
        else:
            raise AttributeError("state supports only 'start' or 'end' as values.")
        output_path = self.type_folder.joinpath(file_name).as_posix()
        pymol_name = file_name[:-4] # remove '.pdb'
        with open('%s/pymol_polymer.tpl' % self._module_folder) as f:
            template = ji.Template(f.read())
        with open(output_path, 'w') as f2:
            f2.write(template.render(pymol=gen))
        return output_path, pymol_name

class Job(PymolPose):
    """docstring for PymolVisJob"""
    def __init__(self, pymol):
        super(Job, self).__init__(pymol)
        pass

    def _create_polymer_pdb(self, state='end'):
        if state in ['start', 'end']:
            gen = self._poly_poses(self.sim, state=state)
            file_name = '%s_%s_poses.pdb' % (self.sim.meta['id'], state)
        else:
            raise AttributeError("state supports only 'start' or 'end' as values.")
        output_path = self.db_folder.joinpath(file_name).as_posix()
        pymol_name = file_name[:-4]
        with open('%s/pymol_polymer.tpl' % self._module_folder) as f:
            template = ji.Template(f.read())
        with open(output_path, 'w') as f2:
            f2.write(template.render(pymol=gen))
        return output_path, pymol_name

class Run(PymolPose):
    """docstring for Run"""
    def __init__(self, pymol):
        super(Run, self).__init__(pymol)
        self.sim = self.pymol.sim
        self.run = self.pymol.run
        
    def _create_polymer_pdb(self, state='end'):
        if state in ['start', 'end', 'full']:
            gen = self._poly_poses(self.sim, state=state)
        else:
            raise AttributeError("state supports only 'start' or 'end' as values.")
        file_name = '%s__run_%d-%s_poses.pdb' % (self.sim.meta['id'], self.run.Id, state)
        output_path = self.pymol.db_folder.joinpath(file_name).as_posix()
        pymol_name = file_name[:-4]
        with open('%s/pymol_polymer.tpl' % self._module_folder) as f:
            template = ji.Template(f.read())
        with open(output_path, 'w') as f2:
            f2.write(template.render(pymol=gen))
        return output_path, pymol_name
    
    def _poly_poses(self, sim, state='end'):
        '''Create generator of data for polymer-pdb via the pymol_polymer.tpl
        Each yield (run) the coordinates of the polymer are updated, while the constant
        information is left unchanged
        '''
        mask = np.in1d(sim._parse.load_traj_type_order(), sim.particle_ids['polymer'])
        pose_data = self._create_pose_array(sim, mask)
        if state in ['start', 'end']:
            for run in [sim[self.run.Id]]:
                np_help.copy_fields(pose_data, self.traj_data(run, state, mask))
                yield pose_data
        if state == 'full':
            for step_data in self.run.polymer_trajectory():
                pose_data['x'] = step_data['xyz'][:,0]
                pose_data['y'] = step_data['xyz'][:,1]
                pose_data['z'] = step_data['xyz'][:,2]
                yield pose_data
        else:
            raise AttributeError("state must have value of 'start', 'end' or 'full'.")

    def traj_data(self, run, state, mask):
            return run.coordinates()[state][mask]
        