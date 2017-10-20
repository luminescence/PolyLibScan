import os as _os
import numpy as np
import jinja2 as ji
import PolyLibScan.helpers.numpy_helpers as np_help
import itertools as it

def _get_pdb_template():
    with open('%s/pymol_polymer.tpl' % _os.path.dirname(__file__)) as f:
        template = ji.Template(f.read())
    return template

class PymolPose(object):
    '''Parent class of pymol Visualisation
    '''
    atom_names = 10*['C', 'N', 'O', 'S', 'H']
    
    template = _get_pdb_template()

    def __init__(self, pymol):
        self.pymol = pymol
        self.pymol_handle = self.pymol.pymol_handle


    def _poses(self, sim, molecule='polymer', state='end'):
        '''Create generator of data for polymer-pdb via the pymol_polymer.tpl
        Each yield (run) the coordinates of the polymer are updated, while the constant
        information is left unchanged
        '''
        if state not in ['start', 'end']:
            raise AttributeError("state must have value of 'start' or 'end'.")
        mask = np.in1d(sim.trajectory_order, sim.particle_ids[molecule])
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
        polymer_pdb_path, name = self._create_pdb(molecule='polymer', state=state)
        self.pymol_handle.load(polymer_pdb_path)
        self.pymol_handle.show_as('sphere', name)
        self.pymol_handle.do('cmd.set("sphere_scale", 2.0, "all")')

    def add_protein(self, state='end'):
        protein_pdb_path, name = self._create_pdb(molecule='protein', state=state)
        self.pymol_handle.load(protein_pdb_path)
        self.pymol_handle.show_as('sphere', name)
        self.pymol_handle.do('cmd.set("sphere_scale", 2.0, "all")')

    def _create_pdb(self, molecule='polymer', state='end'):
        gen = self.pdb_data(molecule=molecule, state=state)
        file_name = self.file_name(molecule, state)
        output_path = self._save_pdb(file_name, gen)
        pymol_name = file_name[:-4] # remove '.pdb'
        return output_path, pymol_name

    def _save_pdb(self, file_name, pdb_data):
        output_path = self.output_folder.joinpath(file_name).as_posix()
        with open(output_path, 'w') as f:
            f.write(self.template.render(pymol=pdb_data))
        return output_path


class Type(PymolPose):
    """docstring for PymolVisPolyType"""
    def __init__(self, pymol):
        super(Type, self).__init__(pymol)
        self.sims = self.pymol.sims
        self.output_folder = self.pymol.type_folder

    def pdb_data(self, molecule='polymer', state='end'):
        return it.chain.from_iterable((self._poses(
                                            sim, 
                                            molecule=molecule, 
                                            state=state) 
                                       for sim in self.sims))

    def file_name(self, molecule, state):
        return '%s_%s_%s_poses.pdb' % (
                    self.sim.poly_type.name, 
                    molecule, 
                    state)


class Job(PymolPose):
    """docstring for PymolVisJob"""
    def __init__(self, pymol):
        super(Job, self).__init__(pymol)
        self.sim = self.pymol.sim
        self.output_folder = self.pymol.db_folder

    def pdb_data(self, molecule='polymer', state='end'):
        return self._poses(self.sim, molecule=molecule, state=state)

    def file_name(self, molecule, state):
        return '%s_%d_%s_%s_poses.pdb' % (
                    self.sim.poly_type.name, 
                    self.sim.Id, 
                    molecule, 
                    state)


class Run(PymolPose):
    """docstring for Run"""
    def __init__(self, pymol):
        super(Run, self).__init__(pymol)
        self.sim = self.pymol.sim
        self.run = self.pymol.run
        self.output_folder = self.pymol.db_folder
    
    def pdb_data(self, molecule='polymer', state='full'):
        return self._poses(self.sim, molecule=molecule, state=state)

    def file_name(self, molecule, state):
        return '%s_%d_%d_%s_%s_poses.pdb' % (
                    self.sim.poly_type.name, 
                    self.sim.Id, 
                    self.run.Id,
                    molecule, 
                    state)
    
    def _poses(self, sim, molecule='polymer', state='full'):
        '''Create generator of data for polymer-pdb via the pymol_polymer.tpl
        Each yield (run) the coordinates of the polymer are updated, while the constant
        information is left unchanged
        '''
        mask = np.in1d(sim.trajectory_order, sim.particle_ids[molecule])
        pose_data = self._create_pose_array(sim, mask)
        if state in ['start', 'end']:
            for run in [sim[self.run.Id]]:
                np_help.copy_fields(pose_data, self.traj_data(run, state, mask), ['x','y', 'z'])
                yield pose_data
        elif state == 'full':
            for step_data in self.run.polymer_trajectory():
                pose_data['x'] = step_data['xyz'][:,0]
                pose_data['y'] = step_data['xyz'][:,1]
                pose_data['z'] = step_data['xyz'][:,2]
                yield pose_data
        else:
            raise AttributeError("state must have value of 'start', 'end' or 'full'.")
