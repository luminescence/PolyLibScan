import os as _os
import numpy as np
import jinja2 as ji
import PolyLibScan.helpers.numpy_helpers as np_help
import itertools as it

from PolyLibScan.Analysis.sim_run import AtomFilter

def _get_pdb_template(name):
    with open('%s/%s' % (_os.path.dirname(__file__),name)) as f:
        template = ji.Template(f.read())
    return template


class PymolPose(object):
    '''Parent class of pymol Visualisation
    '''
    atom_names = 200*['C', 'N', 'O', 'S', 'H']
    
    template = _get_pdb_template('pymol_polymer.tpl')

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
        type_filter = AtomFilter(sim.trajectory_order, sim.sequence, sim.particle_ids[molecule], molecule=molecule)
        pose_data = self._create_pose_array(sim, molecule, type_filter.mask)
        for run in sim:
            np_help.copy_fields(pose_data, self.traj_data(run, state, type_filter.mask), ['x', 'y', 'z'])
            yield pose_data

    def _create_pose_array(self, sim, molecule, mask):
        '''
        '''
        pose_list_type = [('id', '>i2'), ('ele', '|S1'), ('res_name', '|S3'), 
                          ('x', np.float), ('y', np.float), ('z', np.float)]
        elements = dict(zip(sim.particle_ids[molecule], self.atom_names))
        pose_data = np.zeros(mask.sum(), dtype=pose_list_type)
        pose_data['id'] = np.arange(1, mask.sum()+1)
        pose_data['ele'] = map(lambda x:elements[x], sim[0].coordinates()['end'][mask]['atom_type'])
        if molecule == 'polymer':
            seq = sim.sequence['monomer']
        else:
            seq_mask = np.in1d(sim.particle_list['p_id'], sim.particle_ids['protein'])
            seq = sim.particle_list['name'][seq_mask]
        pose_data['res_name'] = seq
        return pose_data

    def traj_data(self, run, state, mask):
        return run.coordinates()[state][mask]

    def add_CG_pose(self, molecule, state='end'):
        pdb_path, name = self._create_pdb(molecule=molecule, state=state)
        self.pymol_handle.load(pdb_path)
        self.pymol_handle.show_as('sphere', name)
        self.set_CG_radii(molecule, name)

    def set_CG_radii(self, molecule, selection_name):
        """assign the correct vdw radii to all visualized beads"""
        sim = self.pymol.sim
        radii = sim.get_particle_parameters(parameter='radius')
        type_filter = AtomFilter(sim.trajectory_order, sim.sequence, sim.particle_ids[molecule], molecule=molecule)
        # only modify vdw radii for beads of the correct molecule,
        # indeces would not fit for the second molecule (usually polymer) if this was not done
        molecule_indeces = np.where(type_filter.mask)[0]
        vdw_radii_for_beads = [radii[i] for i in molecule_indeces]

        # altering vdw radii individually for all beads takes very long; bundle all beads with the same vdw radius
        possible_vdw_radii = set(vdw_radii_for_beads)

        for radius in possible_vdw_radii:
            radius_indeces = np.where(np.array(vdw_radii_for_beads) == radius)[0]
            radius_index_strings = ['index ' + str(index + 1) for index in radius_indeces] # Pymol counts from 1 not 0
            radius_merged_indeces_string = ' or '.join(radius_index_strings)
            radius_selection = '%s and (%s)' % (selection_name, radius_merged_indeces_string)

            self.pymol_handle.do('alter %s, vdw=%s' % (radius_selection, radius))

        self.pymol_handle.do('rebuild')

    def add_polymers(self, state='end'):
        self.add_CG_pose(molecule='polymer', state=state)

    def add_protein(self, state='end'):
        self.add_CG_pose(molecule='protein', state=state)

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
        type_filter = AtomFilter(sim.trajectory_order, sim.sequence, sim.particle_ids[molecule], molecule=molecule)
        pose_data = self._create_pose_array(sim, molecule, type_filter.mask)
        if state in ['start', 'end']:
            for run in [sim[self.run.Id]]:
                np_help.copy_fields(pose_data, self.traj_data(run, state, type_filter.mask), ['x','y', 'z'])
                yield pose_data
        elif state == 'full':
            for step_data in self.run.trajectory(molecule=molecule):
                pose_data['x'] = step_data['xyz'][:,0]
                pose_data['y'] = step_data['xyz'][:,1]
                pose_data['z'] = step_data['xyz'][:,2]
                yield pose_data
        else:
            raise AttributeError("state must have value of 'start', 'end' or 'full'.")
            