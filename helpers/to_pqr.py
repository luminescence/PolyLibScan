"""This module creates a pqr file from a single simulation run object 
of an analysis project.

"""
import os as _os
import numpy as np
import jinja2 as ji
import PolyLibScan.helpers.numpy_helpers as np_help

def _get_pdb_template(name):
    with open('%s/%s' % (_os.path.dirname(__file__),name)) as f:
        template = ji.Template(f.read())
    return template

class to_pqr(object):

    atom_names = 100*['C', 'N', 'O', 'S', 'H']
    
    pqr_template = _get_pdb_template('pqr.tpl')

    """docstring for to_pqr"""
    def __init__(self, run):
        super(to_pqr, self).__init__()
        self.run = run
        self.job = self.run.job

    def save_pqr(self, file_name, molecule='polymer', state='end'):
        data = self.pqr(molecule=molecule, state=state)
        output_path = self._save_data(file_name, data)
        return output_path

    def _save_data(self, file_name, pdb_data):
        with open(file_name, 'w') as f:
            f.write(self.pqr_template.render(pymol=pdb_data))
        return file_name

    def _create_pqr_array(self, molecule, mask):
        '''Create an array to hold the information of the 
        run conformation and also store the static portion 
        in it. This is basically everything except the coordinates.
        '''
        pose_list_type = [('id', '>i2'), ('ele', '|S1'), ('res_name', '|S3'), 
                          ('chain', '|S1'), ('res_id', '>i8'), ('iCode', '|S1'),
                          ('x', np.float), ('y', np.float), ('z', np.float),
                          ('charge', np.float), ('radius', np.float)]
        elements = dict(zip(self.job.particle_ids[molecule], self.atom_names))
        pose_data = np.zeros(mask.sum(), dtype=pose_list_type)
        pose_data['id'] = np.arange(1, mask.sum()+1)
        pose_data['ele'] = map(lambda x:elements[x], self.run.coordinates()['end'][mask]['atom_type'])
        mask = np.in1d(self.job.particle_list['p_id'], self.run.job.particle_ids[molecule])
        res_name = self.job.particle_list['name'][mask].split('_').replace('-', '')
        pose_data['res_name'] = res_name
        pose_data['chain'] = self.job.particle_list['chain'][mask]
        pose_data['charge'] = self.job.particle_list['charge'][mask]
        pose_data['res_id'] = self.job.particle_list['res_id'][mask]
        pose_data['iCode'] = self.job.particle_list['iCode'][mask]
        pose_data['radius'] = np.array([self.job.project.parameters['Atoms'][atom_name][atom_bead]['radius'] 
                                            for atom_name,atom_bead in self.job.particle_list['name'][mask].split('_')])
        return pose_data

    def pqr(self, molecule='polymer', state='end'):
        '''Create generator of data for polymer-pdb via the pymol_polymer.tpl
        Each yield (run) the coordinates of the polymer are updated, while the constant
        information is left unchanged.
        '''
        mask = np.in1d(self.job.trajectory_order, self.job.particle_ids[molecule])
        pose_data = self._create_pqr_array(molecule, mask)
        if state in ['start', 'end']:
            np_help.copy_fields(pose_data, self.run.coordinates()[state][mask], ['x','y', 'z'])
            return pose_data
        else:
            raise AttributeError("state must have value of 'start', 'end'.")
