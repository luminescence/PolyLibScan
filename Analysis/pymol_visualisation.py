import numpy as np
import os as _os
import jinja2 as ji
# own module
import PolyLibScan.helpers.pymol as pym
import density_container as dc

class PymolVisualisation(object):
    atom_names = ['C', 'N', 'O', 'S', 'H']

    def __init__(self, job, protein_path=None):
        self.sim = job
        self._module_folder = _os.path.dirname(__file__)
        self.db_folder = self.sim.db_path.parent.absolute().resolve()
        self.protein_path = self._init_protein_path(protein_path)
        self.pymol_handle = pym.connect(ip='localhost')
    
    def _init_protein_path(self, protein_path):
        if protein_path:
            return protein_path
        else:
            potential_pdb_location = list(self.db_folder.parent.parent.joinpath('static').glob('*.pdb'))
            if len(potential_pdb_location) == 1:
                return potential_pdb_location[0].absolute().resolve()
            else:
                return None

    def _poly_poses(self, state='end'):
        '''Create generator of data for polymer-pdb via the pymol_polymer.tpl
        Each yield (run) the coordinates of the polymer are updated, while the constant
        information is left unchanged
        '''
        if state not in ['start', 'end']:
            raise AttributeError("state must have value of 'start' or 'end'.")
        elements = {id_: name for id_, name in zip(self.sim._particle_ids['polymer'], self.atom_names)}
        mask = np.in1d(self.sim._parse.load_traj_type_order(), self.sim._particle_ids['polymer'])
        pose_list_type = [('id', '>i2'), ('ele', '|S1'), ('res_name', '|S3'), ('x', np.float), ('y', np.float), ('z', np.float)]
        pose_data = np.zeros(mask.sum(), dtype=pose_list_type)
        pose_data['id'] = np.arange(1, mask.sum()+1)
        pose_data['ele'] = map(lambda x:elements[x], self.sim[0].coordinates()[state][mask]['atom_type'])
        pose_data['res_name'] = self.sim.sequence['monomer']
        for run in self.sim:
            coords = run.coordinates()[state][mask]
            pose_data['x'] = coords['x']
            pose_data['y'] = coords['y']
            pose_data['z'] = coords['z']
            yield pose_data

    def setup(self):
        '''Setup a new pymol scene. First, delete any previous setup.
        Load the protein in pymol in cartoon style,
        and highlight the active site.
        '''
        self.pymol_handle.delete('all')
        self.pymol_handle.load(self.protein_path.as_posix())
        self.pymol_handle.show_as('cartoon')
        self.pymol_handle.color('green')
        selection_string = ['(chain %s and resi %d)' % (chain, id_) 
                            for chain,id_ in self.sim.active_site[['chain', 'pdb_id']]]
        self.pymol_handle.select('active_site', ' or '.join(selection_string))
        self.pymol_handle.show('sticks', 'active_site')
        self.pymol_handle.color('red', 'active_site')

    def set_protein_detail(self, level='cg'):
        if not level in ['full', 'cg']:
            raise AttributeError("level must be either 'cg' or 'full'")
        protein_name = _os.path.base_name(self.protein)[:-4]
        if level == 'cg':
            self.pymol_handle.select('backbone', protein_name + ' and name ca')
            self.pymol_handle.hide('(%s)' % protein_name)
            self.pymol_handle.show_as('sphere', 'backbone')
        if level == 'full':
            self.pymol_handle.show_as('cartoon', protein_name)
            self.pymol_handle.show('sticks', 'active_site')

    def add_dx(self, monomer_id, margin=20.0, resolution=1.5):
        density = dc.DensityContainer(self.sim, monomer_id, margin=margin, resolution=resolution)
        dx_path = density._map_path()
        if not dx_path.exists():
            density.create_epitopsy_map(norm='max', save=True)
        self.pymol_handle.load(dx_path.as_posix())
        pymol_name = dx_path.name.rstrip('.dx')
        return monomer_id, pymol_name

    def show_isosurface(self, dx_info, level=0.5):
        monomer_id, pymol_name = dx_info
        surface_name = "%s_%s_density" % (self.sim.meta['polymer_name'], 
                                          'X'.join(map(str, monomer_id)))
        self.pymol_handle.do('cmd.isosurface("%s", "%s", level=%f)' % (surface_name, pymol_name, 
                                                                               level))
        self.pymol_handle.do('color blue, %s' % surface_name)

    def _create_polymer_pdb(self, state='end'):
        if state == 'start':
            gen = self._poly_poses(state=state)
            file_name = 'start_poses.pdb'
        elif state == 'end':
            gen = self._poly_poses(state=state)
            file_name = 'end_poses.pdb'
        else:
            raise AttributeError("state supports only 'start' or 'end' as values.")
        output_path = self.db_folder.joinpath(file_name).as_posix()
        pymol_name = file_name[:-4]
        with open('%s/pymol_polymer.tpl' % self._module_folder) as f:
            template = ji.Template(f.read())
        with open(output_path, 'w') as f2:
            f2.write(template.render(pymol=gen))
        return output_path, pymol_name

    def add_polymers(self, state='end'):
        polymer_pdb_path, name = self._create_polymer_pdb(state=state)
        self.pymol_handle.load(polymer_pdb_path)
        self.pymol_handle.show_as('sphere', name)
        self.pymol_handle.do('cmd.set("sphere_scale", 2.0, "all")')
