import numpy as np
import tqdm
import itertools as it
import epitopsy.DXFile as dx
import PolyLibScan.helpers.pymol as pym
import numerics as numerics
import os as _os
import jinja2 as ji

class DensityContainer(object):
    """docstring for DensityContainer"""
    def __init__(self, simulation, monomer_id='all', margin=20.0, resolution=2.0):
        super(DensityContainer, self).__init__()
        self.sim = simulation
        self.monomer_id = monomer_id
        self.margin = margin
        self.resolution = resolution
        self.db_folder = self.sim.db_path.parent.absolute().resolve()
        
        #self.box, self.box_size, self.map = self._create_empty_map()

    @property
    def monomer_id(self):
        return self._monomer_id

    @monomer_id.setter
    def monomer_id(self, value):
        if isinstance(value, int):
            m_id = [value]
        elif isinstance(value, np.ndarray):
            m_id = list(value)
        elif isinstance(value, list):
            m_id = value
        elif value == 'all':
            m_id = list(self.sim._particle_ids['polymer'])
        else:
            raise ValueError('monomer_id must be int or list or "all".')

        if not set(m_id) <= set(self.sim._particle_ids['polymer']):
            raise ValueError("Some Ids are not ids of monomer. Choose from %s" % self.sim._particle_ids['polymer'])
        else:
            self._monomer_id = m_id

    @property
    def margin(self):
        return self._margin

    @margin.setter
    def margin(self, value):
        if not hasattr(self, '_margin'):
            self._margin = value
        else:
            raise UserWarning('This attribute should not be changed after object creation.')
    
    @property
    def resolution(self):
        return self._resolution

    @resolution.setter
    def resolution(self, value):
        if not hasattr(self, '_resolution'):
            self._resolution = value
        else:
            raise UserWarning('This attribute should not be changed after object creation.')


    def _create_empty_map(self):
        '''Set up an empty array for the density and get the size and offset 
        of the protein model.
        '''
        box = self.sim._calc_protein_box(self.margin)
        self.box_size = (box[1] - box[0])
        epi_map = np.zeros(np.ceil(self.box_size/float(self.resolution))+1, np.int32)
        self.box, self.box_size, self.map = box, box_size, epi_map
        
    def _map_path(self):
        '''Create a unique path name for the 3d-density map.
        '''
        prefix = 'map_'
        mono_str = 'ID-%s_' % 'X'.join(map(str, self.monomer_id))
        margin_str = 'M-%.02f_' % (round(self.margin,2))
        resolution_str = 'R-%.02f_' % (round(self.resolution,2))
        suffix = '.dx'
        name = ''.join([prefix, mono_str, margin_str, resolution_str, suffix])
        return self.db_folder.joinpath(name).absolute().resolve()

    def _create_atom_type_filter(self, particle_order, monomer_id=None):
        if not monomer_id:
            monomer_id = self.monomer_id
        mask = np.in1d(particle_order, monomer_id)
        iterator = it.cycle(mask)
        def atom_type_filter(whatever):
            return iterator.next()
        return atom_type_filter

    def _add_run_to_epitopsy_map(self, run_id, monomer_id=None):
        """Add all mapped coordinates to epitopsy grid.
        """
        if not monomer_id:
            monomer_id = self.monomer_id
        particle_order = self.sim._parse.load_traj_type_order()
        type_filter = self._create_atom_type_filter(particle_order)
        offset = self.box[0]
        traj_iterator = self.sim._parse.trajectory(run_id)
        monomer_coords = it.ifilter(type_filter, traj_iterator)
        for coord in it.ifilter(lambda x:numerics.in_box(x, self.box), monomer_coords):
            idx = np.around(((coord - offset) / self.resolution)).astype(np.int)
            self.map[idx[0],idx[1],idx[2]] += 1
        return self.map
    
    def _create_normalized(self, norm_type='max'):
        if norm_type == 'max':
            self.normed_map = self.map / self.map.max().astype(np.float)
        else:
            self.normed_map = self.map

    def create_epitopsy_map(self, norm='max', save=False):
        '''Create a 3-dimensional density map of one monomers. This density map 
        
        input:
        monomer_id: 
        norm: Options: 'max', None [string]
        '''
        self._create_empty_map()
        for run in tqdm.tqdm_notebook(self.sim, leave=False, desc='Adding Job data'):
            self._add_run_to_epitopsy_map(run.Id, self.monomer_id)
        self._create_normalized(norm_type=norm)
        if save:
            self.save()
        return self.box, self.normed_map

    def save(self,):
        '''Save data to the dx format via the DXBox object of epitopsy.
        '''
        dx_obj = dx.DXBox(box=self.normed_map, meshsize=self.box_size/self.normed_map.shape, offset=self.box[0])
        path = self._map_path().as_posix()
        dx_obj.write(path)
        return path

class PymolVisualisation(object):
    atom_names = ['C', 'N', 'O', 'S', 'H']

    def __init__(self, job, protein_path=None):
        self.sim = job
        self._module_folder = _os.path.dirname(__file__)
        self.db_folder = self.sim.db_path.parent.absolute().resolve()
        if protein_path:
            self.protein_path = protein_path
        else:
            potential_pdb_location = list(self.db_folder.parent.parent.joinpath('static').glob('*.pdb'))
            if len(potential_pdb_location) == 1:
                self.protein_path = potential_pdb_location[0].absolute().resolve()
            else:
                self.protein_path = None
        self.pymol_handle = pym.connect(ip='localhost')
    
    def _poly_end_poses(self):
        '''Create generator of data for polymer-pdb via the pymol_polymer.tpl
        Each yield (run) the coordinates of the polymer are updated, while the constant
        information is left unchanged
        '''
        elements = {id_: name for id_, name in zip(self.sim._particle_ids['polymer'], self.atom_names)}
        mask = np.in1d(self.sim._parse.load_traj_type_order(), self.sim._particle_ids['polymer'])
        pose_list_type = [('id', '>i2'), ('ele', '|S1'), ('res_name', '|S3'), ('x', np.float), ('y', np.float), ('z', np.float)]
        pose_data = np.zeros(mask.sum(), dtype=pose_list_type)
        pose_data['id'] = np.arange(1, mask.sum()+1)
        pose_data['ele'] = map(lambda x:elements[x], self.sim[0].coordinates()['end'][mask]['atom_type'])
        pose_data['res_name'] = self.sim.sequence['monomer']
        for run in self.sim:
            coords = run.coordinates()['end'][mask]
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

    def add_dx(self, monomer_id, margin=20.0, resolution=1.5):
        density = DensityContainer(self.sim, monomer_id, margin=margin, resolution=resolution)
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

    def _create_polymer_pdb(self):
        gen = self._poly_end_poses()
        output_path = self.db_folder.joinpath('polymer_.pdb').as_posix()
        with open('%s/pymol_polymer.tpl' % self._module_folder) as f:
            template = ji.Template(f.read())
        with open(output_path, 'w') as f2:
            f2.write(template.render(pymol=self))
        return output_path

    def add_polymers(self):
        polymer_pdb_path = self._create_polymer_pdb()
        self.pymol_handle.load(polymer_pdb_path)
        self.pymol_handle.show_as('sphere', 'polymer_')
        self.pymol_handle.do('cmd.set("sphere_scale", 2.0, "all")')
