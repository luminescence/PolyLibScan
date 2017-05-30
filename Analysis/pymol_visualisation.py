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
        self.db_folder = self.sim.db_path.parent.absolute().resolve()
        self.margin = margin
        self.resolution = resolution
        
        self.box, self.box_size, self.map = self._create_empty_map()

    @property
    def monomer_id(self):
        return self._monomer_id

    @monomer_id.setter
    def monomer_id(self, value):
        if isinstance(value, int):
            m_id = [monomer_id]
        elif isinstance(value, np.ndarray):
            m_id = list(value)
        elif isinstance(value, list):
            m_id = value
        elif value == 'all':
            m_id = list(self.sim._particle_ids['polymer'])
        else:
            raise ValueError('monomer_id must be int or list or "all".')
        if len(set(m_id) - set(self.sim._particle_ids['polymer'])) > 0:
            raise ValueError("Some Ids are not ids of monomer. Choose from %s" % self.sim._particle_ids['polymer'])
        else:
            self._monomer_id = m_id

    @property
    def margin(self):
        return self._margin

    @margin.setter
    def margin(self):
        raise UserWarning('This attribute should not be changed after object creation.')
    
    @property
    def resolution(self):
        return self._resolution

    @resolution.setter
    def resolution(self):
        raise UserWarning('This attribute should not be changed after object creation.')


    def _create_empty_map(self):
        '''Set up an empty array for the density and get the size and offset 
        of the protein model.
        '''
        box = self.sim._calc_protein_box(self.margin)
        self.box_size = (box[1] - box[0])
        epi_map = np.zeros(np.ceil(self.box_size/float(self.resolution))+1, np.int32)
        return box, box_size, epi_map
        
    def _map_path(self, monomer_id):
        '''Create a unique path name for the 3d-density map.
        '''
        prefix = 'map_'
        mono_str = 'ID-%s_' % 'X'.join(map(str, monomer_id))
        margin_str = 'M-%.02f_' % (round(self.margin,2))
        resolution_str = 'R-%.02f_' % (round(self.resolution,2))
        suffix = '.dx'
        name = ''.join([prefix, mono_str, margin_str, resolution_str, suffix])
        return self.db_folder.joinpath(name).absolute().resolve()

    def _create_atom_type_filter(self, particle_order, monomer_id):
        mask = np.in1d(particle_order, monomer_id)
        iterator = it.cycle(mask)
        def atom_type_filter(whatever):
            return iterator.next()
        return atom_type_filter

    def _add_run_to_epitopsy_map(self, run_id, monomer_id):
        """Add all mapped coordinates to epitopsy grid.
        """
        particle_order = self.sim._parse.load_traj_type_order()
        type_filter = self._create_atom_type_filter(particle_order, self.monomer_id)
        offset = self.box[0]
        traj_iterator = self.sim._parse.trajectory(run_id)
        monomer_coords = it.ifilter(type_filter, traj_iterator)
        for coord in it.ifilter(lambda x:numerics.in_box(x, self.box), monomer_coords):
            idx = np.around(((coord - offset) / self.resolution)).astype(np.int)
            self.map[idx[0],idx[1],idx[2]] += 1
        return self.map
    
    def create_epitopsy_map(self, monomer_id, margin=20, resolution=1.0, norm=None, save=False):
        '''Create a 3-dimensional density map of one monomers. This density map 
        
        input:
        monomer_id: 
        norm: if true, rescales all values by deviding through the maximum value of the map. [bool]
        '''
        box, epi_box = self._create_empty_epitopsy_map(self.margin, self.resolution)
        for run in tqdm.tqdm_notebook(self.sim, leave=False, desc='Adding Job data'):
            self._add_run_to_epitopsy_map(run.Id, monomer_id)
        if norm:
            normed_epi_box = epi_box / epi_box.max().astype(np.float)
        else:
            normed_epi_box = epi_box
        if save:
            box_size = box[1] - box[0]
            dx_obj = dx.DXBox(box=normed_epi_box, meshsize=box_size/normed_epi_box.shape, offset=box[0])
            path = self._map_path(monomer_id, margin, resolution)
            dx_obj.write(path.as_posix())
        return box, normed_epi_box

    def save(self):
        '''Save data to the dx format via the DXBox object of epitopsy.
        '''
        box_size = (self.box[1] - self.box[0])
        dx_obj = dx.DXBox(box=self.map, meshsize=box_size/self.map.shape, offset=self.box[0])
        dx_obj.write(self._map_path().as_posix())

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
    
    def save_density_maps(self, margin=20, resolution=1.0):
        complete_box = None
        complete_path = self._map_path(self.sim._particle_ids['polymer'], margin, resolution)
        for monomer_id in tqdm.tqdm_notebook(self.sim._particle_ids['polymer'], desc='Mono Maps'):
            # generate standardised file path
            path = self._map_path(monomer_id, margin, resolution)
            if path.exists() and complete_path.exists():
                continue
            box, map_ = self.create_epitopsy_map(monomer_id, margin=margin, resolution=resolution, norm=True, save=False)
            if complete_box == None:
                complete_box = box
                complete_map = np.zeros_like(map_)
            complete_map += map_
        # generating map of all monomers
        if not complete_path.exists():
            complete_box_size = complete_box[1] - complete_box[0]
            dx_obj = dx.DXBox(box=complete_map, meshsize=complete_box_size/complete_map.shape, offset=complete_box[0])
            dx_obj.write(complete_path.as_posix())

    def _poly_end_poses(self):
        elements = {id_: name for id_, name in zip(self.sim._particle_ids['polymer'], self.atom_names)}
        mask = np.in1d(self.sim._parse.load_traj_type_order(), self.sim._particle_ids['polymer'])
        pose_list_type = [('id', '>i2'), ('ele', '|S1'), ('x', np.float), ('y', np.float), ('z', np.float)]
        pose_data = np.zeros(mask.sum(), dtype=pose_list_type)
        pose_data['id'] = np.arange(1, mask.sum()+1)
        pose_data['ele'] = map(lambda x:elements[x], self.sim[0].coordinates()['end'][mask]['atom_type'])
        for run in self.sim:
            coords = run.coordinates()['end'][mask]
            pose_data['x'] = coords['x']
            pose_data['y'] = coords['y']
            pose_data['z'] = coords['z']
            yield pose_data

    def setup(self):
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
        if monomer_id == 'all':
            monomer_id = self.sim._particle_ids['polymer']
        elif not (isinstance(monomer_id, list) or isinstance(monomer_id, np.ndarray)):
            monomer_id = [monomer_id]
        path = self._map_path(monomer_id, margin=margin, resolution=resolution)
        if not path.exists():
            self.create_epitopsy_map(monomer_id, margin=margin, resolution=resolution, norm=True, save=True)
        self.pymol_handle.load(path.as_posix())
        pymol_name = path.name.rstrip('.dx')
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