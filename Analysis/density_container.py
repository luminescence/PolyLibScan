import numpy as np
import tqdm
import itertools as it
import epitopsy.DXFile as dx
import numerics as numerics


class DensityContainer(object):
    """docstring for DensityContainer"""
    def __init__(self, simulation, monomer_id='all', margin=20.0, resolution=2.0):
        super(DensityContainer, self).__init__()
        self.sim = simulation
        self.monomer_id = monomer_id
        self.margin = margin
        self.resolution = resolution
        self.db_folder = self.sim.db_path.parent.absolute().resolve()
        
        self._create_empty_map()

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

    def __eq__(self, other):
        if self.margin != other.margin: 
            return False
            raise AttributeError("Margins in maps differ: %f and %f" % (self.margin, other.margin))
        if self.resolution != other.resolution:
            return False
            raise AttributeError("Resolution in maps differ: %f and %f" % (self.resolution, other.resolution))
        if self.monomer_id != other.monomer_id:
            return False
            raise AttributeError("Monomers in maps differ: %s and %s" % (self.monomer_id, other.monomer_id))
        if (self.box!=other.box).any():
            return False
            raise AttributeError("Boxes in maps differ: %s and %s" % (self.box, other.box))
        return True

    def __ne__(self, other):
        bool_val = (self == other)
        return not bool_val

    def __add__(self, other):
        if self != other:
            raise ValueError(self._check_for_difference(other))
        new_map = DensityContainer(self.sim, self.monomer_id, self.margin, self.resolution)
        new_map.map = self.map + other.map
        new_map.box = self.box
        return new_map

    def __iadd__(self, other):
        if self != other:
            raise ValueError(self._check_for_difference(other))
        self.map += other.map
        return self

    def _check_for_difference(self, other):
        if self.margin != other.margin: 
            return "Margins in maps differ: %.1f and %.1f" % (self.margin, other.margin)
        if self.resolution != other.resolution:
            return "Resolution in maps differ: %.1f and %.1f" % (self.resolution, other.resolution)
        if self.monomer_id != other.monomer_id:
            return "Monomers in maps differ: %s and %s" % (self.monomer_id, other.monomer_id)
        if (self.box!=other.box).any():
            return "Boxes in maps differ: %s and %s" % (self.box, other.box)
        else:
        	return 'No differences found.'

    def _create_empty_map(self):
        '''Set up an empty array for the density and get the size and offset 
        of the protein model.
        '''
        box = self.sim._calc_protein_box(self.margin)
        box_size = (box[1] - box[0])
        epi_map = np.zeros(np.ceil(box_size/float(self.resolution))+1, np.int32)
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
    
    def normalized_map(self, norm_type='max'):
    	'''public instance of map normalisation'''
        self._create_normalized(norm_type)
        return self.normed_map

    def _create_normalized(self, norm_type='max'):
        if norm_type == 'max':
            self.normed_map = self.map / self.map.max().astype(np.float)
        elif nor_type == 'probability':
            self.normed_map = self.map / self.map.sum().astype(np.float)
        else:
            self.normed_map = self.map

    def create_epitopsy_map(self, norm='max', save=False):
        '''Create a 3-dimensional density map of one monomers. This density map 
        
        input:
        monomer_id: 
        norm: Options: 'max', None [string]
        '''
        self._create_empty_map()
        for run in tqdm.tqdm_notebook(self.sim, leave=False, desc='Adding Run'):
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