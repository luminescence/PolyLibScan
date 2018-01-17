import numpy as np
import tqdm
import itertools as it
import epitopsy.DXFile as dx
import numerics as numerics
import PolyLibScan
from .PymolPose import generate_mask

class DensityContainer(object):
    """docstring for DensityContainer"""
    def __init__(self, simulations, monomer_id='all', margin=20.0, resolution=3.0, 
                 norm_type='max'):
        super(DensityContainer, self).__init__()
        self.sims = simulations
        self.poly_type = self.sims[0].meta['polymer_name']
        self.sim_number = self.sims
        self.monomer_id = monomer_id
        self.margin = margin
        self.resolution = resolution
        self.norm_type = norm_type
        
        self._create_empty_map()

    @property
    def sim_number(self):
        return self._sim_number

    @sim_number.setter
    def sim_number(self, list_):
        if len(list_) == 1:
            self._sim_number = list_[0].Id
        elif len(list_) > 1:
            self._sim_number = -1
        else:
            raise AttributeError('Empty list.')

    @property
    def sims(self):
        return self._sims

    @sims.setter
    def sims(self, value):
        if not isinstance(value, list):
            self._sims = [value]
        else:
            sorted_sims = sorted(value, key=lambda x:x.meta['id'])
            self._sims = sorted_sims

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
            m_id = list(self.sims[0].particle_ids['polymer'])
        else:
            raise ValueError('monomer_id must be int or list or "all".')

        if not set(m_id) <= set(self.sims[0].particle_ids['polymer']):
            raise ValueError("Some Ids are not ids of monomer. Choose from %s" % self.sims[0].particle_ids['polymer'])
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

    def __eq__(self, other, return_differences=False):
        difference_list = []
        if self.margin != other.margin:
            difference_list.append("Margins in maps differ: %.1f and %.1f" % (self.margin, other.margin))
        if self.resolution != other.resolution:
            difference_list.append("Resolution in maps differ: %.1f and %.1f" % (self.resolution, other.resolution))
        if self.monomer_id != other.monomer_id:
            difference_list.append("Monomers in maps differ: %s and %s" % (self.monomer_id, other.monomer_id))
        if (self.box!=other.box).any():
            difference_list.append("Boxes in maps differ: %s and %s" % (self.box, other.box))

        if return_differences:
            return difference_list
        if len(difference_list) == 0:
            return True
        else:
            return False

    def _check_for_difference(self, other):
        differences = self.__eq__(other, return_differences=True)
        if len(differences) == 0:
        	return 'No differences found.'
        else:
            return '\n'.join(differences)

    def __ne__(self, other):
        bool_val = (self == other)
        return not bool_val

    # def __add__(self, other):
    #     if self != other:
    #         raise ValueError(self._check_for_difference(other))
    #     new_map = DensityContainer(self.sims+other.sims, self.monomer_id, self.margin, self.resolution)
    #     return new_map

    # def __iadd__(self, other):
    #     if self != other:
    #         raise ValueError(self._check_for_difference(other))
    #     self.map += other.map
    #     return self

    def _create_empty_map(self):
        '''Set up an empty array for the density and get the size and offset 
        of the protein model.
        '''
        boxes = np.zeros([2,3,len(self.sims)])
        for i, sim in enumerate(self.sims):
            boxes[:,:,i] = sim._calc_protein_box(margin=self.margin)
        # create largest box that surrounds all boxes
        box = np.array([boxes[0,:,:].min(axis=1),boxes[1,:,:].max(axis=1)])
        box_size = (box[1] - box[0])
        epi_map = np.zeros(np.ceil(box_size/float(self.resolution))+1, np.int32)
        self.box, self.box_size, self.map = box, box_size, epi_map

    def _map_path(self, root):
        '''Create a unique path name for the 3d-density map.
        '''
        prefix = 'map'
        p_type = 'PT-%s' % self.poly_type
        s_no = 'sims-%d' % self.sim_number
        mono_str = 'ID-%s' % 'X'.join(map(str, self.monomer_id))
        margin_str = 'M-%.02f' % (round(self.margin,2))
        resolution_str = 'R-%.02f' % (round(self.resolution,2))
        suffix = '.dx'
        name = '_'.join([prefix, p_type, s_no, mono_str, margin_str, 
                         resolution_str, self.norm_type, suffix])
        return root.joinpath(name).absolute().resolve()

    def _create_atom_type_filter(self, particle_order, sim, monomer_id=None):
        if not monomer_id:
            monomer_id = self.monomer_id
        mask = generate_mask(particle_order, monomer_id)
        iterator = it.cycle(mask)
        def atom_type_filter(dummy_variable=True):
            # one variable is required but there is nothing to be passed
            return iterator.next()
        return atom_type_filter

    def _add_run_to_epitopsy_map(self, sim, run_id, monomer_id=None):
        """Add all mapped coordinates to epitopsy grid.
        """
        if not monomer_id:
            monomer_id = self.monomer_id
        particle_order = sim.trajectory_order
        type_filter = self._create_atom_type_filter(particle_order, monomer_id=monomer_id)
        offset = self.box[0]
        traj_iterator = sim._parse.trajectory_load(run_id)
        monomer_coords = it.ifilter(type_filter, traj_iterator)
        for coord in it.ifilter(lambda x:numerics.in_box(x, self.box), monomer_coords):
            idx = np.around(((coord - offset) / self.resolution)).astype(np.int)
            self.map[idx[0],idx[1],idx[2]] += 1
        return self.map
    
    def normalized_map(self):
    	'''public instance of map normalisation'''
        self._create_normalized(norm_type)
        return self.normed_map

    def _create_normalized(self, norm_type='max'):
        if self.norm_type == 'max':
            self.normed_map = self.map / self.map.max().astype(np.float)
        elif self.norm_type == 'probability':
            self.normed_map = self.map / self.map.sum().astype(np.float)
        else:
            self.normed_map = self.map

    def create_epitopsy_map(self, norm='max'):
        '''Create a 3-dimensional density map of one monomers. This density map 
        
        input:
        monomer_id: 
        norm: Options: 'max', None [string]
        '''
        self._create_empty_map()
        number_of_runs = len(self.sims) * len(self.sims[0])
        for sim in tqdm.tqdm_notebook(self.sims, leave=False, desc='Sim:'):
            for run in tqdm.tqdm_notebook(sim, leave=False, desc='Adding Run'):
                self._add_run_to_epitopsy_map(sim, run.Id, self.monomer_id)
        self._create_normalized(norm_type=norm)
        return self.box, self.normed_map

    def save(self, root_folder=None):
        '''Save data to the dx format via the DXBox object of epitopsy.
        '''
        if not hasattr(self, 'normed_map'):
            raise AttributeError("Create a map with create_epitopsy_map()")
        dx_obj = dx.DXBox(box=self.normed_map, meshsize=self.box_size/self.normed_map.shape, offset=self.box[0])
        path = self._map_path(root_folder).as_posix()
        dx_obj.write(path)
        return path

class DensityMap(object):

    def __init__(self, poly_type, sim_no, monomer_id, margin, resolution, path, norm):
        self.poly_type = poly_type
        self.sim_no = sim_no
        self.monomer_id = monomer_id
        self.margin = margin
        self.resolution = resolution
        self.norm = norm
        self.path = path

    @property
    def pymol_map_name(self):
        return self.path.name.rstrip('.dx')

    @property
    def name(self):
        return "%s_density" % ('X'.join(map(str, [
            self.poly_type,
            self.sim_no,
            self.monomer_id,
            self.margin,
            self.resolution,
            self.norm])))

    def to_tuple(self):
        return (self.poly_type, self.sim_no, self.monomer_id, self.margin, self.resolution, self.norm)

    def __eq__(self, other):
        if self.to_tuple() == other.to_tuple():
            return True
        else:
            return False

    def __hash__(self):
        return self.to_tuple().__hash__()

    def __repr__(self):
        return str(self.to_tuple())
