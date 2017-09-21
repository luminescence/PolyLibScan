import numpy as np
import os as _os
import jinja2 as ji
# own module
import PolyLibScan.helpers.pymol as pym
import PolyLibScan.helpers.numpy as np_help
import density_container as dc
import pathlib2 as pl
import itertools as it


class PymolVisualisation(object):
    '''Parent class of pymol Visualisation
    '''
    atom_names = ['C', 'N', 'O', 'S', 'H']
    surface_colors = ['blue', 'brightorange', 'cyan', 'forest', 'orange', 
                      'lightmagenta', 'deepteal', 'oxygen', 'purple', 'yellow']
    def __init__(self, protein_path=None):
        self.protein_path = protein_path
        self._module_folder = _os.path.dirname(__file__)
        self.pymol_handle = pym.connect(ip='localhost')
        self.map_data = set([])

    def _init_protein_path(self, protein_path, search_path=None):
        '''Set the path of the pdb model to protein_path or
        if protein_path is set to None look for the most likely 
        path relative to the job's database path. 
        If that fails, return None.
        '''
        if protein_path:
            return protein_path
        elif search_path:
            potential_pdb_location = list(search_path.glob('*.pdb'))
            if len(potential_pdb_location) == 1:
                return potential_pdb_location[0].absolute().resolve()
            else:
                raise AttributeError('Could not find pdb file at %s.' % search_path)
        else:
            raise AttributeError('No pdb file or file location given')

    def _poly_poses(self, sim, state='end'):
        '''Create generator of data for polymer-pdb via the pymol_polymer.tpl
        Each yield (run) the coordinates of the polymer are updated, while the constant
        information is left unchanged
        '''
        if state not in ['start', 'end']:
            raise AttributeError("state must have value of 'start' or 'end'.")
        mask = np.in1d(sim._parse.load_traj_type_order(), sim.particle_ids['polymer'])
        pose_data = self._create_poly_pose_array(sim, mask)
        for run in sim:
            np_help.copy_fields(pose_data, self.traj_data(run, state, mask), ['x','y', 'z'])
            yield pose_data

    def _create_poly_pose_array(self, sim, mask):
        '''
        '''
        pose_list_type = [('id', '>i2'), ('ele', '|S1'), ('res_name', '|S3'), 
                          ('x', np.float), ('y', np.float), ('z', np.float)]
        elements = dict(zip(sim.particle_ids['polymer'], self.atom_names))
        pose_data = np.zeros(mask.sum(), dtype=pose_list_type)
        pose_data['id'] = np.arange(1, mask.sum()+1)
        pose_data['ele'] = map(lambda x:elements[x], sim[0].coordinates()['end'][mask]['atom_type'])
        pose_data['res_name'] = sim.sequence['monomer']


    def traj_data(self, run, state, mask):
        return run.coordinates()[state][mask]

    def add_polymers(self, state='end'):
        polymer_pdb_path, name = self._create_polymer_pdb(state=state)
        self.pymol_handle.load(polymer_pdb_path)
        self.pymol_handle.show_as('sphere', name)
        self.pymol_handle.do('cmd.set("sphere_scale", 2.0, "all")')


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

    def isosurface(self, lvl=0.5, dx_info=None, color=None):
        if dx_info:
            if not color: color = self.surface_colors[0]
            self.change_isosurface(dx_info, level=lvl, color=color)
        else:
            for color, density_map in it.izip(self.surface_colors, self.map_data):
                self.change_isosurface(density_map, level=lvl, color=color)

    def change_isosurface(self, density_map, level=0.5, color='blue'):
        self.pymol_handle.do('cmd.isosurface("%s", "%s", level=%f)' % (
            density_map.name, 
            density_map.pymol_map_name, 
            level))    
        self.pymol_handle.do('color %s, %s' % (color, density_map.name))


class PymolVisProject(PymolVisualisation):
    """docstring for PymolVisPolyType"""
    def __init__(self, project, protein_path=None):
        self.project = project
        self._default_pdb_folder = self.project.path.joinpath('static')
        protein_path = self._init_protein_path(
                        protein_path, 
                        search_path=self._default_pdb_folder)
        super(PymolVisProject, self).__init__(protein_path=protein_path)

    
class PymolVisPolyType(PymolVisualisation):
    """docstring for PymolVisPolyType"""
    def __init__(self, poly_type, protein_path=None):
        self.poly_type = poly_type
        self.sims = sorted(poly_type.sims, key=lambda x:x.Id)
        self.sim = self.sims[0]
        first_sim_folder = self.sim.db_path.parent
        self._default_pdb_folder = self.poly_type.project.path.joinpath('static')
        protein_path = self._init_protein_path(
                        protein_path, 
                        search_path=self._default_pdb_folder)
        self.type_folder = self.create_polytype_folder(first_sim_folder)
        super(PymolVisPolyType, self).__init__(protein_path)

    def create_polytype_folder(self, sim_path):
        '''Create the polytype folder in the first sim folder
        (order lexilographic)
        '''
        polytype_folder = sim_path.joinpath(self.poly_type.name)
        polytype_folder.mkdir(exist_ok=True)
        return polytype_folder

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

    def add_dx(self, monomer_id='all', margin=20.0, resolution=1.5, norm='max'):
        density = dc.DensityContainer(self.sims, monomer_id, margin=margin, 
            resolution=resolution, norm_type=norm)
        dx_path = density._map_path(self.type_folder)
        if not dx_path.exists():
            density.create_epitopsy_map(norm=norm)
            density.save(self.type_folder)
        self.pymol_handle.load(dx_path.as_posix())
        density_obj = DensityMap(self.poly_type.name, len(self.poly_type.sims), 
                                 monomer_id, margin, resolution, dx_path, norm)
        self.map_data.add(density_obj)
        self.poly_type.project.pymol.map_data.add(density_obj)
        return density_obj

class PymolVisJob(PymolVisualisation):
    """docstring for PymolVisJob"""
    def __init__(self, job, protein_path=None):
        self.sim = job
        self._default_pdb_folder = self.sim.project.path.joinpath('static')
        protein_path = self._init_protein_path(
                        protein_path, 
                        search_path=self._default_pdb_folder)
        self.db_folder = self.sim.db_path.parent.absolute().resolve()
        super(PymolVisJob, self).__init__(protein_path)
        
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

    def add_dx(self, monomer_id='all', margin=20.0, resolution=1.5, norm='max'):
        density = dc.DensityContainer(self.sim, monomer_id, margin=margin, 
            resolution=resolution, norm_type=norm)
        dx_path = density._map_path(self.db_folder)
        if not dx_path.exists():
            density.create_epitopsy_map(norm=norm)
            density.save(self.db_folder)
        self.pymol_handle.load(dx_path.as_posix())
        density_obj = DensityMap(self.sim.poly_type.name, self.sim.Id, 
                                 monomer_id, margin, resolution, dx_path, norm)
        self.map_data.add(density_obj)
        return density_obj

class PymolVisRun(PymolVisualisation):

    def __init__(self, run, protein_path=None):
        self.run = run
        self.sim = run.job
        self._default_pdb_folder = self.sim.project.path.joinpath('static')
        protein_path = self._init_protein_path(
                        protein_path, 
                        search_path=self._default_pdb_folder)
        self.db_folder = self.sim.db_path.parent.absolute().resolve()
        super(PymolVisRun, self).__init__(protein_path)


    def _create_polymer_pdb(self, state='end'):
        if state in ['start', 'end', 'full']:
            gen = self._poly_poses(self.sim, state=state)
        else:
            raise AttributeError("state supports only 'start' or 'end' as values.")
        file_name = '%s__run_%d-%s_poses.pdb' % (self.sim.meta['id'], self.run.Id, state)
        output_path = self.db_folder.joinpath(file_name).as_posix()
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
        pose_data = self._create_poly_pose_array(sim, mask)
        if state in ['start', 'end']:
            for run in [sim[self.run.Id]]:
                np_help.copy_fields(pose_data, self.traj_data(run, state, mask))
                yield pose_data
        if state == 'full':
            for step_data in self.run.trajectory():
                pose_data['x'] = step_data['xyz'][0]
                pose_data['y'] = step_data['xyz'][1]
                pose_data['z'] = step_data['xyz'][2]
                yield pose_data
        else:
            raise AttributeError("state must have value of 'start' or 'end'.")

    def traj_data(self, run, state, mask):
            return run.coordinates()[state][mask]
        

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
