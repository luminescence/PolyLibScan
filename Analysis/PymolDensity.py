import density_container as dc
import itertools as it


class PymolDensity(object):
    """docstring for PymolDensity"""
    surface_colors = ['blue', 'brightorange', 'cyan', 'forest', 'orange', 
                      'lightmagenta', 'deepteal', 'oxygen', 'purple', 'yellow']
    def __init__(self, pymol_vis):
        super(PymolDensity, self).__init__()
        self.pymol = pymol_vis
        self.pymol_handle = pymol_vis.pymol_handle
        self.map_data = set([])

    def isosurface(self, lvl=0.5, dx_info=None, color=None):
        '''Show isosurfaces of the created densities.

        if dx_info is not given, then all available 
        densities are visualized in pymol.
        if there dx_info is given or there is only one, 
        the color argument is used, otherwise the colors 
        are taken from the surface_color list.
        '''
        if dx_info:
            if not color: color = self.surface_colors[0]
            self.change_isosurface(dx_info, level=lvl, color=color)
        elif len(self.map_data) == 1:
            if not color: color = self.surface_colors[0]
            first_dx = next(iter(self.map_data))
            self.change_isosurface(first_dx, level=lvl, color=color)
        else:
            for color, density_map in it.izip(self.surface_colors, self.map_data):
                self.change_isosurface(density_map, level=lvl, color=color)

    def change_isosurface(self, density_map, level=0.5, color='blue'):
        self.pymol_handle.do('cmd.isosurface("%s", "%s", level=%f)' % (
            density_map.name, 
            density_map.pymol_map_name, 
            level))    
        self.pymol_handle.do('color %s, %s' % (color, density_map.name))

class Project(PymolDensity):
    """docstring for PymolVisPolyType"""
    def __init__(self, pymol_vis):
        super(Project, self).__init__(pymol_vis)


class Type(PymolDensity):
    """docstring for Type"""
    def __init__(self, pymol_vis):
        super(Type, self).__init__(pymol_vis)
        self.sims = pymol_vis.sims
        
    def add_dx(self, monomer_id='all', filter_specification='type', margin=20.0, resolution=1.5, norm='max'):
        '''Add cumulative densities of all simulations of this type.

        input:
            margin:      margin around the protein in Angstrom [Float]
            resolution:  bin-size in angstrom [Float]
            norm:        options: ['max', 'probability']
        '''
        density = dc.DensityContainer(self.sims, monomer_id, filter_specification=filter_specification, margin=margin,
            resolution=resolution, norm_type=norm)
        dx_path = density._map_path(self.pymol.type_folder)
        if not dx_path.exists():
            density.create_epitopsy_map(norm=norm)
            density.save(self.pymol.type_folder)
        self.pymol_handle.load(dx_path.as_posix())
        density_obj = dc.DensityMap(self.pymol.poly_type.name, len(self.pymol.poly_type.sims), 
                                 monomer_id, margin, resolution, dx_path, norm)
        self.map_data.add(density_obj)
        self.poly_type.project.pymol.map_data.add(density_obj)
        return density_obj


class Job(PymolDensity):
    """docstring for Job"""
    def __init__(self, pymol_vis):
        super(Job, self).__init__(pymol_vis)
        self.sim = pymol_vis.sim

    def add_dx(self, monomer_id='all', filter_specification='type', margin=20.0, resolution=1.5, norm='max'):
        density = dc.DensityContainer(self.sim, monomer_id, filter_specification=filter_specification, margin=margin,
                                      resolution=resolution, norm_type=norm)
        dx_path = density._map_path(self.pymol.db_folder)
        if not dx_path.exists():
            density.create_epitopsy_map(norm=norm)
            density.save(self.pymol.db_folder)
        self.pymol_handle.load(dx_path.as_posix())
        density_obj = dc.DensityMap(self.pymol.sim.poly_type.name, self.pymol.sim.Id, 
                                 monomer_id, margin, resolution, dx_path, norm)
        self.map_data.add(density_obj)
        return density_obj

