import os as _os
import pathlib2 as pl
# own module
import PolyLibScan.helpers.pymol as pym
import PymolDensity
import PymolPose

class PymolVisualisation(object):
    '''Parent class of pymol Visualisation
    '''
    def __init__(self, protein_path=None):
        self.protein_path = protein_path
        self.pymol_handle = pym.connect(ip='localhost')

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
                raise AttributeError('Found more that one pdb file at %s.' % search_path)
        else:
            raise AttributeError('No pdb file or file location given')

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
        """Visualize the .pdb file of the protein as spheres. Note that the visualization represents the crystal
        structure, not the conformation of the protein during the MD simulation!"""
        if not level in ['full', 'cg']:
            raise AttributeError("level must be either 'cg' or 'full'")
        protein_name = self.protein_path.name[:-4]
        if level == 'cg':
            self.pymol_handle.select('backbone', protein_name + ' and name ca')
            self.pymol_handle.hide('(%s)' % protein_name)
            self.pymol_handle.show_as('sphere', 'backbone')
        if level == 'full':
            self.pymol_handle.show_as('cartoon', protein_name)
            self.pymol_handle.show('sticks', 'active_site')


class PymolVisProject(PymolVisualisation):
    """docstring for PymolVisProject"""
    def __init__(self, project, protein_path=None):
        self.project = project
        job_example = self.project.jobs[0]
        self.sim = job_example
        self._default_pdb_folder = self.project.path.joinpath('static')
        if not protein_path:
            pdb_name = job_example._parse.misc['protein']
            protein_path = self._default_pdb_folder.joinpath('%s.pdb' % pdb_name)
        protein_path = self._init_protein_path(
                        protein_path, 
                        search_path=self._default_pdb_folder)
        super(PymolVisProject, self).__init__(protein_path=protein_path)
        self.density = PymolDensity.PymolDensity(self)

    
class PymolVisPolyType(PymolVisualisation):
    """docstring for PymolVisPolyType"""
    def __init__(self, poly_type, protein_path=None):
        self.poly_type = poly_type
        self.sims = sorted(poly_type.sims, key=lambda x:x.Id)
        self.sim = self.sims[0]
        first_sim_folder = self.sim.db_path.parent
        self._default_pdb_folder = self.poly_type.project.path.joinpath('static')
        if not protein_path:
            pdb_name = self.sim._parse.misc['protein']
            protein_path = self._default_pdb_folder.joinpath('%s.pdb' % pdb_name)
        protein_path = self._init_protein_path(
                        protein_path, 
                        search_path=self._default_pdb_folder)
        self.type_folder = self.create_polytype_folder(first_sim_folder)
        super(PymolVisPolyType, self).__init__(protein_path)
        self.density = PymolDensity.Type(self)
        self.pose = PymolPose.Type(self)

    def create_polytype_folder(self, sim_path):
        '''Create the polytype folder in the first sim folder
        (order lexilographic)
        '''
        polytype_folder = sim_path.joinpath(self.poly_type.name)
        polytype_folder.mkdir(exist_ok=True)
        return polytype_folder


class PymolVisJob(PymolVisualisation):
    """docstring for PymolVisJob"""
    def __init__(self, job, protein_path=None):
        self.sim = job
        self._default_pdb_folder = self.sim.project.path.joinpath('static')
        if not protein_path:
            pdb_name = self.sim._parse.misc['protein']
            protein_path = self._default_pdb_folder.joinpath('%s.pdb' % pdb_name)
        protein_path = self._init_protein_path(
                        protein_path, 
                        search_path=self._default_pdb_folder)
        self.db_folder = self.sim.db_path.parent.absolute().resolve()
        super(PymolVisJob, self).__init__(protein_path)
        self.density = PymolDensity.Job(self)
        self.pose = PymolPose.Job(self)


class PymolVisRun(PymolVisualisation):

    def __init__(self, run, protein_path=None):
        self.run = run
        self.sim = run.job
        self._default_pdb_folder = self.sim.project.path.joinpath('static')
        if not protein_path:
            pdb_name = self.sim._parse.misc['protein']
            protein_path = self._default_pdb_folder.joinpath('%s.pdb' % pdb_name)
        protein_path = self._init_protein_path(
                        protein_path, 
                        search_path=self._default_pdb_folder)
        self.db_folder = self.sim.db_path.parent.absolute().resolve()
        super(PymolVisRun, self).__init__(protein_path)
        self.pose = PymolPose.Run(self)
