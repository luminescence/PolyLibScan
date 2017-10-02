import yaml
from numpy import zeros, array
# Own modules
from lmp_types import *
import lmp_force_field as force_field
import container
import PolyLibScan.helpers.idGenerator as idGenerator


class Environment(object):
    '''Acts as a Container for then entire 
    simulation. It stores all particles, interactions
    and relations between all particles.
    '''

    def __init__(self, config_file):
        self.new_id = idGenerator.IdGen()

        config = self._get_config(config_file)
        # Pair Type Container
        self.pair_type = container.Container('pair')
        self.load_pair_types(config['Pairs'])
        
        # Force Field Container
        self.ff = self.create_force_fields(self.pair_type)
        
        self.atom_type = container.Container('atom')
        for field in self.ff.values():
            self.atom_type.register(field)
        self.load_atom_types(config['Atoms'])
        
        self.bond_type = container.Container('bond')
        self.load_bond_types(config['Bonds'])

        self.angle_type = container.Container('angle')
        self.load_angle_types(config['Angles'])

        self.dihedral_type = container.Container('dihedral')
        self.load_dihedral_types(config['Dihedrals'])        
        # Read in Types and Globals

        self.globals = dict()
        self.load_globals(config['globals'])
        
        self.box = zeros([3,2])

        self.molecules = container.Container('molecules')

    def _get_config(self, config_path):
        with open(config_path) as f:
            config = yaml.load(f)
        return config

    def create_force_fields(self, pair_types):
        '''only force fields that are single_typed need to be addressed,
        in the lammps data file.
        '''
        ff = {}
        for name, field_type in pair_types.items():
            if field_type.single_type_parametrised:
                ff[name] = force_field.Interaction(self, field_type)
        return ff

    def load_pair_types(self, config):
        for name, pair_parameters in config.items():
            self.pair_type[name] = PairType.from_dict(pair_parameters, name=name)

    def load_atom_types(self, config):
        for name, particle in config.items():
            for subname, parameter in particle.items():
                combinedname = name + '_' + subname
                self.atom_type.define_type(combinedname, AtomType(combinedname, parameter))

    def load_bond_types(self, config):
        for name, values in config.items():
            self.bond_type.define_type(name, BondType.from_dict(values, name=name))
    def load_angle_types(self, config):
        for name, values in config.items():
            self.angle_type.define_type(name, AngleType.from_dict(values, name=name))

    def load_dihedral_types(self, config):
        for name, values in config.items():
            self.dihedral_type.define_type(name, DihedralType.from_dict(values, name=name))

    def load_globals(self, parameters):
        valid_keys = ['angle_style', 'atom_style', 'bond_style', 
                      'box_margin', 'affinity_file', 'pair_style']
        # Check if all required parameters are present
        for key in valid_keys:
            if key not in parameters.keys():
                print 'Warning: %s parameter is not listed in parameter file.' % key
        
        for key, value in parameters.items():
            self.globals[key] = value
        

    def add_molecule(self, molecule):
        self.molecules[molecule.name] = molecule
        self.molecules.register(self)

    def calc_box(self, add_margin=False):
        margin = 0
        if add_margin:
            margin = float(self.globals['box_margin'])
        if len(self.molecules) > 0:
            boxes = array([mol.box for mol in self.molecules.values()])
        else:
            boxes = array([zeros([3,2])])
        for i in range(3):
            self.box[i,0] = boxes[:,i,0].min()-margin
            self.box[i,1] = boxes[:,i,1].max()+margin

    def update(self, value):
        self.calc_box()
