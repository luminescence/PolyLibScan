from lmp_types import *
import lmp_force_field as force_field
import yaml
from numpy import zeros, array

class Container(object):
    """Stores types very similar to a 
    dictionary with the added features of 
    automatic id generation for every item and
    hooks for updating dependent objects."""

    def __init__(self, type_name):
        super(Container, self).__init__()
        self.name = type_name
        self.new_type_id = self._idGenerator()
        self._types = {}
        self._defined_types = {}
        self._registered = []
        
    def __repr__(self):
        out = ['Type Category: %s' % self.name]
        for key, type_ in self._types.iteritems():
            out.append('%s: %s' % (key, repr(type_)))
        return '\n'.join(out)

    def __getitem__(self, key):
        try:
            return self._types[key]
        except KeyError:
            if key in self._defined_types.keys():
                # adding defined type to available types using 
                # __setitem__ method.
                self[key] = self._defined_types[key]
                return self._types[key]
            else:
                raise KeyError('Type %s not defined. Add type to config file.' % key)


    def __contains__(self, key):
        if key in self._types.keys() or key in self._defined_types.keys():
            return True
        else:
            return False

    def __len__(self):
        return len(self._types)

    def __iter__(self):
        return iter(self._types)

    def __delitem__(self, key):
        del self._types[key]

    def register(self, obj):
        self._registered.append(obj)

    def _update_hook(self, value):
        '''update registered objects via their update function
        when new ids are added to the atom type list.
        '''
        for obj in self._registered:
            obj.update(value)

    def __setitem__(self, key, value):
        type_id = self.new_type_id.next()
        if key in self._types.keys():
            unique_key = key + str(type_id)
        else:
            unique_key = key
        self._types[unique_key] = value
        self._types[unique_key].Id = type_id
        # update hooked object
        self._update_hook(value)

    def define_type(self, key, value):
        self._defined_types[key] = value

    def keys(self):
        return self._types.keys()

    def values(self):
        return self._types.values()

    def items(self):
        return self._types.items()

    def iteritems(self):
        return key,value in self._types.iteritems()

    def _idGenerator(self):
        '''generate ids of an infinite loop
        starting with 1.
        '''
        counter = 1
        while True:
            yield counter
            counter += 1

class Environment(object):
    '''Acts as a Container for then entire 
    simulation. It stores all particles, interactions
    and relations between all particles.
    '''

    def __init__(self, config_file):
        self.new_id = IdGen()

        config = self._get_config(config_file)
        # Pair Type Container
        self.pair_type = Container('pair')
        self.load_pair_types(config['Pairs'])
        
        # Force Field Container
        self.ff = self.create_force_fields(self.pair_type)
        
        self.atom_type = Container('atom')
        for field in self.ff.values():
            self.atom_type.register(field)
        self.load_atom_types(config['Atoms'])
        
        self.bond_type = Container('bond')
        self.load_bond_types(config['Bonds'])

        self.angle_type = Container('angle')
        self.load_angle_types(config['Angles'])

        self.dihedral_type = Container('dihedral')
        self.load_dihedral_types(config['Dihedrals'])        
        # Read in Types and Globals

        self.globals = dict()
        self.load_globals(config['globals'])
        
        self.box = zeros([3,2])

        self.molecules = Container('molecules')

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
        for name, values in config.items():
            self.atom_type.define_type(name, AtomType(name, values))

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

class IdGen(object):
    '''generates ids for every key you throw at it.
    If a key is not registered, it will generate a new
    counter that starts with one.
    '''
    def __init__(self):
        self._counters = {}
    
    def __getitem__(self, key):
        try:
            return self._counters[key].next()
        except KeyError:
            return self.__missing__(key)
    
    def __contains__(self, key):
        if key in self._counters.keys():
            return True
        else:
            return False

    def __missing__(self, key):
        self._counters[key] = self._id_generator()
        return self._counters[key].next()
    
    def keys(self):
        return self._counters.keys()

    def _id_generator(self):
        counter = 1
        while True:
            yield counter
            counter += 1