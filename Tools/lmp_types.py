class AtomType(object):
    
    def __init__(self, name, parameters, unique=False):
        self.group_id = 1
        self.name = name
        self.Id = None
        self._mass = 10
        self._radius = 2.0
        self.charge = 0.0
        self.hydrophobicity = 0.0
        self._interacting = True
        self.unique = unique
        self.set_parameters(parameters)

    def set_parameters(self, parameters):
        parameter_names = ['mass', 'radius', 'interacting', 'charge', 'hydrophobicity']
        for key in filter(lambda x:x in parameters.keys(), parameter_names):
            setattr(self, key, parameters[key])

    def mass():
        doc = "The mass property."
        def fget(self):
            return self._mass
        def fset(self, value):
            self._mass = float(value)
        def fdel(self):
            del self._mass
        return locals()
    mass = property(**mass())

    def interacting():
        doc = "The interacting property."
        def fget(self):
            return self._interacting
        def fset(self, value):
            if isinstance(value, basestring):
                if value == 'True':
                    self._interacting = True
                elif value == 'False':
                    self._interacting = False
            elif isinstance(value, bool):
                self._interacting = value
            else:
                raise Exception('interacting parameter must be True or False, in strings or as bool')
        def fdel(self):
            del self._interacting
        return locals()
    interacting = property(**interacting())

    def radius():
        doc = "The radius property."
        def fget(self):
            return self._radius
        def fset(self, value):
            self._radius = float(value)
        def fdel(self):
            del self._radius
        return locals()
    radius = property(**radius())

    def charge():
        doc = "The charge property."
        def fget(self):
            return self._charge
        def fset(self, value):
            self._charge = float(value)
        def fdel(self):
            del self._charge
        return locals()
    charge = property(**charge())

    def __eq__(self, other):
        if self.name == other.name and self.Id == other.Id:
            if self.group_id == other.group_id and self.mass == other.mass:
                return True
            else:
                return False
        else:
            return False

    def __ne__(self, other):
        if self.name != other.name or self.Id != other.Id:
            if self.group_id != other.group_id or self.mass != other.mass:
                return True
            else:
                return False
        else:
            return False

    def __repr__(self):
        return 'AtomType | Name: %s - Atom-Id: %d - Mass: %f' % (self.name, self.Id, self.mass)


class ArchType(object):
    def __init__(self, name, kind, coeffs):
        # self.type_ = {'name': '', 'no': 1}
        self.Id = None
        self.name = name
        self.parameters = {}
        self.kind = kind
        self.parameters['coeffs'] = coeffs

    def __eq__(self, other):
        if self.name == other.name and self.Id == other.Id:
            if self.kind == other.kind and self.parameters['coeffs'] == other.parameters['coeffs']:
                return True
            else:
                return False
        else:
            return False

    def __ne__(self, other):
        if self.name != other.name or self.Id != other.Id:
            if self.kind != other.kind or self.parameters['coeffs'] != other.parameters['coeffs']:
                return True
            else:
                return False
        else:
            return False

class PairType(ArchType):
    def __init__(self, name, kind, rep_only, coeffs, cutoff, single_type_parametrised):
        super(PairType, self).__init__(name, kind, coeffs)
        self.single_type_parametrised = single_type_parametrised
        self.repulsive_only = rep_only
        self.cutoff = cutoff
        
    def __repr__(self):
        repr_str  = 'Pair Type | Name: {} - Pair-Id: {:> 3d} - Kind: {} \nCoefficients: '
        repr_str += '{:> 2f} '*len(self.parameters['coeffs'])
        return  repr_str.format(self.name, self.Id, self.kind, *self.parameters['coeffs'])

class BondType(ArchType):
    def __init__(self, name, kind, coeffs):
        super(BondType, self).__init__(name, kind, coeffs)
        
    def __repr__(self):
    	repr_str  = 'Bond Type | Name: {} - Bond-Id: {:> 3d} - Kind: {} \nCoefficients: ' 
        repr_str += '{:> 2f} '*len(self.parameters['coeffs'])
        return  repr_str.format(self.name, self.Id, self.kind, *self.parameters['coeffs'])

class AngleType(ArchType):
    
    def __init__(self, name, kind, coeffs):
        super(AngleType, self).__init__(name, kind, coeffs)

    def __repr__(self):
        repr_str  = 'Angle Type | Name: {} - Bond-Id: {:> 3d} - Kind: {} \nCoefficients: ' 
        repr_str += '{:> 2f} '*len(self.parameters['coeffs'])
        return  repr_str.format(self.name, self.Id, self.kind, *self.parameters['coeffs'])


class DihedralType(ArchType):
    
    def __init__(self, name, kind, coeffs):
        super(DihedralType, self).__init__(name, kind, coeffs)

    def __repr__(self):
        repr_str  = 'Dihedral Type | Name: {} - Bond-Id: {:> 3d} - Kind: {} \nCoefficients: ' 
        repr_str += '{:> 2f} '*len(self.parameters['coeffs'])
        return  repr_str.format(self.name, self.Id, self.kind, *self.parameters['coeffs'])


class affinity(object):
    '''container for the affinity between 
    two interaction partners
    '''
    def __init__(self, type1, type2, interaction):
        self.type1 = type1
        self.type2 = type2
        self.strength = interaction

    def __repr__(self):
        return '%s | %s | strength: %f'% (self.type1, self.type2, self.strength)