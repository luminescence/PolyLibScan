import collections as col

Particle_id = col.namedtuple('Particle_id', 'name, chain, id')

class Particle(object):
    
    def __init__(self, molecule, id_, atom_type, res_id, position):
        self.Id = id_ 
        self.type_ = atom_type
        self.position = position
        self.mol = molecule
        self.residue = res_id
    
    @property
    def residue(self):
        return self._residue

    @residue.setter
    def residue(self, values):
        self._residue = Particle_id._make(values)

    @property
    def mol_name(self):
        return self.mol.name

    @property
    def mol_id(self):
        return self.mol.Id

    def __repr__(self):
        return 'Particle | Id: %d - %s' % (self.Id, self.type_)

class Bond(object):

    def __init__(self, id_, bond_type, members):
        self.Id = id_ 
        self.type_ = bond_type
        self.members = members
        if len(self.members) != 2:
            raise Exception('Bonds must have 2 members')

    def __repr__(self):
        return 'Bond | Id: %d - %s' % (self.Id, self.type_)

class Angle(object):

    def __init__(self, id_, angle_type, members):
        self.Id = id_ 
        self.type_ = angle_type
        self.members = members
        if len(self.members) != 3:
            raise Exception('Angle must have 3 members')

    def __repr__(self):
        return 'Angle | Id: %d - %s' % (self.Id, self.type_)

class Dihedral(object):

    def __init__(self, id_, dihedral_type, members):
        self.Id = id_ 
        self.type_ = dihedral_type
        self.members = members
        if len(self.members) != 4:
            raise Exception('Dihedral must have 4 members')

    def __repr__(self):
        return 'Dihedral | Id: %d - %s' % (self.Id, self.type_)
