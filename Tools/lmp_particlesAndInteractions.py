import collections as col
import warnings
from lmp_types import AtomType

Particle_id = col.namedtuple('Particle_id', 'name, chain, id')

class LmpObj(object):
    def __init__(self, id_, type_):
        self.Id = id_ 
        self.type_ = type_

    def __eq__(self, other):
        return id(self) == id(other)

    def __hash__(self):
        return id(self)


class Particle(LmpObj):
    
    def __init__(self, molecule, id_, atom_type, res_id, position):
        super(Particle, self).__init__(id_, atom_type)
        self.position = position
        self.mol = molecule
        self.residue = res_id
        self.charge = self.type_.charge

    @property
    def type_(self):
        return self._type_

    @type_.setter
    def type_(self, new_type):
        """raise Warning if the atom type is truly changed, i.e. if it's not just about making the particle unique or changing from backbone to amino acid"""

        particle_already_had_type_ = hasattr(self, 'type_')
        if particle_already_had_type_:
            allowed_type_changes = ['BB_bb', new_type.name, new_type.name.split('|')[0]]
            if self.type_.name not in allowed_type_changes:
                warnings.warn('Trying to change type, this would require redoing the protonation and checking the charge! Implement it!')

        self._type_ = new_type


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
        return 'Particle | Id: %d - %s | Charge: %d' % (self.Id, self.type_, self.charge)


class Interaction(LmpObj):
    def __init__(self, id_, type_, members):
        super(Interaction, self).__init__(id_, type_)
        self.members = members


class Bond(Interaction):

    def __init__(self, id_, bond_type, members):
        super(Bond, self).__init__(id_, bond_type, members)
        if len(self.members) != 2:
            raise Exception('Bonds must have 2 members')

    def __repr__(self):
        return 'Bond | Id: %d - %s' % (self.Id, self.type_)


class Angle(Interaction):

    def __init__(self, id_, angle_type, members):
        super(Angle, self).__init__(id_, angle_type, members)
        if len(self.members) != 3:
            raise Exception('Angle must have 3 members')

    def __repr__(self):
        return 'Angle | Id: %d - %s' % (self.Id, self.type_)


class Dihedral(Interaction):

    def __init__(self, id_, dihedral_type, members):
        super(Dihedral, self).__init__(id_, dihedral_type, members)
        if len(self.members) != 4:
            raise Exception('Dihedral must have 4 members')

    def __repr__(self):
        return 'Dihedral | Id: %d - %s' % (self.Id, self.type_)


class particle_methods_bundled:
    def _find_particle_by_pdb_id(self, pdb_residue_id, molecule, consider_ghost_atoms=False):
        '''Finds the id of the particle in the molecule object based on the
        pdb internal id.
        '''
        #reformat pdb_residue_id
        pdb_residue_id = list(pdb_residue_id)
        pdb_residue_chain = pdb_residue_id.pop(0)
        pdb_residue_id = [' '] + pdb_residue_id # 'real' particles have ' ' in the first position, 'ghost' atoms have 'ghost'
        pdb_residue_id = tuple(pdb_residue_id)

        if consider_ghost_atoms:
            NotImplementedError('Finding ghost particles by Id is not implemented yet!')

        def has_right_id(search_particle):
            return (search_particle.residue != None
                    and search_particle.residue.chain== pdb_residue_chain  # chain
                    and search_particle.residue.id == pdb_residue_id)   # all other information

        particle = filter(has_right_id, molecule.data['particles'])
        if len(particle)>1:
            raise Exception('Found more than one particle: %s \nCheck your pdb file for duplicate entries.' % (
                [p.residue for p in particle]))
        elif len(particle)==0:
            raise Exception("There is no particle with chain %s, id %d, iCode '%s' and ghost status %s" % (
                                    pdb_residue_chain, 
                                    pdb_residue_id[1], 
                                    pdb_residue_id[2], 
                                    pdb_residue_id[0]))
        return particle[0]

    def _make_particle_unique(self, particle):
        '''assigns new, unique particle-type to given particle.
        The new particle type has the same properties as the old one,
        but is uniquely assigned to just this particle.
        Returns the newly created atom type.

        Input:
            particle:  [Particle Obj]

        Output:
            [atom_type]
        '''

        new_type_name = '%s|%s|%d|%d' % (particle.residue.name, particle.residue.chain,
                                         particle.residue.id[1], particle.type_.Id)
        old_type = self.env.atom_type[particle.type_.name]
        self.env.atom_type[new_type_name] = AtomType(new_type_name,
                                                 {'mass': old_type.mass,
                                                  'radius': old_type.radius,
                                                  'charge': old_type.charge,
                                                  'hydrophobicity': old_type.hydrophobicity,
                                                  'surface_energy': old_type.surface_energy,
                                                  'interacting': old_type.interacting},
                                                  unique=True)

        # change particle type to newly created
        particle.type_ = self.env.atom_type[new_type_name]
        return self.env.atom_type[new_type_name]
