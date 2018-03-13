import numpy as np
import lmp_particlesAndInteractions as PaI

class Monomer(object):
    def __init__(self, position, name, type_, env, molecule):
        self.name = name
        self.position = np.array(position)
        self.type_ = type_
        self.env = env
        self.molecule = molecule
        self.particles = self.create_particles()
        self.bonds = self.make_internal_bond()
        self.angles = self.make_internal_angles()
        self.dihedrals = self.make_internal_dihedrals()

    def __repr__(self):
        return '%s monomer | Position %s | Molecule: %s' % (
                            self.type_.name, 
                            self.position, 
                            self.molecule.mol_info())
        
    def combined_name(self, part_name):
        return self.name + '_' + part_name
        
    def create_particles(self):
        part = {}
        for mono_id, name in self.type_.particles.iteritems():
            atom_parameters = self.env.atom_type[name]
            id_ = self.env.new_id['particle']
            part[mono_id] = PaI.Particle(self.molecule, id_, atom_parameters, 
                                      (name, 'A', (' ',id_,' ')), 
                                      self.position + self.type_.positions[mono_id])
            self.molecule.data['particles'].append(part[mono_id])
        return part
        
    def make_internal_bond(self):
        bonds = []
        for bond_type,elem1,elem2 in self.type_.bonds:
            bonds.append(PaI.Bond(self.env.new_id['bond'], self.env.bond_type[bond_type], 
                                 [self.particles[elem1], self.particles[elem2]]))
        self.molecule.data['bonds'] += bonds
        return bonds
    
    def bind_with(self, monomer):
        for bond_type,p1,p2 in self.molecule.type.bonds:
            bond = PaI.Bond(self.env.new_id['bond'], self.env.bond_type[bond_type], 
                                  [self.particles[p1], monomer.particles[p2]])
            self.bonds.append(bond)
            monomer.bonds.append(bond)
            self.molecule.data['bonds'].append(bond)

    def make_internal_angles(self):
        angle = []
        for angle_type,elem1,elem2,elem3 in self.type_.angles:
            angle.append(PaI.Angle(self.env.new_id['angle'], self.env.angle_type[angle_type], 
                                 [self.particles[elem1], self.particles[elem2], 
                                  self.particles[elem3]]))
        return angle

    def angles_with_one(self, partner, parameters):
        for type_name,own,other in parameters:
            particles = [self.particles[id_] for id_ in own]
            particles+= [partner.particles[id_] for id_ in other]

            angle = PaI.Angle(self.env.new_id['angle'], self.env.angle_type[type_name], 
                                [particles[0], particles[1],particles[2]])
            self.angles.append(angle)
            partner.angles.append(angle)
            self.molecule.data['angles'].append(angle)

    def angles_with_two(self, prev_mono, next_mono, parameters):
        for type_name,prev,own,next_ in parameters:
            particles = [prev_mono.particles[id_] for id_ in prev]
            particles+= [self.particles[id_] for id_ in own]
            particles+= [next_mono.particles[id_] for id_ in next_]

            angle = PaI.Angle(self.env.new_id['angle'], self.env.angle_type[type_name], 
                                [particles[0], particles[1],particles[2]])
            self.angles.append(angle)
            prev_mono.angles.append(angle)
            next_mono.angles.append(angle)
            self.molecule.data['angles'].append(angle)
        
    def make_internal_dihedrals(self):
        dihedral = []
        for elem1,elem2,elem3,elem4 in self.type_.dihedrals:
            dihedral.append(PaI.Dihedral(self.env.new_id['dihedral'], self.env.dihedral_type['polymer'], 
                                 [self.particles[elem1], self.particles[elem2], 
                                  self.particles[elem3],
                                 self.particles[elem4]]))
        return dihedral
            
    def dihedral_with_two(self, partner, parameters):
        for type_name,own,other in parameters:
            particles = [self.particles[id_] for id_ in own]
            particles+= [partner.particles[id_] for id_ in other]

            dihedral = PaI.Dihedral(self.env.new_id['dihedral'], 
                                       self.env.dihedral_type[type_name],
                                       particles)
            self.dihedrals.append(dihedral)
            partner.dihedrals.append(dihedral)
            self.molecule.data['dihedrals'].append(dihedral)
            
    def dihedral_with_three(self, neighbor, next_, parameters):
        for type_name,own_p_ids,neighbor_p_ids,next_neighbor_p_ids in parameters:
            particles = [self.particles[id_] for id_ in own_p_ids]
            particles+= [neighbor.particles[id_] for id_ in neighbor_p_ids]
            particles+= [next_.particles[id_] for id_ in next_neighbor_p_ids]

            dihedral = PaI.Dihedral(self.env.new_id['dihedral'], 
                                       self.env.dihedral_type[type_name],
                                       particles)
            self.dihedrals.append(dihedral)
            neighbor.dihedrals.append(dihedral)
            next_.dihedrals.append(dihedral)
            self.molecule.data['dihedrals'].append(dihedral)
    
    def dihedral_with_four(self, neighbor, next_, next_next, parameters):
        for type_name,own_p_ids,neighbor_p_ids,next_neighbor_p_ids,next_next_p_ids in parameters:
            particles = [self.particles[id_] for id_ in own_p_ids]
            particles+= [neighbor.particles[id_] for id_ in neighbor_p_ids]
            particles+= [next_.particles[id_] for id_ in next_neighbor_p_ids]
            particles+= [next_next.particles[id_] for id_ in next_next_p_ids]

            dihedral = PaI.Dihedral(self.env.new_id['dihedral'], 
                                       self.env.dihedral_type[type_name],
                                       particles)
            self.dihedrals.append(dihedral)
            neighbor.dihedrals.append(dihedral)
            next_.dihedrals.append(dihedral)
            next_next.dihedrals.append(dihedral)
            self.molecule.data['dihedrals'].append(dihedral)