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
        
    def combined_name(self, part_name):
        return self.name + '_' + part_name
        
    def create_particles(self):
        part = {}
        for name in self.type_.particles:
            atom_parameters = self.env.atom_type[name]
            id_ = self.env.new_id['particle']
            part[name] = PaI.Particle(self.molecule, id_, atom_parameters, 
                                      (name, 'A', (' ',id_,' ')), 
                                      self.position + atom_parameters.position)
        return part
        
    def make_internal_bond(self):
        bond = []
        for elem1,elem2 in self.type_.bonds:
            bond.append(PaI.Bond(self.env.new_id['bond'], self.env.bond_type[self.name], 
                                 [self.particles[elem1], self.particles[elem2]]))
        return bond
    
    def bind_with(self, monomer):
        if len(monomer.particles) > 1:
            bb_bond = PaI.Bond(self.env.new_id['bond'], self.env.bond_type['polymer'], 
                                      [self.particles['mono_bb'], monomer.particles['mono_bb']])
        else:
            bb_bond = PaI.Bond(self.env.new_id['bond'], self.env.bond_type['polymer'],
                                       [self.particles[self.combined_name('bb')],
                                        monomer.particles[monomer.combined_name('bb')]])
        self.bonds.append(bb_bond)
        monomer.bonds.append(bb_bond)

    def make_internal_angles(self):
        angle = []
        for elem1,elem2,elem3 in self.type_.angles:
            angle.append(PaI.Angle(self.env.new_id['angle'], self.env.angle_type['polymer'], 
                                 [self.particles[elem1], self.particles[elem2], 
                                  self.particles[elem3]]))
        return angle

    def bb_angle_with(self, monomer, monomer2):
        if len(monomer.particles) > 1:
            bb_angle = PaI.Angle(self.env.new_id['angle'], self.env.angle_type['polymer'],
                                [monomer.particles['mono_bb'], self.particles['mono_bb'], 
                                 monomer2.particles['mono_bb']])
            self.angles.append(bb_angle)
            monomer.angles.append(bb_angle)
            monomer2.angles.append(bb_angle)
        else:
            pass

    def angle_with(self, monomer):
    	if len(monomer.particles) > 1:
        	new_angle = PaI.Angle(self.env.new_id['angle'], self.env.angle_type['polymer'], 
                            	[self.particles[self.combined_name('sc')], self.particles['mono_bb'],
                             	monomer.particles['mono_bb']])
        	self.angles.append(new_angle)
        	monomer.angles.append(new_angle)
        else:
        	pass
        
    def make_internal_dihedrals(self):
        dihedral = []
        for elem1,elem2,elem3,elem4 in self.type_.dihedrals:
            dihedral.append(PaI.Dihedral(self.env.new_id['dihedral'], self.env.dihedral_type['polymer'], 
                                 [self.particles[elem1], self.particles[elem2], 
                                  self.particles[elem3],
                                 self.particles[elem4]]))
        return dihedral
        
    def bb_dihedral_with(self, monomer, monomer2, monomer3):
        if len(monomer.particles) > 1:
            bb_dihedral = PaI.Dihedral(self.env.new_id['dihedral'], self.env.dihedral_type['polymer'],
                                [self.particles['mono_bb'], monomer.particles['mono_bb'], 
                                 monomer2.particles['mono_bb'], monomer3.particles['mono_bb']])
        else:
            bb_dihedral = PaI.Dihedral(self.env.new_id['dihedral'], self.env.angle_type['polymer'],
                            [self.particles[self.combined_name('bb')], monomer.particles[monomer.combined_name('bb')], 
                             monomer2.particles[monomer2.combined_name('bb')], monomer3.particles[monomer3.combined_name('bb')]])
        self.dihedrals.append(bb_dihedral)
        monomer.dihedrals.append(bb_dihedral)
        monomer2.dihedrals.append(bb_dihedral)
        monomer3.dihedrals.append(bb_dihedral)
            
    def dihedral_with(self, monomer):
    	if len(monomer.particles) > 1:
        	new_dihedral = PaI.Dihedral(self.env.new_id['dihedral'], self.env.dihedral_type['polymer'], 
                            	[self.particles[self.combined_name('sc')], self.particles['mono_bb'], 
                             	monomer.particles['mono_bb'], monomer.particles[monomer.combined_name('sc')]])
        	self.dihedrals.append(new_dihedral)
        	monomer.dihedrals.append(new_dihedral)
        else:
        	pass