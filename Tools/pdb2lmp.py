import itertools as it

import Bio.PDB as PDB
import PolyLibScan.Database.db as DB
import numpy as np
import yaml

from lmp_creator import LmpCreator
from lmp_particlesAndInteractions import *


class protonation_methods_bundled(particle_methods_bundled):
    def add_protonation(self, molecule, protonation_data, ph=8.0):
        for res_name, atom_id,pka in it.izip(protonation_data['resn'],
                                             protonation_data[['chain', 'resi', 'iCode']],
                                             protonation_data['pka']):
            protonation_change = self._calc_protonation_change(res_name, ph, pka)
            if protonation_change != 0:
                particle = self._find_particle_by_pdb_id(atom_id, molecule)
                self._change_protonation(particle, protonation_change)

    @staticmethod
    def _calc_protonation_change(res_name, ph, pka):
        # Residues that get protonated below their pKa
        if res_name in ['N+', 'GLU', 'ASP', 'HIS']:
            change = (1,0)
        # Residues that get deprotonated above their pKa
        elif res_name in ['LYS', 'ARG', 'TYR', 'CYS', 'C-']:
            change = (0,-1)
        # Weird case
        else:
            raise ValueError("for some reason you want to change the protonation of %s" % res_name)
        if ph < pka:
            return change[0]
        elif ph > pka:
            return change[1]
        else:
            return 0

    def _change_protonation(self, particle, protonation_change):
        particle.charge = particle.type_.charge + protonation_change


class surface_methods_bundled(particle_methods_bundled):
    def add_surface_energy(self, atom_id, surface_energy, molecule):
        particle = self._find_particle_by_pdb_id(atom_id, molecule)
        if not particle.type_.unique:
            protein_particle_type = self._make_particle_unique(particle)
        particle.type_.surface_energy = surface_energy

    def add_surface_energies(self, molecule, surface_data, area_energy=False):
        if area_energy:
            e_id = 'areaParameter'
        else:
            #singleParameter column
            e_id = 'singleParameter'

        for atom_id,surface_energy in it.izip(surface_data[['chain', 'id', 'iCode']], surface_data[e_id]):
            self.add_surface_energy(atom_id, surface_energy, molecule)


class ProteinCreator(LmpCreator, protonation_methods_bundled, surface_methods_bundled):
    '''
    input:
        pdb_path:            string  | path to pdb file 
                                       make sure the filename contains the pdb id.
        config:              string  | path to config file
        mol_id:              int     | id of the molecule
        cg_lvl:              string  | specifies the level of coarse graining
                                       options: full, bb+side, backbone, geometric_center
        with_ions:           bool    | if true, ions are added to the protein and kept in 
                                       place like all protein particles.
        with_ghosts:         bool    | if true, ghost particles are added to hold real 
                                       ones in place.
        add_protein_binding: bool    | if true, particles are bound together as 
                                       a polypeptide string
    '''
    def __init__(self, environment, pdb_path, cg_lvl='backbone', 
                 with_ions=False, with_ghosts=True, add_protein_binding=True,
                 surface_file='', protonation_file='', 
                 use_amino_acids=True, ph=8.0):
        super(ProteinCreator, self).__init__(environment)
        self.mol_type = 'protein'
        self.pdb_path = pdb_path
        self.pdb_structure = self._read_pdb(self.pdb_path)
        self.options = {}
        self._cg_lvl = cg_lvl
        self._add_protein_binding = add_protein_binding
        self._with_ions = with_ions
        self._with_ghosts = with_ghosts
        self._with_amino_acids = use_amino_acids
        self.surface_file = surface_file
        self.protonation_file = protonation_file
        self.ph = ph

    def __str__(self):
        string  = ['Molecule Factory creates: %s' % (self.mol_type.upper())]
        string += ['Atom style:  %s' % self.globals['atom_style']]
        string += ['Bond style:  %s' % self.globals['bond_style']]
        string += ['Angle style: %s' % self.globals['angle_style']]
        string += ['Pdb Structure: %s' % self.pdb_path]
        string += ['Coarse Grain Level: %s' % self._cg_lvl.upper()]
        return '\n'.join(string)

    ### CREATE FUNCTIONS
    def _create_config(self, molecule):
        molecule.mol_type = self.mol_type

    # _create_types() is implemented in parent class

    def _create_particles_and_interactions(self, molecule):
        molecule.data['particles'] = self.coarse_grain_particles(molecule)
        molecule.pdb_structure = self.pdb_structure
        molecule.pdb_id = self.pdb_path.split('/')[-1].split('.')[0] # this has to be changed in the long run.

    def _create_custom(self, molecule):
        if self.with_ghosts:
            self.add_ghost_particles(molecule)
        if self._with_amino_acids:
            self.change_to_res_based(molecule)
        if self.surface_file:
            surface_data = self.get_surface_data(self.surface_file, molecule.pdb_id.upper())
            self.add_surface_energies(molecule, surface_data)
        if self.protonation_file:
            protonation_data = self.get_protonation_data(self.protonation_file, molecule.pdb_id.upper())
            self.add_protonation(molecule, protonation_data, self.ph)
            
        
        if self.add_protein_binding:
            self.add_protein_properties(molecule)

    ### HELPER FUNCTIONS
    def get_surface_data(self, path, pdb_id):
        surface_db = DB.Database(path)
        return surface_db._load_table('/', pdb_id)

    def get_protonation_data(self, path, pdb_id):
        surface_db = DB.Database(path)
        return surface_db._load_table('/protonation', pdb_id)

    def coarse_grain_particles(self, lmp_obj):
        '''
        input
            cg_lvl: string | options: full, bb+side, backbone
        '''
        if self.cg_lvl == 'backbone':
            particle_data_gen = self._extract_Calphas()
        elif self.cg_lvl == 'geometric_center':
            particle_data_gen = self._extract_resi_centered()
	elif self.cg_lvl == 'bb+sc':
	    particle_data_gen = self._extract_bb_sc()
        else:
            raise ValueError('%s-option is not implemented yet.' % self.cg_lvl)
        particles = self._add_particles_to_molecule(lmp_obj, particle_data_gen)
        
        if self._with_ions:
            ions_gen = self._extract_ions()
            ions = self._add_particles_to_molecule(lmp_obj, ions_gen)
            particles = np.concatenate((particles, ions))

        return particles

    def _extract_ions(self):
        filter_non_aa = lambda x: x.get_full_id()[3][0]!=' '
        filter_out_water = lambda x: x.get_full_id()[3][0]!='W'
        ions = filter(lambda x: filter_non_aa(x) and filter_out_water(x), 
                      self.pdb_structure.get_atoms())

        for i, atom in enumerate(ions):
            yield atom.coord.copy(), atom.id, atom.parent.parent.id , atom.parent.id
            
    def _read_pdb(self, filename):
        parser = PDB.PDBParser()
        return parser.get_structure('receptor', filename)

    def _extract_resi_centered(self):
        '''Extract the geometrical centers of all residues
        and places particles in those positions
        '''
        # filter out HETATM
        filter_out_non_aa = lambda x: x.id[0]==' '
        residues = filter(filter_out_non_aa, self.pdb_structure.get_residues())
        resi_centers = np.zeros((len(residues),3), dtype=np.dtype('Float64'))
        for i, resi in enumerate(residues):
            atom_coords = np.array([atom.coord for atom in resi.child_list])
            center = atom_coords.mean(axis=0)
            resi_centers[i] = center
            yield center, residues[i].resname, residues[i].parent.id , residues[i].id

    def _add_particles_to_molecule(self, lmp_obj, particle_data):
        particles = []
        for i, p in enumerate(particle_data):
            particles.append(Particle(lmp_obj, self.env.new_id['particle'], 
                                      lmp_obj.env.atom_type['BB_bb'], 
                                      (p[1], p[2], p[3]),
                                      p[0].copy()))
        return particles        

    def _get_config(self, res_type_config):
        with open(res_type_config) as f:
            type_data = yaml.load(f)
        if 'Atoms' in type_data.keys():
            return type_data['Atoms']
        else:
            return type_data

    def _extract_bb_sc(self):
        '''Extract position of backbone and sidechain
	and create generate particle-information 
	of yield backbone and geometrically centered 
	sidechain''' 
	residues = self.pdb_structure.get_residues()
	
	for resi in residues:
            atoms = resi.child_dict
            c_alpha = atoms['CA']
            yield c_alpha.coord.copy(), resi.resname, resi.parent.id , resi.id

	    side_chain = [atom for name,atom in atoms.items() if not name in ['CA','N', 'C', 'O']] 
	    if not len(side_chain) == 0:
 	        sc_center = np.array([atom.coord for atom in side_chain]).mean(axis=0)
	        raise NotImplemented('sidechain particle needs unique signiture, but has same as bb.')
                yield sc_center, resi.resname, resi.parent.id, resi.id

    def _extract_Calphas(self):
        filter_c_alphas = lambda x: x.id=='CA' and x.element=='C'
        atoms = filter(filter_c_alphas, self.pdb_structure.get_atoms())

        for i, atom in enumerate(atoms,0):
            yield atom.coord.copy(), atom.get_parent().resname, atom.parent.parent.id , atom.parent.id

    def add_ghost_particles(self, lmp_obj):
        '''Adds the same particle to the system but with a different atom type 
        that does not participate in pair interaction. Molecule types are added dynamically.
        '''
        
        # create particles
        
        existing_bonds = len(lmp_obj.data['bonds'])
        g_particles = np.empty(len(lmp_obj.data['particles']), dtype=object)
        g_bonds = []
        for i, real_particle in enumerate(lmp_obj.data['particles']):
            res_id = ('ghost', real_particle.residue.chain, ('ghost', real_particle.residue.id[1], real_particle.residue.id[2]))
            g_particles[i] = Particle(lmp_obj, self.env.new_id['particle'],
                                    lmp_obj.env.atom_type['BB_ghost_bb'],
                                    res_id,
                                    real_particle.position.copy())
        
            # create bonds
            if 'ghost' not in lmp_obj.env.bond_type:
                raise Exception('add bond type with name ghost to the config file.')

            g_bonds += [Bond(self.env.new_id['bond'],
                           lmp_obj.env.bond_type['ghost'],
                           [g_particles[i], real_particle])]
        lmp_obj.data['particles'] = np.append(
                lmp_obj.data['particles'], g_particles)
        lmp_obj.data['bonds'] += g_bonds

    def add_protein_properties(self, protein):
        '''Adds the bonds, angles and dihedrals to the particles that
        are typical for the protein, i.e. it's a chain.
        '''
        particle_no = len(protein.data['particles'])
        bonds = []
        angles = []
        dihedrals = []
        filter_non_aa = lambda x: x.residue.id[0]==' '
        amino_acids = filter(filter_non_aa, protein.data['particles'])

        for particle1,particle2 in it.izip(amino_acids[:-1], amino_acids[1:]):
                if particle1.mol_id == particle2.mol_id:
                    bonds.append(Bond(self.env.new_id['bond'], protein.env.bond_type['peptide'], 
                                 [particle1, particle2]))
        for particle1,particle2,particle3 in it.izip(amino_acids[:-2], amino_acids[1:-1], amino_acids[2:]):
                if particle1.mol_id == particle2.mol_id and particle2.mol_id == particle3.mol_id:
                    angles.append(Angle(self.env.new_id['angle'], protein.env.angle_type['peptide'], 
                                         [particle1, particle2, particle3]))
        # for particle1,particle2,particle3,particle4 in it.izip(amino_acids[:-3], amino_acids[1:-2], amino_acids[2:-1], amino_acids[3:]):
        #         if particle1.mol_id == particle2.mol_id and particle1.mol_id == particle3.mol_id and particle1.mol_id == particle4.mol_id:
        #             dihedrals.append(Dihedral(self.env.new_id['dihedral'], protein.env.dihedral_type['peptide'], 
        #                                     [particle1, particle2, particle3, particle4]))
        
        protein.data['bonds']     += bonds
        protein.data['angles']    +=  angles
        protein.data['dihedrals'] +=  dihedrals

    def change_to_seqIdBased(self, res_config):
        '''change all amino-acid beads to be of sequence id based
        type.
        '''
        type_data = self.get_config(res_config)
        type_data = {typ['name']: typ for typ in type_data.values()}
        ## add types
        #get max id
        max_id = max([a_type.Id for a_type in self.types['particles'].values()])

        real_particle = np.extract([x.type_.name !='BB_ghost' for x in self.data['particles']], self.data['particles'])
        for i, particle in enumerate(real_particle):
            type_id = max_id + 1 + i
            #print type_id
            type_name = '%s|%s%d|%d' % (particle.residue.name, particle.residue.chain, particle.residue.id[1], type_id)
            type_data[particle.residue.name]['name'] = type_name
            self.env.atom_type[type_name] = AtomType(type_id, type_data[particle.residue.name])
            # setting type for the particle
            real_particle[i].type_ = self.types['particles'][type_name]

    def change_to_res_based(self, molecule, res_type_config=None):
        if res_type_config:
            type_data = self._get_config(res_type_config)
            ## add types
            for name, config_type in type_data.items():
                molecule.env.atom_type.define_type(name, AtomType(name, config_type))
        ## read particle res and change restype
        real_particle = np.extract([x.type_.name !='BB_ghost_bb' for x in molecule.data['particles']], molecule.data['particles'])
        for i, particle in enumerate(real_particle):
            if particle.residue.name + '_bb' in molecule.env.atom_type:
                if particle.type_.unique:
                    # in case the particle type is already unique, the resi
                    # specific properties are copied over.
                    res_based_type = molecule.env.atom_type[particle.residue.name]
                    unique_particle_type = particle.type_
                    properties = ['mass', 'radius', 'charge', 'hydrophobicity']
                    for prop in properties:
                        setattr(unique_particle_type, prop, getattr(res_based_type, prop))
                else:
                    new_type = molecule.env.atom_type[particle.residue.name + '_bb']
                    particle.type_ = new_type
                    particle.charge = new_type.charge
            else:
                raise ValueError('particle type %s not found in types.' % particle.residue.name)

    def cg_lvl():
        doc = "The cg_lvl property."
        def fget(self):
            return self._cg_lvl
        def fset(self, value):
            if value in ['full', 'bb+sc', 'backbone', 'geometric_center']:
                self._cg_lvl = value
            else:
                raise ValueError('%s option not possible.'% value)
        return locals()
    cg_lvl = property(**cg_lvl())


    def with_ghosts():
        doc = "The with_ghosts property."
        def fget(self):
            return self._with_ghosts
        def fset(self, value):
            if value in [True, False]:
                self._with_ghosts = value
            else:
                raise ValueError('attribute must be True or False.')
        return locals()
    with_ghosts = property(**with_ghosts())

    def add_protein_binding():
        doc = "The add_protein_binding property."
        def fget(self):
            return self._add_protein_binding
        def fset(self, value):
            if value in [True, False]:
                self._add_protein_binding = value
            else:
                raise ValueError('attribute must be True or False.')
        return locals()
    add_protein_binding = property(**add_protein_binding())
