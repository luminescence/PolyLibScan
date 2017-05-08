import numpy as np
import Bio.PDB as PDB
import yaml
from lmp_types import *
from lmp_particlesAndInteractions import *
from lmpObj import LmpObj
from lmp_creator import LmpCreator

class ProteinCreator(LmpCreator):
    '''
    input:
        pdb_path:           string  | path to pdb file 
                                      make sure the filename contains the pdb id.
        config:             string  | path to config file
        mol_id:             int     | id of the molecule
        cg_lvl:             string  | specifies the level of coarse graining
                                      options: full, bb+side, backbone, geometric_center
        with_ghosts:        bool    | if true, ghost particles are added to hold real 
                                      ones in place.
        add_protein_binding:bool    | if true, particles are bound together as 
                                      a polypeptide string
    '''
    def __init__(self, environment, pdb_path, cg_lvl='backbone', with_ghosts=True, add_protein_binding=True):
        super(ProteinCreator, self).__init__(environment)
        self.mol_type = 'protein'
        self.pdb_path = pdb_path
        self.pdb_structure = self._read_pdb(self.pdb_path)
        self.options = {}
        self._cg_lvl = cg_lvl
        self._add_protein_binding = add_protein_binding
        self._with_ghosts = with_ghosts

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
            g_particles, g_bonds = self.add_ghost_particles(molecule)
            molecule.data['particles'] = np.append(
                molecule.data['particles'], g_particles)
            molecule.data['bonds'] += g_bonds
        
        if self.add_protein_binding:
            self.add_protein_properties(molecule)



    ### HELPER FUNCTIONS
    
    def coarse_grain_particles(self, lmp_obj):
        '''
        input
            cg_lvl: string | options: full, bb+side, backbone
        '''
        if self.cg_lvl == 'backbone':
            particles = self._extract_Calphas(lmp_obj)
        elif self.cg_lvl == 'geometric_center':
            particles = self._extract_resi_centered(lmp_obj)
        else:
            raise ValueError('%s-option is not implemented yet.' % self.cg_lvl)
        return particles

    def _read_pdb(self, filename):
        parser = PDB.PDBParser()
        return parser.get_structure('receptor', filename)

    def _extract_resi_centered(self, lmp_obj):
        '''Extract the geometrical centers of all residues
        and places particles in those positions
        '''
        # filter out HETATM
        residues = filter(lambda x:x.get_full_id()[3][0] == ' ', 
                       [resi for resi in self.pdb_structure.get_residues()])
        resi_centers = np.zeros((len(residues),3), dtype=np.dtype('Float64'))
        for i, resi in enumerate(residues):
            atom_coords = np.array([atom.coord for atom in resi.child_list])
            center = atom_coords.mean(axis=0)
            resi_centers[i] = center
        particles = np.empty(len(residues), dtype=object)

        for i, resi in enumerate(residues,0):
            particles[i] = Particle(lmp_obj, self.env.new_id['particle'], lmp_obj.env.atom_type['CA'], 
                                    resi_centers[i].copy())
            particles[i].residue = (residues[i].resname, residues[i].parent.id , residues[i].id)

        return particles

    def _get_config(self, res_type_config):
        with open(res_type_config) as f:
            type_data = yaml.load(f)
        if 'Atoms' in type_data.keys():
            return type_data['Atoms']
        else:
            return type_data

    def _extract_Calphas(self, lmp_obj):
        atoms = filter(lambda x:'CA' in x.fullname, 
                       [atom for atom in self.pdb_structure.get_atoms()])
        particles = np.empty(len(atoms), dtype=object)

        for i, atom in enumerate(atoms,0):
            particles[i] = Particle(lmp_obj, self.env.new_id['particle'], lmp_obj.env.atom_type['CA'], 
                                    atom.coord.copy())
            particles[i].residue = (atom.get_parent().resname, atom.parent.parent.id , atom.parent.id)

        return particles

    def add_ghost_particles(self, lmp_obj):
        '''Adds the same particle to the system but with a different atom type 
        that does not participate in pair interaction. Molecule types are added dynamically.
        '''
        
        # create particles
        
        existing_bonds = len(lmp_obj.data['bonds'])
        g_particles = np.empty(len(lmp_obj.data['particles']), dtype=object)
        g_bonds = []
        for i, real_particle in enumerate(lmp_obj.data['particles']):
            g_particles[i] = Particle(lmp_obj, self.env.new_id['particle'],
                                    lmp_obj.env.atom_type['CA_ghost'],
                                    real_particle.position.copy())
        
            # create bonds
            if 'ghost' not in lmp_obj.env.bond_type:
                raise Exception('add bond type with name ghost to the config file.')

            g_bonds += [Bond(self.env.new_id['bond'],
                           lmp_obj.env.bond_type['ghost'],
                           [g_particles[i], real_particle])]
        return g_particles, g_bonds


    def add_protein_properties(self, protein):
        '''Adds the bonds, angles and dihedrals to the particles that
        are typical for the protein, i.e. its a chain.
        '''
        particle_no = len(protein.data['particles'])
        bonds = []
        angles = []    
        dihedrals = []
        for i, particle in enumerate(protein.data['particles']):
            if i+1 < particle_no:
                if particle.mol_id == protein.data['particles'][i+1].mol_id:
                    bonds.append(Bond(self.env.new_id['bond'], protein.env.bond_type['peptide'], 
                                    [particle, protein.data['particles'][i+1]] ) )
                    
            if i+2 < particle_no:
                if particle.mol_id == protein.data['particles'][i+1].mol_id and particle.mol_id == protein.data['particles'][i+2].mol_id:
                    angles.append( Angle( self.env.new_id['angle'], protein.env.angle_type['peptide'], 
                                        [ particle, protein.data['particles'][i+1], protein.data['particles'][i+2]] ) )
            if i+3 < particle_no:
                if particle.mol_id == protein.data['particles'][i+1].mol_id and particle.mol_id == protein.data['particles'][i+2].mol_id and particle.mol_id == protein.data['particles'][i+3].mol_id:
                    dihedrals.append(Dihedral(self.env.new_id['dihedral'], protein.env.dihedral_type['peptide'], 
                                            [ particle, protein.data['particles'][i+1], protein.data['particles'][i+2], protein.data['particles'][i+3]]))
        
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

        real_particle = np.extract([x.type_.name !='CA_ghost' for x in self.data['particles']], self.data['particles'])
        for i, particle in enumerate(real_particle):
            type_id = max_id + 1 + i
            #print type_id
            type_name = '%s|%s%d|%d' % (particle.residue[0], particle.residue[1], particle.residue[2][1], type_id)
            type_data[particle.residue[0]]['name'] = type_name
            self.env.atom_type[type_name] = AtomType(type_id, 
                                                      type_data[particle.residue[0]])
            # setting type for the particle
            real_particle[i].type_ = self.types['particles'][type_name]

    def change_to_res_based(self, molecule, res_type_config=None):
        if res_type_config:
            type_data = self._get_config(res_type_config)
            ## add types

            for name, config_type in type_data.items():
                molecule.env.atom_type.define_type(name, AtomType(name, config_type))
        ## read particle res and change restype
        real_particle = np.extract([x.type_.name !='CA_ghost' for x in molecule.data['particles']], molecule.data['particles'])
        for i, particle in enumerate(real_particle):
            if particle.residue[0] in molecule.env.atom_type:
                particle.type_ = molecule.env.atom_type[particle.residue[0]]
            else:
                raise ValueError('particle type %s not found in types.' % particle.residue[0])

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
