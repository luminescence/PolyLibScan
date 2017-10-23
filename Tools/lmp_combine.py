from lmp_types import *
from lmp_particlesAndInteractions import *
from lmp_creator import LmpCreator
from lmpObj import LmpObj
import lmp_helpers as helpers
import numpy as np
import itertools as it
import sets
import pandas as pd
import PolyLibScan.Database.db as db


class EnvManipulator(object, particle_methods_bundled):
    '''LmpObj constructor is not used here!!
    '''
    def __init__(self, environment, randomize_positions=True, auto_repulsion=True, 
                 include_affinities=False, affinity_file=None):
        self.env = environment
        self.molecules = self.env.molecules
        self.option = {}
        self.option['random_positions'] = randomize_positions
        self.option['include_affinities'] = include_affinities
        self.option['auto_repulsion'] = auto_repulsion
        self.affinity_file = affinity_file
        self.repulsion_value = -50
        
        self.distribute_molecules()

        self.randomize_positions = randomize_positions
        self._create_custom()

    def _list_all_monomers(self):
        polymers = filter(lambda x:x.mol_type == 'polymer', 
                 self.molecules.values())
        monomers_used = sets.Set([])
        for poly in polymers:
            mono_set = sets.Set([mono.type_.name for mono in poly.data['particles']])
            monomers_used = monomers_used.union(mono_set)
        return list(monomers_used)

    def _create_custom(self):
        if self.option['random_positions']:
            self.create_random_start_positions()
        if self.option['auto_repulsion']:
            for mol in self.molecules.values():
                if mol.mol_type == 'polymer':
                    self.add_auto_repulsion(mol, self.repulsion_value)
        if self.option['include_affinities']:
            #proteins = filter(lambda x:x.mol_type == 'protein', self.molecules)
            #pdb_ids = [mol.pdb_id for mol in proteins]
            self.add_affinities()

    ### HELPER FUNCTIONS

    ## RANDOMIZE STARTING POSITIONS
    def create_random_start_positions(self):
        '''randomly distributes all polymers around the first 
        protein in the list.
        '''
        i = 0
        # pick first protein molecule
        protein = filter(lambda x:x.mol_type=='protein', self.molecules.values())[0]
        for polymer in filter(lambda x:x.mol_type=='polymer', self.molecules.values()):
            shift = self._random_shift(polymer, protein)
            polymer.move(shift)
            axis, angle = self._random_rotation()
            polymer.rotate(axis, angle)
            reduced_mol_list = filter(lambda x:x.name != polymer.name, self.molecules.values())
            while self.overlapping_with_any(polymer, reduced_mol_list):
                shift = self._random_shift(polymer, protein)
                polymer.move(shift)
                axis, angle = self._random_rotation()
                polymer.rotate(axis, angle)
        self.env.calc_box(add_margin=True)

    def _random_shift(self, shift_mol, center_mol):
        '''shifts a box with the size of the shift_mol 
        around the center_mol.

        input:
            shift_mol : molecule that shifts [LmpObj Type]
            center_mol: molecule around the shift occurs[LmpObj Type]

        output:
            deltas: np.ndarray
            
        First, the shift is done by the particles' hypothetical box.
        The box is subsequently used to measure and see if there is overlap to existing 
        boxes.
        '''
        c_mol_radius_c = (center_mol.box[:,1] - center_mol.box[:,0])/2.0
        # polymer_vectors
        s_mol_radius_c = (shift_mol.box[:,1] - shift_mol.box[:,0])/2.0
        # transform to spherical coordinates
        c_sph_radius_s = helpers.Spherical_np(c_mol_radius_c)[0]
        s_sph_radius_s = helpers.Spherical_np(s_mol_radius_c)[0]
        # create suitable new random coordinates
        vector_center_new_s = np.concatenate((np.random.uniform(0.0, s_sph_radius_s[0]/2, 1) + c_sph_radius_s[0],
                                              np.random.uniform(0.0, 2 * np.pi, 1),
                                              np.random.uniform(0.0,     np.pi, 1)) )
        new_shift_center_c = helpers.Cartesian_np(vector_center_new_s)[0] + center_mol.center()
        # technically, we calculate the difference of the centers here and will use 
        # that later on for the transformation of the particles. since these are only the 
        # differences it is still correct.
        delta_center = new_shift_center_c - shift_mol.center()
        return delta_center

    def _random_rotation(self):
        '''creates random values for rotation.
        '''
        rot_angle =  np.random.uniform(0.0, 2 * np.pi, 1)[0]
        rot_axis = np.concatenate((np.array([1.0]),
                                    np.random.uniform(0.0, 2 * np.pi, 1),
                                    np.random.uniform(0.0,     np.pi, 1)) )
        return helpers.Cartesian_np(rot_axis)[0], rot_angle
        # molecule.rotate(helpers.Cartesian_np(rot_axis)[0], rot_angle)
        
            
    def distribute_molecules(self):
        '''makes sure no molecules are overlapping.
        '''
        molecules = self.molecules.items()
        
        self.env.box = molecules[0][1].box.copy()
        for name, mol in molecules[1:]:
            if self.overlapping(mol.data['particles'], 
                [p for p in mol.data['particles'] 
                    for mol in self.env.molecules.values()]):
                self.shift_mol(mol, self.env.box)
            self.env.box = self.add_to_box(mol, self.env)

    def add_to_box(self, mol, environment):
        '''Extends the global box by including the box given in the argument.
        If it was done differently, all the molecules would be included at once and 
        the repositioning when overlapping would include the molecule itself.
        Still, this should be changed.
        '''
        big_box = np.zeros([3,2])
        for dim, dim_interval in enumerate(mol.box):
            big_box[dim,0] = min(environment.box[dim,0], dim_interval[0])
            big_box[dim,1] = max(environment.box[dim,1], dim_interval[1])
        return big_box

        
    def shift_mol(self, mol, box, margin=5.0):
        '''
        '''
        mol_box = mol.box.copy()
        global_box_dim = np.array([dim[1]-dim[0] for dim in box])
        small_dim = global_box_dim.argsort()[0]
        shift = np.zeros(3)
        shift[small_dim] = box[small_dim,1] - mol_box[small_dim,0] + margin
        # print 'Shifting Molecule %d by %f in dimension %d'% (mol_number, shift[small_dim], small_dim)
        mol.move(shift)
        
    def overlapping_with_any(self, mol1, other_molecules):
        for mol2 in other_molecules:
            if self.overlapping(mol1.data['particles'], 
                                mol2.data['particles']):
                return True
        return False

    def overlapping(self, mol1_particles, mol2_particles, margin=5):
        '''Determines, if two boxes overlap

        input:
            mol1_particles:   lmp_obj
            mol2_particles:   lmp_obj
            margin:     integer [Angstrom]

        The margin should make sure that no atoms come so close 
        as to blow up the system.
        '''
        mol1_coords = np.array([p.position for p in mol1_particles])
        mol2_coords = np.array([p.position for p in mol2_particles])

        for c1 in mol1_coords:
            for c2 in mol2_coords:
                if np.linalg.norm((c1-c2)) < margin:
                    return True
        return False

    ### FEATURE: ADD AFFINITIES


    def activeSiteParticles(self, protein, active_site_path):
        '''returns the information of the residues of the active site.
        
        input:
            protein: lammps molecule [lmpObj]
            pdb_id: pdb id of the protein [string]
            active_site_database: pandas database [pandas.io.pytables.HDFStore]
        
        return numpy array
        
        From the Information of the database, the corresponding particles 
        are identified, and the particle id, the chain Id and the residues
        id of the pdb file are returned as a numpy array. This format makes it 
        easy to write to a HDF5 table.
        '''
        pdb_id = protein.pdb_id
        active_site_database = pd.HDFStore(active_site_path)
        AS_residues = active_site_database['/%s' % pdb_id.upper()]
        active_site_database.close()
        particle_ids = [self._find_particle_by_pdb_id((resi.chain, resi.ID, resi.iCode), protein).Id
                         for resi in AS_residues.itertuples()]
        results = np.array(zip(particle_ids, AS_residues['chain'], AS_residues['ID'], AS_residues['iCode']), 
                 dtype=[('xyz', '<i2'), ('chain', '|S1'), ('pdb_id', '<i2'), ('iCode', '|S1')])

        return results

    def add_affinities(self, affinity_file=None):
        '''reads in all docking and epitopsy data.
        The files have to end with either, '.pdb' or '.dx'
        in order to be used correctly.

        input:
            affinity_file:      String      path to hdf5 formatted file 
                                            (created with pandas.HDFStore)
        '''

        # checking input
        if affinity_file:
            pass
        elif self.affinity_file:
            affinity_file = self.affinity_file
        else:
            affinity_file = self.env.globals['affinity_file']
        store = db.Database(affinity_file, mode='r')

        
        pdbIds = self.get_protein_ids(self.molecules)
    
        monomer_names = self._list_all_monomers()

        for protein_id, monomer in it.product(pdbIds, monomer_names):
            affinities = self.get_affinities(store, protein_id, monomer, 50)
            self.add_interaction(affinities, protein_id, monomer)
        store.close()

    def get_protein_ids(self, molecules):
        proteins = filter(lambda x:x.mol_type == 'protein', molecules.values())
        return [mol.pdb_id for mol in proteins]

    def get_affinities(self, database, pdb_id, monomer, strength):
        '''read interaction values from database
        '''
        affinities = database._load_table('/'+pdb_id, monomer)
        return affinities

    def add_interaction(self, affinities, pdb_id, monomer_name):
        proteins = filter(lambda x:x.mol_type == 'protein', self.molecules.values())
        for protein in filter(lambda x:x.pdb_id == pdb_id, proteins):
            for chain,res_id,iCode,type_,value in affinities: 
                residue_id = (chain, res_id, iCode)
                particle = self._find_particle_by_pdb_id(residue_id, protein) 
                if not particle.type_.unique: 
                    protein_particle_type = self._make_particle_unique(particle)
     
                self.env.ff['affinity'][(particle.type_, self.env.atom_type[monomer_name])].epsilon = value 

    def add_auto_repulsion(self, lmpObj, repulsion_value):
        if not 'affinity' in self.env.ff:
            raise Exception('Affinity force field is not available, which is needed for repulsion.')
        # get Id List
        typeList = self._list_all_monomers()
        repulsionList = sorted([p for p in self.env.atom_type.values() 
                                    if p.name in typeList], key=lambda x:x.Id, reverse=True)
        for i, type1 in enumerate(repulsionList):
            for j, type2 in enumerate(repulsionList[i:]):
                self.env.ff['affinity'][(type1,type2)].epsilon = self.repulsion_value
