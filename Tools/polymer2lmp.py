import itertools as it
import numpy as np
from lmp_particlesAndInteractions import *
import monomers as Mono
from lmpObj import LmpObj
from lmp_creator import LmpCreator
import lmp_helpers as helpers
import sets

class PolymerCreator(LmpCreator):
    
    def __init__(self, environment, pattern, weights=None, 
                 auto_repulsion=False, mode='random', length=None):
        super(PolymerCreator, self).__init__(environment)
        
        self.mol_type = 'polymer'
        self.pattern = pattern
        self.length = length


        # Mode initialized
        self._mode = None
        self.mode = mode

        #Weights initialized
        self.weights = weights

        self.auto_repulsion = auto_repulsion
        self.repulsion_value = -50


    def __str__(self):
        string  = ['Molecule Factory creates: %s' % (self.mol_type.upper())]
        string += ['Atom style:  %s' % self.env.globals['atom_style']]
        string += ['Bond style:  %s' % self.env.globals['bond_style']]
        string += ['Angle style: %s' % self.env.globals['angle_style']]
        string += ['Polymer length: %d' % self.length]
        string += ['Polymer Creation Mode: %s' % self.mode.upper()]
        string += ['Monomer repulsion: %s' % self.auto_repulsion]
        return '\n'.join(string)

    ### MANDATORY CREATE FUNCTIONS
    def _create_config(self, molecule):
        molecule.mol_type = self.mol_type

    # _create_types() is implemented in parent class
    
    def _create_particles_and_interactions(self, molecule):
        molecule.data['monomers'] = self.create_polymer(molecule)
        molecule.data['particles'] = self.create_particles(molecule)
        molecule.data['bonds'] = self.create_polymer_bonds(molecule)
        molecule.data['angles'] = self.create_polymer_angles(molecule)
        molecule.data['dihedrals'] = self.create_polymer_dihedrals(molecule)

    def _create_custom(self, molecule):
        if self.auto_repulsion:
            self.add_auto_repulsion(molecule)


    ### HELPER FUNCTIONS
    def _generator_random(self):
        '''creates a generator that yields a random series of 
        the specified monomers. The length of the series is determined
        by self.length or by the number of monomers.
        '''
        if type(self.weights) == list:
            weights = helpers.normalize(self._weights)
        elif self.weights == None:
            weights = None
        else:
            raise AttributeError('weights type not "list" or None.')
            
        if self.length != None:
            for monomer in np.random.choice(self.pattern, size=self.length, p=weights):
                yield monomer
        else:
            for monomer in np.random.choice(self.pattern, size=len(self.pattern), p=weights):
                yield monomer

    def _generator_cyclic(self):
        '''creates a generator that yields a repeating series of 
        the specified monomers. The length of the series is determined 
        by self.length or, if length is not specified, the pattern list is 
        returned once.
        '''
        if self.length != None:
            polymer_pattern = it.cycle(self.pattern)
            for i in range(self.length):
                yield polymer_pattern.next()
        else:
            for monomer in self.pattern:
                yield monomer

    def polymer_list(self):
        '''if length is set, the pattern will be continued until the length 
        is reached
        '''
        raise NotImplementedError('Something clearly went wrong. Set mode to "cycle" or "random"')

    def create_particles(self, molecule):
        particles = set([])
        for mono in molecule.data['monomers']:
            for particle, value in mono.particles.items():
                 particles.add(value)
        return sorted(list(particles), key=lambda x:x.Id)
    
    def create_polymer(self, molecule):
        polymer = []
        for element in self.polymer_list():
            polymer.append(self.add_monomer(element, polymer, molecule))
        return polymer
    
    def add_monomer(self, element, polymer, lmpObj):
        if len(polymer) == 0:
            coords = np.array([0.0, 0.0, 0.0])
        else:
            coords = polymer[-1].position + np.array([4.0, 0.0, 0.0])
        return Mono.Monomer(coords, element, self.env.monomer_type[element] ,self.env, lmpObj)

    @staticmethod
    def unique_ordered_DOF(monomers, DOF):
        all_DOF = []
        for x in monomers:
            all_DOF += getattr(x, DOF)
        ordered_unique_bonds = sorted(list(set(tuple(all_DOF))), key=lambda x: x.Id)
        return ordered_unique_bonds

    def create_polymer_bonds(self, molecule):
        '''links all the monomers to 
        one string together.
        '''
        monomers = molecule.data['monomers']
        if len(monomers) < 2:
            return []
        for current_mono,next_mono in zip(monomers[0:-1], monomers[1:]):
            current_mono.bind_with(next_mono)

        ordered_unique_bonds = self.unique_ordered_DOF(monomers, 'bonds')
        return ordered_unique_bonds

    def create_polymer_angles(self, molecule):
        monomers = molecule.data['monomers']
        if len(monomers) >= 2:
            monomers[0].angle_with(monomers[1])
            monomers[-1].angle_with(monomers[-2])
        if len(monomers) >= 3:
            for previous_element, current_element, next_element in zip(
                                                    monomers[:-2], 
                                                    monomers[1:-1], 
                                                    monomers[2:]):
                current_element.angle_with(next_element)
                current_element.bb_angle_with(previous_element, next_element)

        ordered_unique_angles = self.unique_ordered_DOF(monomers, 'angles')
        return ordered_unique_angles

    def create_polymer_dihedrals(self, molecule):
        monomers = molecule.data['monomers']
        if len(monomers) >= 2:
            molecule.data['monomers'][-1].dihedral_with(molecule.data['monomers'][-2])
        if len(monomers) >= 3:
            molecule.data['monomers'][-2].dihedral_with(molecule.data['monomers'][-3])
        if len(monomers) >= 4:
            for previous_element, current_element, next_element, after_next_element in zip(
                molecule.data['monomers'][:-3], 
                molecule.data['monomers'][1:-2],
                molecule.data['monomers'][2:-1], 
                molecule.data['monomers'][3:]):
                current_element.dihedral_with(previous_element)
                current_element.bb_dihedral_with(previous_element, next_element, after_next_element)

        ordered_unique_dihedrals = self.unique_ordered_DOF(monomers, 'dihedrals')
        return ordered_unique_dihedrals

    def add_auto_repulsion(self, lmpObj, repulsion_value=None):
        '''Make all monomers repulse each other.

        Input:
            lmpObj: Polymer object
            repulsion_value: [float]

        Output:
            None
        '''
        if repulsion_value:
            epsilon = repulsion_value
        else:
            epsilon = self.repulsion_value
        # get Id List
        monos = list(set([p.type_ for p in lmpObj.data['particles']]))
        repulsionList = sorted(monos, key=lambda x:x.Id, reverse=True)
        for i, type1 in enumerate(repulsionList):
            for j, type2 in enumerate(repulsionList[i:]):
                self.env.ff['affinity'][(type1,type2)].epsilon = epsilon

    # PROPERTIES
    def mode():
        doc = "The mode property."
        def fget(self):
            return self._mode
        def fset(self, value):
            if value not in ['random', 'cycle']:
                raise ValueError('mode must be set to "random" or "cycle".')
            self._mode = value
            if value == 'random':
                self.polymer_list = self._generator_random
            elif value == 'cycle':
                self.polymer_list = self._generator_cyclic
        return locals()
    mode = property(**mode())

    def weights():
        doc = "The weights property."
        def fget(self):
            return self._weights
        def fset(self, value):
            if value == None:
                self._weights = None
            elif len(value) == len(self.pattern):
                self._weights = value
            else:
                raise ValueError('weights must be a list and must have the same length as pattern.')
        return locals()
    weights = property(**weights())
