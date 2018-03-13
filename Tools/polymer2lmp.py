import itertools as it
import numpy as np
from lmp_particlesAndInteractions import *
import monomers as Mono
from lmpObj import LmpObj
from lmp_creator import LmpCreator
import lmp_helpers as helpers
import sets

class PolymerCreator(LmpCreator):
    
    def __init__(self, environment, pattern, poly_type=None, weights=None, 
                 auto_repulsion=False, mode='random', length=None):
        super(PolymerCreator, self).__init__(environment)
        
        self.mol_type = 'polymer'
        self.pattern = pattern
        self.length = length
        if poly_type is None:
            # take first
            poly_type = self.env.polymer_type._defined_types.keys()[0]
        self.type = self.env.polymer_type[poly_type]

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
        molecule.type = self.type

    # _create_types() is implemented in parent class
    
    def _create_particles_and_interactions(self, molecule):
        molecule.data['monomers'] = self.create_polymer(molecule)
        self.create_polymer_bonds(molecule)
        self.create_polymer_angles(molecule)
        self.create_polymer_dihedrals(molecule)

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
    
    def create_polymer(self, molecule):
        polymer = []
        for element in self.polymer_list():
            polymer.append(self.add_monomer(element, polymer, molecule))
        return polymer
    
    def add_monomer(self, element, polymer, lmpObj):
        if len(polymer) == 0:
            coords = np.array([0.0, 0.0, 0.0])
        else:
            coords = polymer[-1].position + np.array(self.type.offset)
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

    def create_polymer_angles(self, molecule):
        '''Create angles between monomers in the polymer.

        There are three possible types of angles between a monomer and its neighbors.
        We enumerate them with 1,2,3 where 1 is the monomer to the left,
                                           2 is the central monomer and 
                                           3 is the monomer to the right:
            - between a monomer and the monomer to the left (1-2)
            - between a monomer and the monomer to the right (2-3)
            - between a monomer and the monomers to the left and right (1-2-3)

        Angles involving two monomers are subsumed into one iteration and 
        angles involving three monomers have its own loop
        '''
        monomers = molecule.data['monomers']
        # 2-Partner loop
        angle_info = filter(lambda x: len(x)==3, molecule.type.angles)
        for current_element,next_element in zip(
                                                monomers[0:-1], 
                                                monomers[1:]):
            current_element.angles_with_one(next_element, angle_info)
        # 1-2-3 Angles
        angle_info = filter(lambda x: len(x)==4, molecule.type.angles)
        for previous_element, current_element, next_element in zip(
                                                monomers[0:-2], 
                                                monomers[1:-1], 
                                                monomers[2:]):
            current_element.angles_with_two(previous_element, next_element, angle_info)

    def create_polymer_dihedrals(self, molecule):
        monomers = molecule.data['monomers']
        # 1-2 Dihedrals
        parameters = filter(lambda x: len(x)==3, molecule.type.dihedrals)
        print 'pars:', parameters
        for current_element,next_element in zip(
                                                monomers[0:-1], 
                                                monomers[1:]):
            current_element.dihedral_with_two(next_element, parameters)
        # 1-2-3 Dihedrals
        parameters = filter(lambda x: len(x)==4, molecule.type.dihedrals)
        for current_element,next_element,next_next in zip(
                                                monomers[0:-2], 
                                                monomers[1:-1], 
                                                monomers[2:]):
            current_element.dihedral_with_three(next_element, next_next, parameters)
        # 1-2-3-4 Dihedrals
        parameters = filter(lambda x: len(x)==5, molecule.type.dihedrals)
        for current_element, next_element,next_next,next_next_next in zip(
                                                monomers[0:-3], 
                                                monomers[1:-2], 
                                                monomers[2:-1], 
                                                monomers[3:]):
            current_element.dihedral_with_four(next_element, next_next, next_next_next, parameters)

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
