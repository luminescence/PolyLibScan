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
        molecule.data['particles'] =  self.create_polymer(molecule)
        molecule.data['bonds'] = self.create_polymer_bonds(molecule)
        molecule.data['angles'] = self.create_polymer_angles(molecule)
        molecule.data['dihedrals'] = [] # not implemented

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
            polymer.append(self.add_particle(element, polymer, molecule))
        return polymer
    
    def add_particle(self, element, polymer, lmpObj):
        if len(polymer) == 0:
            coords = np.array([0.0, 0.0, 0.0])
        else:
            coords = polymer[-1].position + np.array([4.0, 0.0, 0.0])
        return Mono.Monomer(coords, element, self.env.monomer_type[element] ,self.env, lmpObj)

    def create_polymer_bonds(self, molecule):
        '''links all the particles to 
        one string together.
        '''
        bonds = []
        for i, monomer in enumerate(molecule.data['particles'][:-1], 1):
            bonds.append(Bond(self.env.new_id['bond'], molecule.env.bond_type['polymer'], [monomer, molecule.data['particles'][i]]))
        return bonds
    
    def create_polymer_angles(self, molecule):
        angles = []
        for i, monomer in enumerate(molecule.data['particles'][:-2], 1):
            angles += [Angle(self.env.new_id['angle'], molecule.env.angle_type['polymer'], 
                            [monomer, molecule.data['particles'][i], molecule.data['particles'][i+1]])]
        return angles


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
