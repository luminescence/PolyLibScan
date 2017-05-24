from lmp_particlesAndInteractions import *
from lmp_helpers import *
from lmpObj import LmpObj

class LmpCreator(object):

    def __init__(self, environment):
        self.env = environment
        
        self.mol_type = 'unknown'
        
    def __repr__(self):
        description = 'Molecule Factory creates %s' % (self.mol_type)
        return description

    def __str__(self):
        string  = ['Molecule Factory creates %s:' % (self.mol_type)]
        string += ['Atom style: %s' % self.env.globals['atom_style']]
        string += ['Bond style: %s' % self.env.globals['bond_style']]
        string += ['Angle style: %s' % self.env.globals['angle_style']]
        return '\n'.join(string)

    def create(self):
        '''Creates the molecule in the default order that is identical for
        proteins and polymers. The internal methods
        - _create_config
        - _create_particles_and_interactions
        - _create_custom 
        are different for protein and polymers.
        The molecule in embedded into the environment of the creator and 
        also returned at the end.
        '''
        molecule = LmpObj(self.env)
        self._create_config(molecule)
        self._create_particles_and_interactions(molecule)
        self._create_custom(molecule)
        molecule.calc_box()
        # register molecule in environment
        self.env.add_molecule(molecule)
        return molecule
        
    def _create_config(self, molecule):
        raise NotImplementedError('Each subclass has to implement this on its own.')

    def _create_particles_and_interactions(self, molecule):
        raise NotImplementedError('Each subclass has to implement this on its own.')

    def _create_custom(self, molecule):
        raise NotImplementedError('Each subclass has to implement this on its own.')
