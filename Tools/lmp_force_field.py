
class Interaction(object):
    '''Container object that stored all interaction parameters in 
    a matrix and takes care of all the paperwork.
    '''
    def __init__(self, environment, pair_type):
        self._env = environment
        self.pair_type = pair_type
        self._matrix = {}
        #self._create_matrix()

    def _add(self, type_):
        '''This function is only one newly added type, 
        not for the creation of the interaction matrix.
        '''
        types_existing_before = filter(lambda x:x.Id!=type_.Id, self._env.atom_type.values())
        self._matrix[type_.Id] = {old_type.Id: Pair(self.pair_type, type_, old_type) 
                                  for old_type in types_existing_before}
        self._matrix[type_.Id][type_.Id] = Pair(self.pair_type, type_, type_)

    def __getitem__(self, types):
        '''return the interaction information of type_a 
        and type_b.

        Input:
            types: [type_a, type_b]
        
        Output:
            Pair Object
        '''
        try:
            return self._matrix[types[0].Id][types[1].Id]
        except KeyError:
            return self._matrix[types[1].Id][types[0].Id]

    def _create_matrix(self):
        '''creates a dictionary with the different pair_coeffs.
        '''
        # create list sorted by id numbers; largest first
        atomTypeList = sorted(self._env.atom_type.values(), key=lambda x:x.Id, reverse=True)
        for i, type1 in enumerate(atomTypeList):
            self._matrix[type1.Id] = {}
            for type2 in atomTypeList[i+1:]:
                self._matrix[type1.Id][type2.Id] = Pair(self.pair_type, type1, type2)

    def __contains__(self, types):
        '''Check if the interaction between both types are present.

        This is checked via the IDs.

        Input:
            types = [Type_a, Type_b]

        Output:
            bool
        '''
        if types[0].Id in self._matrix.keys() and types[1].Id in self._matrix.keys():
            return True
        else:
            return False

    def update(self, new_type):
        '''Update hook for adding additional types to the system. 
        '''
        self._add(new_type)


class Pair(object):

    def __init__(self, pair_type, atom_type1, atom_type2, epsilon=None, alpha=None):
        '''
        '''
        self.pair_type = pair_type
        self.atom_type1 = atom_type1
        self.atom_type2 = atom_type2

        self.epsilon = epsilon
        self.alpha = alpha
        
        self.set_pair_type_functions(self.pair_type.kind)
        self.epsilon_fct()
        self.alpha_fct()

    def epsilon_fct(self):
        NotImplementedError('This method should be overridden.')

    def alpha_fct(self):
        NotImplementedError('This method should be overridden.')

    def set_pair_type_functions(self, pair_type):
        if pair_type == 'soft':
            self.epsilon_fct = self.soft_epsilon
            self.alpha_fct = self.soft_alpha
        elif pair_type == 'lj/cut':
            self.epsilon_fct = self.lj_epsilon
            self.alpha_fct = self.lj_alpha
        elif pair_type == 'morse':
            self.epsilon_fct = self.morse_epsilon
            self.alpha_fct = self.morse_alpha

    ## MORSE parameter functions
    #
    def morse_epsilon(self):
        '''Set default value, if epsilon is not set at initialisation.
        '''
        if self.epsilon == None:
            self.epsilon = self.pair_type.parameters['coeffs'][0]

    def morse_alpha(self):
        if self.alpha == None:
            self.alpha = self.pair_type.parameters['coeffs'][1]

    ## SOFT parameter functions without surface
    #
    def soft_epsilon(self):
        if self.epsilon == None:
            if self.atom_type1.hydrophobicity > 0 and self.atom_type2.hydrophobicity > 0:
                self.epsilon = self.pair_type.parameters['coeffs'][0] * 0.5 * (self.atom_type1.hydrophobicity + self.atom_type2.hydrophobicity)
            else:
                self.epsilon = 0.0

    def soft_alpha(self):
        self.alpha = None

    ## SOFT parameter functions with surface energy
    #
    def soft_epsilon(self):
        if self.epsilon == None:
            if self.atom_type1.surface_energy > 0.0 or self.atom_type2.surface_energy > 0.0:
                surface_energy = max(self.atom_type1.surface_energy, self.atom_type2.surface_energy)
                self.epsilon = surface_energy * (self.atom_type1.hydrophobicity + self.atom_type2.hydrophobicity)/2.0
            elif self.atom_type1.hydrophobicity > 0 and self.atom_type2.hydrophobicity > 0:
                self.epsilon = self.pair_type.parameters['coeffs'][0] * (self.atom_type1.hydrophobicity + self.atom_type2.hydrophobicity)/2.0
            else:
                self.epsilon = 0.0

    def soft_alpha(self):
        self.alpha = None

    ## LJ parameter functions
    #
    def lj_epsilon(self):
        '''Set default value, if epsilon is not set at initialisation.
        '''
        if self.epsilon == None:
            self._epsilon = self.pair_type.parameters['coeffs'][0]

    def lj_alpha(self):
        self.alpha = None

    def cutoff():
        doc = "The cutoff property."
        def fget(self):
            if self.pair_type.repulsive_only == 1:
                return self.sigma
            else:
                return self.pair_type.cutoff
            return self._cutoff
        def fset(self, value):
            raise Exception('Cannot be set')
        def fdel(self):
            del self._cutoff
        return locals()
    cutoff = property(**cutoff())

    def epsilon():
        doc = "The epsilon property."
        def fget(self):
            if self.atom_type1.interacting and self.atom_type2.interacting:
                return self._epsilon
            else:
                return 0.0
        def fset(self, value):
            self._epsilon = value
        return locals()
    epsilon = property(**epsilon())

    def sigma():
        doc = "The sigma property."
        def fget(self):
            return self.atom_type1.radius + self.atom_type2.radius
        def fset(self, value):
            raise Exception('Cannot be set')
        return locals()
    sigma = property(**sigma())

    def __repr__(self):
        return '%s | %s | %s | a: %f, e: %f, s: %f' % (self.pair_type, self.atom_type1.name, 
            self.atom_type2.name, self.alpha, self.epsilon, self.sigma)
