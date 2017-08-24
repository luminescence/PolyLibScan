
class Interaction(object):

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
        if epsilon == None:
            if self.pair_type.kind == 'soft':
                if self.atom_type1.hydrophobicity > 0 or self.atom_type2.hydrophobicity > 0:
                    self._epsilon = self.pair_type.parameters['coeffs'][0] * 0.5 * (self.atom_type1.hydrophobicity + self.atom_type2.hydrophobicity)
                else:
                    self._epsilon = 0.0
            else:
                self._epsilon = self.pair_type.parameters['coeffs'][0]
        else:
            self._epsilon = epsilon
        if self.pair_type.kind == 'morse' and alpha == None:
            self.alpha = self.pair_type.parameters['coeffs'][1]
        else:
            self.alpha = alpha

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
