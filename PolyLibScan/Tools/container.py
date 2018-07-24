
class Container(object):
    """Stores types very similar to a 
    dictionary with the added features of 
    automatic id generation for every item and
    hooks for updating dependent objects."""

    def __init__(self, type_name):
        super(Container, self).__init__()
        self.name = type_name
        self.new_type_id = self._idGenerator()
        self._types = {}
        self._defined_types = {}
        self._registered = []
        
    def __repr__(self):
        out = ['Type Category: %s' % self.name]
        for key, type_ in self._types.iteritems():
            out.append('%s: %s' % (key, repr(type_)))
        return '\n'.join(out)

    def __getitem__(self, key):
        try:
            return self._types[key]
        except KeyError:
            if key in self._defined_types.keys():
                # adding defined type to available types using 
                # __setitem__ method.
                self[key] = self._defined_types[key]
                return self._types[key]
            else:
                raise KeyError('Type %s not defined. Add type to config file.' % key)


    def __contains__(self, key):
        if key in self._types.keys() or key in self._defined_types.keys():
            return True
        else:
            return False

    def __len__(self):
        return len(self._types)

    def __iter__(self):
        return iter(self._types)

    def __delitem__(self, key):
        del self._types[key]

    def register(self, obj):
        self._registered.append(obj)

    def _update_hook(self, value):
        '''update registered objects via their update function
        when new ids are added to the atom type list.
        '''
        for obj in self._registered:
            obj.update(value)

    def __setitem__(self, key, value):
        type_id = self.new_type_id.next()
        if key in self._types.keys():
            unique_key = key + str(type_id)
        else:
            unique_key = key
        self._types[unique_key] = value
        self._types[unique_key].Id = type_id
        # update hooked object
        self._update_hook(value)

    def define_type(self, key, value):
        self._defined_types[key] = value

    def keys(self):
        return self._types.keys()

    def values(self):
        return self._types.values()

    def items(self):
        return self._types.items()

    def iteritems(self):
        return key,value in self._types.iteritems()

    def _idGenerator(self):
        '''generate ids of an infinite loop
        starting with 1.
        '''
        counter = 1
        while True:
            yield counter
            counter += 1
