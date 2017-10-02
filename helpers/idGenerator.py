class IdGen(object):
    '''generates ids for every key you throw at it.
    If a key is not registered, it will generate a new
    counter that starts with one.
    '''
    def __init__(self):
        self._counters = {}
    
    def __getitem__(self, key):
        try:
            return self._counters[key].next()
        except KeyError:
            return self.__missing__(key)
    
    def __contains__(self, key):
        if key in self._counters.keys():
            return True
        else:
            return False

    def __missing__(self, key):
        self._counters[key] = self._id_generator()
        return self._counters[key].next()
    
    def keys(self):
        return self._counters.keys()

    def _id_generator(self):
        counter = 1
        while True:
            yield counter
            counter += 1
