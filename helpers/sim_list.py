import os
from . import time as tm

class SimState(object):
    """docstring for SimState"""
    def __init__(self, line):
        super(SimState, self).__init__()
        self.line = line
        self.process_line(self.line)

    @classmethod
    def create(cls, sim_id, path):
        id_str = '%05d'% sim_id
        hostname = os.uname()[1] # localhost string
        current_time = tm.time_string()
        line = ';'.join([id_str,hostname,path,current_time])
        return cls(line)
        
    def process_line(self, line):
        data = line.split(';')
        
        self.id, self.node, self.path, self.time = data

    def __repr__(self):
        return self.line
        
    @property
    def id(self):
        return self._id

    @id.setter
    def id(self, value):
        self._id = int(value)

    @property
    def time(self):
        return self._time

    @time.setter
    def time(self, value):
        self._time = tm.conversion(value)

class SimList(object):
    '''This class manages the simulation runs.
    '''
    
    def __init__(self, path, data_folder, sampling_rate):
        self.list_path = path
        self.data_folder = data_folder
        self.sampling_rate = sampling_rate
        self.sims = []
        self._read_list()
        
    def _read_list(self):
        if not os.path.exists(self.list_path):
            open(self.list_path, 'a').close()
        else:
            with open(self.list_path, 'r') as f:
                for line in f.read().splitlines():
                    self.sims.append(SimState(line))
        
    def mark_complete(self, index):
        new_completed_sim = SimState.create(index, self.data_folder)
        self.sims.append(new_completed_sim)
        with open(self.list_path, 'a') as f:
            info = '%s\n' % new_completed_sim
            f.write(info)

    def __iter__(self):
        full_id_list = range(0, self.sampling_rate)
        complete_sim_ids = [sim.id for sim in self.sims]
        remaining_ids = (i for i in full_id_list if i not in complete_sim_ids)
        return remaining_ids
