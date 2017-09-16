import tables as tb
import pathlib2 as pl

class JobDataBase(object):
    def __init__(self, path, mode='a'):
        super(JobDataBase, self).__init__()
        self.db = Database(path, mode='a')

    def close(self):
        self.db.close()

    def open(self):
        self.db.open()

    def create_groups(self):
        '''Creates the default groups in the database.
        '''
        self.db._create_group('energies', self.db._handle.root)
        self.db._create_group('distances', self.db._handle.root)
        self.db._create_group('start_trajectories', self.db._handle.root)
        self.db._create_group('end_trajectories', self.db._handle.root)
        self.db._create_group('meta', self.db._handle.root)
        self.db._create_group('trajectory', '/meta')
        self.db._create_group('trajectories', self.db._handle.root)

    @property
    def active_site(self):
        return self.db._load_table('/meta', 'sequence')
    
    @active_site.setter
    def active_site(self, value):
        self.db._save_table(value, '/meta', 'active_site')

    @property
    def sequence(self):
        return self.db._load_table('/meta', 'sequence')
    
    @sequence.setter
    def sequence(self, value):
        self.db._save_table(value, '/meta', 'sequence')
       
    @property
    def weights(self):
        return self.db._load_table('/meta', 'weights')
    
    @weights.setter
    def weights(self, value):
        self.db._save_table(value, '/meta', 'weights')

    @property
    def misc(self):
        return self.db._load_table('/meta', 'misc')
    
    @misc.setter
    def misc(self, value):
        self.db._save_table(value, '/meta', 'misc')

    @property
    def parameter(self):
        return self.db._load_table('/meta', 'parameter')
    
    @parameter.setter
    def parameter(self, value):
        self.db._save_table(value, '/meta', 'parameter')

    @property
    def versions(self):
        return self.db._load_table('/meta', 'versions')
    
    @versions.setter
    def versions(self, value):
        self.db._save_table(value, '/meta', 'versions')

    @property
    def histogramm(self):
        return self.db._load_table('/', 'histogramm')
    
    @histogramm.setter
    def histogramm(self, value):
        self.db._save_table(value, '/', 'histogramm')

    @property
    def traj_type_order(self):
        return self.db._load_table('/meta/trajectory/', 'type_order')

    @traj_type_order.setter
    def traj_type_order(self, value):
        self.db._save_table(value, '/meta/trajectory/', 'type_order')

    @property
    def traj_info(self):
        return self.db._load_table('/meta/trajectory/', 'info')

    @traj_info.setter
    def traj_info(self, value):
        self.db._save_table(value, '/meta/trajectory/', 'info')

    @property
    def end_states(self):
        return self.db._load_table('/', 'end_state')

    @end_states.setter
    def end_states(self, value):
        self.db._save_table(value, '/', 'end_state')

    def distance_ts_load(self, ID):
        return self.db._load_table('/distances', 'd%d' % ID)

    def distance_ts_save(self, distance_ts, ID):
        self.db._save_table(distance_ts, '/distances', 'd%d' % ID)

    def energie_ts_load(self, ID):
        return self.db._load_ctable('/energies', 'd%d' % ID)

    def energie_ts_save(self, energie_ts, ID):
        self.db._save_array(energie_ts, '/energies', 'd%d' % ID)

    def start_trajectories_load(self, ID):
        return self.db._load_table('/start_trajectories', 't%d' % ID)

    def start_trajectories_save(self, trajectorie, ID):
        self.db._save_table(trajectorie, '/start_trajectories', 't%d' % ID)

    def end_trajectories_load(self, ID):
        return self.db._load_table('/end_trajectories', 't%d' % ID)

    def end_trajectories_save(self, trajectorie, ID):
        self.db._save_table(trajectorie, '/end_trajectories', 't%d' % ID)

    def trajectory_load(self):
        return self.db._load_ctable('/trajectories', 'traj_%d' % Id, as_iterator=True)

    def trajectory_save(self, data, ID):
        self.db._save_array(data, '/trajectories', 'traj_%d' % ID, compress=True)


class Database(object):
    """docstring for Database"""
    def __init__(self, path, mode='a'):
        super(Database, self).__init__()
        self.path = path
        self.mode = mode
        self._handle = None
        self.open()
        
    def __repr__(self):
        return repr(self._handle)

    def open(self, path=None, mode=None):
        # use local parameters if not set
        if not path:
            path = self.path
        if not mode:
            mode = self.mode

        # Make sure path is not pathlib
        if isinstance(path, pl.PosixPath):
            path = path.as_posix()
        f = tb.open_file(path, mode)
        self._handle = f

    def is_open(self):
        if self._handle.isopen == 0:
            return False
        else:
            return True

    def close(self):
        self._handle.close()

    def _save_table(self, data, group, table_name):
        if not self.is_open():
            self.open()
        table = self._handle.create_table(group, table_name , data.dtype,
                                          expectedrows=data.shape[0])
        for d in data:
            for name, item in zip(data.dtype.names, d):
                table.row[name] = item
            table.row.append()
        table.flush()

    def _save_array(self, data, group, table_name, compress=False):
        if not self.is_open():
            self.open()
        if compress:
            compressor = tb.Filters(complevel=5, complib='blosc')
            self._handle.create_carray(group, table_name , obj=data,
                                       filters=compressor)
        else:
            self._handle.create_carray(group, table_name , obj=data)
        
    def _create_group(self, group_name, location='/'):
        self._handle.create_group(location, group_name)
        
    def _load_table(self, group, table_name, as_iterator=False):
        if not self.is_open():
            self.open()
        if as_iterator:
            return self._handle.get_node(group, table_name).iterrows()
        else:
            return self._handle.get_node(group, table_name)[:]

    def _load_ctable(self, group, table_name, col=None, as_iterator=False):
        """read carray from the database and return the data as numpy array 
        or iterator. If the data is returned as an iterator, col is ignored.

        Arguments:
            group -- location of table [string]
            table_name -- name of the table [string]
            col -- column of the array [integer]
            as_iterator -- switch between returning array or iterator [bool]

        Return:
            numpy array or iterator
        """
        if not self.is_open():
            self.open()
        if as_iterator:
            return self._handle.get_node(group, table_name).iterrows()
        if col is not None:
            return self._handle.get_node(group, table_name)[:,col]
        else:
            return self._handle.get_node(group, table_name)[:]
