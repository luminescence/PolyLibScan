import tables as tb
import pathlib2 as pl

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
        table = self._handle.create_table(group, table_name , data.dtype,
                                          expectedrows=data.shape[0])
        for d in data:
            for name, item in zip(data.dtype.names, d):
                table.row[name] = item
            table.row.append()
        table.flush()

    def _save_array(self, data, group, table_name, compress=False):
        if compress:
            compressor = tb.Filters(complevel=5, complib='blosc')
            self._handle.create_carray(group, table_name , obj=data,
                                       filters=compressor)
        else:
            self._handle.create_carray(group, table_name , obj=data)
        
    def _create_group(self, group_name, location='/'):
        self._handle.create_group(location, group_name)
        
    def _load_table(self, group, table_name, as_iterator=False):
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
        if as_iterator:
            return self._handle.get_node(group, table_name).iterrows()
        if col is not None:
            return self._handle.get_node(group, table_name)[:,col]
        else:
            return self._handle.get_node(group, table_name)[:]
