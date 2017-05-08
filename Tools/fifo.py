import os
import subprocess as sp
import tempfile as temp

class FiFo(object):
    '''Generic Fifo class implementing default functionality.
    '''
    def __init__(self, fifo_path, script, out_file, step_size):
        self._fifo_path = ''
        self.fifo_path = fifo_path
        self.script = script
        self.out_file = out_file
        self.step_size = step_size
        self._process = None

    def __repr__(self):
        return 'Script: %s - Path: %s' % (os.path.basename(self.script), 
                                          os.path.basename(self.fifo_path))

    def create(self):
        os.mkfifo(self.fifo_path)
        if not os.path.exists(self.fifo_path):
            raise OSError("FIFO could not be created. Maybe you need to set local to 1")

    def terminate(self):
        os.remove(self.fifo_path)
        if os.path.exists(self.fifo_path):
            raise OSError('Could not delete fifo %s' % fifo['path'])

    def _script_exists(self):
        if not os.path.exists(self.script):
            raise Exception('Script %s not found!' % self.script)

    def _fifo_exists(self):
        '''Checks if the fifo file exists.
        '''
        if os.path.exists(self.fifo_path):
            return True
        else:
            return False

    def create_temp_folder(self):
        return temp.mkdtemp()

    def activate(self, index):
        out_file = self.output_path(index)
        self._start(out_file)

    def _start(self, out_path):
        '''the args argument can be used by scripts 
        to shovel additional parameters into it that
        are needed for calculation.
        '''
        self._process = sp.Popen('python "%s" "%s" "%s" "%s"' % (self.script, self.fifo_path, 
                                                                out_path,    self.args     ), 
                     stdin=None, stdout=None, stderr=None, shell=True, close_fds=True)

    def additional_arguments(self):
        raise NotImplementedError('still needs implementation')

    def fifo_path():
        doc = """The fifo_path property automates the creation and 
        removal of the fifo file."""
        def fget(self):
            return self._fifo_path
        def fset(self, value):
            if self._fifo_exists():
                self.terminate()
            self._fifo_path = value
            if not self._fifo_exists():
                self.create()
        def fdel(self):
            self.terminate()
            self._fifo_path = ''
        return locals()
    fifo_path = property(**fifo_path())