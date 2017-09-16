import os
import subprocess as sp
import tempfile as temp

class BaseFiFo(object):
    '''Generic Fifo class implementing default functionality.
    '''
    def __init__(self, job, fifo_path, script, out_file, step_size):
        self._fifo_path = ''
        self.parent = job
        self.fifo_path = fifo_path
        self.script = script
        self.out_file = out_file
        self.step_size = step_size
        self._process = None

    @classmethod
    def from_dict(cls, job, fifo_dict):
        return cls(job, 
            fifo_dict['path'], 
            fifo_dict['script'], 
            fifo_dict['out'], 
            fifo_dict['stepsize'])

    def __repr__(self):
        return 'Script: %s - Path: %s' % (os.path.basename(self.script), 
                                          os.path.basename(self.fifo_path))

    def generate_path(self, file_name):
        '''fifo files should reside in the fifo folder.
        '''
        if 'fifo' in self.parent.config.lmp_path:
            fifo_folder = self.parent.config.lmp_path['fifo']
        else:
            fifo_folder = self.create_temp_folder()
        return os.path.join(fifo_folder, file_name)

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
                                                                out_path, self.additional_arguments()
                                                                ), 
                     stdin=None, stdout=None, stderr=None, shell=True, close_fds=True)

    def additional_arguments(self):
        raise NotImplementedError('still needs implementation')

    def lammps_string(self):
        raise NotImplementedError('still needs implementation')

    def output_path(self, index):
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



class FiFo(object):

    @staticmethod
    def create(fifo_type, job, fifo_path, script, out_file, step_size):
        if fifo_type == 'distance_fifo':
            return DistanceFifo(job, data['path'], data['script'], data['out'], data['stepsize'])
        elif fifo_type == 'traj_compression':
            return TrajCompressionFifo(job, data['path'], data['script'], data['out'], data['stepsize'])
        else:
            raise NotImplementedError('this fifo based postproduction is not yet implemented.')

    @staticmethod
    def create_from_dict(fifo_type, job, fifo_dict):
        if fifo_type == 'distance_fifo':
            return DistanceFifo.from_dict(job, fifo_dict)
        elif fifo_type == 'traj_compression':
            return TrajCompressionFifo.from_dict(job, fifo_dict)
        else:
            raise NotImplementedError('this fifo based postproduction is not yet implemented.')


class TrajCompressionFifo(BaseFiFo):
    '''This class feeds the output of LAMMPS to 
    the python script that calculates the distance between 
    the polymer and the active site.
    '''
    def __init__(self, job, fifo_name, script, file_name, steps_size=100):
        super(TrajCompressionFifo, self).__init__(fifo_path, script, file_name, steps_size)

    def additional_arguments(self):
        '''The monomer IDs are needed to be able to discern
        between active site and polymer.
        '''
        return ' '

    def output_path(self, index):
        file_name = '%s%05d.xyz.gz' % (self.out_file, index)
        return os.path.join(self.parent.config.lmp_path['output'], file_name)

    def lammps_string(self):
        return 'dump fifo_traj solid xyz %d "%s"' % (self.step_size, self.fifo_path)


class DistanceFifo(BaseFiFo):
    '''This class feeds the output of LAMMPS to 
    the python script that calculates the distance between 
    the polymer and the active site.
    '''
    def __init__(self, job, fifo_name, script, file_name, steps_size=100):
        super(DistanceFifo, self).__init__(fifo_path, script, file_name, steps_size)

    def additional_arguments(self):
        '''The monomer IDs are needed to be able to discern
        between active site and polymer.
        '''
        return '-'.join(map(str, self.parent.config.lmp_parameter['monomer_ids']))

    def lammps_string(self):
        return 'dump fifo_distance distance_group xyz %d "%s"' % (self.step_size, self.fifo_path)

    def output_path(self, index):
        file_name = '%s%05d' % (self.out_file, index)
        return os.path.join(self.parent.config.lmp_path['output'], file_name)
