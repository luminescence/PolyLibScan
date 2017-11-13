from lammps import lammps


class lmp_controler(object):
    """initiate and control a lammps instance"""
    @staticmethod
    def set_variables(vars_, lmp_instance):
        for name, val in vars_.items():
            # lists are converted to strings
            if isinstance(val, list):
                val = ' '.join(map(str, val))
            lmp_instance.command('variable %s string "%s"'% (name,val))

    @staticmethod
    def set_paths(paths, lmp_instance):
        for name, path in paths.items():
            lmp_instance.command('variable %s string "%s"'% (name, path))

    @staticmethod
    def set_fifos(fifos, lmp_instance):
        for name,fifo in fifos.items():
            lmp_instance.command(fifo.lammps_string())

    @staticmethod
    def _start_fifo_capture(fifos, index):
        for name, fifo in fifos.items():
            fifo.activate(index)

    def lmps_run(self, Id, parameters, paths, fifos={}):
        self._start_fifo_capture(fifos, Id)
        lammps_sim = lammps()
        # submitting parameters
        self.set_variables(parameters, lammps_sim)
        # submitting paths
        self.set_paths(paths, lammps_sim)
        # submitting run Id
        lammps_sim.command('variable num string %05d'% Id)
        # starting script
        lammps_sim.file(paths['script'])
        # specify fifo dumps
        self.set_fifos(fifos, lammps_sim)
        lammps_sim.command('run ${time_steps}')
        # write snapshot of end-comformation
        lammps_sim.command('write_dump solid xyz ${end_xyz}.xyz')
        lammps_sim.close()