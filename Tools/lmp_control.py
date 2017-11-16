from collections import OrderedDict
from lammps import PyLammps


class lmp_controler(object):
    """Initialize and control a lammps instance"""

    def __init__(self, Id, parameters, paths, fifos={}, previous_instance=''):
        if isinstance(previous_instance, PyLammps):
            self.lmp_instance = previous_instance
        elif previous_instance != '':
            raise ValueError('The previous_instance passed to lmp_controler is not a PyLammps instance!')
        else:
            self.lmp_instance = PyLammps()

            self.Id = Id
            self.parameters = parameters
            self.paths = paths
            self.fifos = fifos

    def set_dictionary_as_lammps_variables(self, dictionary, var_style):
        """dictionary keys will be variable names, dictionary values will be variable values"""
        for name, val in dictionary.items():
            # lists are converted to strings
            if isinstance(val, list):
                val = self.convert_python_list_to_lammps_list(val)
            self.lmp_variable(name, var_style, val)

    def lmp_variable(self, name, style, val):
        """
        >>> test_controler = lmp_controler()
        LAMMPS output is captured by PyLammps wrapper
        >>> test_controler.lmp_variable('a', 'string', '2')
        >>> test_controler.lmp_instance.variables['a'].value == 2.0
        True"""
        self.lmp_instance.variable('%s %s %s' % (name, style, val))

    @staticmethod
    def convert_python_list_to_lammps_list(val, with_quotes=True):
        """convert from [a,b,c] to "a b c" which is required in lammps"""
        val = ' '.join(map(str, val))
        if with_quotes:
            val = '"' + val + '"'
        return val

    @staticmethod
    def variable_sigil_for_lammps(varname):
        """correctly wrap ${...} around variable name"""
        return "${%s}" % varname

    def set_fifos(self, fifos):
        for name,fifo in fifos.items():
            self.lmp_instance.command(fifo.lammps_string())

    @staticmethod
    def _start_fifo_capture(fifos, index):
        for name, fifo in fifos.items():
            fifo.activate(index)

    def default_settings(self):

        #initialization
        initialization_cmds = ['units		real',
                               'boundary	p p p',
                               'atom_style	hybrid angle charge',
                               'bond_style  harmonic',
                               'angle_style harmonic']
        self.lmp_execute_list_of_commands(initialization_cmds)

        #reading coordinate file
        self.lmp_instance.command('read_data	%s/%05d' % (self.paths['input'], self.Id))

    def lmp_minimization(self):
        # minimization
        self.lmp_instance.command('minimize 0.0 1.0e-8 1000 100000')

    def configure_output(self):
        # output & computes
        self.lmp_instance.command('compute 2 contacts group/group polymer pair yes')
        self.lmp_instance.command('compute SolidTemp solid temp')
        output_vars = OrderedDict([('group0', 'c_2'),
                                   ('energy_', 'etotal'),
                                   ('potential_e', 'pe'),
                                   ('kinetic_e', 'ke'),
                                   ('sim_temp', 'c_SolidTemp'),
                                   ('time_step', 'step')])
        self.set_dictionary_as_lammps_variables(output_vars, 'equal')
        self.lmp_instance.command('fix 		5 all print 100 %s file %s/Energy%s screen no' % (
        self.convert_python_list_to_lammps_list([self.variable_sigil_for_lammps(x) for x in output_vars.keys()]), self.paths['output'], self.Id))

    def lmp_production_MD(self, temperature_production):
        # production
        self.lmp_instance.command('reset_timestep 0')
        self.langevin_dynamics(T_0=temperature_production, T_end=temperature_production)
        self.lmp_instance.run(self.parameters['time_steps'])
        # write snapshot of end-comformation
        output_file_final_conformation = self.paths['output'] + '/trajectoryE' + '%05d' % self.Id
        self.lmp_instance.command('write_dump solid xyz %s.xyz' % output_file_final_conformation)

    def lmp_equilibration_MD(self, temperature_production, temperature_start):
        output_file_start_conformation = self.paths['output'] + '/trajectoryS' + '%05d' % self.Id

        equilibration_cmds = ['velocity 	solid create %s 1321 dist gaussian' % temperature_start,
                              'write_dump	solid xyz %s.xyz' % output_file_start_conformation]
        self.lmp_execute_list_of_commands(equilibration_cmds)
        self.langevin_dynamics(T_0=temperature_start, T_end=temperature_production)
        equilibration_timesteps = 170
        self.lmp_instance.run(equilibration_timesteps)

    def group_declaration(self):
        # grouping
        groups = ['group protein type %s' % self.parameters['bb_id'],
                  'group ghost_protein type %s' % self.parameters['ghost_id'],
                  'group polymer type %s' % self.convert_python_list_to_lammps_list(self.parameters['monomer_ids'], with_quotes=False),
                  'group contacts subtract all protein ghost_protein polymer',
                  'group interactions union polymer contacts',
                  'group solid subtract all ghost_protein',
                  'group activesite id %s' % self.convert_python_list_to_lammps_list(self.parameters['active_site_ids'], with_quotes=False),
                  'group distance_group union polymer activesite']
        self.lmp_execute_list_of_commands(groups)

    def langevin_dynamics(self, T_0, T_end):
        cmds = ['fix 		3 solid  langevin %s %s 100.0 699483 gjf yes' % (T_0, T_end),
                'fix 		4 solid nve']
        self.lmp_execute_list_of_commands(cmds)

    def model_specifications(self):
        # specify the model
        pair_coeffs = ['timestep 	%s' % self.parameters['timestep'],
                       'pair_style hybrid/overlay lj96/cut 10.5 coul/debye %s 25.0' % self.parameters['debye_kappa'],
                       'dielectric %s' % self.parameters['dielectric_par'],
                       'pair_coeff * * lj96/cut 0.14   4.00   8.50',
                       'pair_coeff * * coul/debye']
        self.lmp_execute_list_of_commands(pair_coeffs)

    def lmp_execute_list_of_commands(self, list_of_commands):
        for line in list_of_commands:
            self.lmp_instance.command(line)

    def lmps_run(self, run_script=False):
        self._start_fifo_capture(self.fifos, self.Id)
        # submitting run Id
        self.lmp_variable('num', 'string','%05d' % self.Id)
        # use default settings
        self.default_settings()
        # starting script if desired
        if run_script:
            self.lmp_instance.file(self.paths['script'])
        self.model_specifications()
        self.group_declaration()
        #modify neighbor list
        self.lmp_instance.command('neigh_modify exclude group ghost_protein all')
        self.lmp_minimization()
        #equilibrate
        temperature_start = 100.0
        temperature_production = 300.0
        self.lmp_equilibration_MD(temperature_production, temperature_start)
        self.configure_output()
        self.lmp_production_MD(temperature_production)
        # specify fifo dumps
        self.set_fifos(self.fifos)
        self.lmp_instance.close()



if __name__ == "__main__":
    import doctest

    doctest.testmod()