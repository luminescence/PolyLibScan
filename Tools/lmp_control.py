from lammps import PyLammps


class lmp_controler(object):
    """Initialize and control a lammps instance"""

    def __init__(self, previous_instance=''):
        if isinstance(previous_instance, PyLammps):
            self.lmp_instance = previous_instance
        elif previous_instance != '':
            raise ValueError('The previous instance passed to lmp_controler is not a PyLammps instance!')
        else:
            self.lmp_instance = PyLammps()

    def set_dictionary_as_lammps_variables(self, dictionary, var_style):
        """dictionary keys will be variable names, dictionary values will be variable values"""
        for name, val in dictionary.items():
            # lists are converted to strings
            if isinstance(val, list):
                val = self.convert_python_list_to_lammps_list(val)
            self.lmp_variable(name, var_style, val)

    def lmp_variable(self, name, style, val):
        self.lmp_instance.variable('%s %s %s' % (name, style, val))

    @staticmethod
    def convert_python_list_to_lammps_list(val):
        """convert from [a,b,c] to "a b c" which is required in lammps"""
        val = ' '.join(map(str, val))
        val = '"' + val + '"'
        return val

    def set_fifos(self, fifos):
        for name,fifo in fifos.items():
            self.lmp_instance.command(fifo.lammps_string())

    @staticmethod
    def _start_fifo_capture(fifos, index):
        for name, fifo in fifos.items():
            fifo.activate(index)

    def default_settings(self):

        #variables
        file_paths = {'dat_file': '${input}/${num}',
                      'start_xyz': '${output}/trajectoryS${num}',
                      'end_xyz': '${output}/trajectoryE${num}',
                      'group_output': '${output}/Energy${num}'}
        self.set_dictionary_as_lammps_variables(file_paths, 'string')

        #initialization
        initialization_cmds = ['units		real',
                               'timestep 	${timestep}',
                               'boundary	p p p',
                               'atom_style	hybrid angle charge',
                               'bond_style  harmonic',
                               'angle_style harmonic',
                               'pair_style hybrid/overlay lj96/cut 10.5 coul/debye ${debye_kappa} 25.0',
                               'dielectric ${dielectric_par}']
        self.lmp_execute_list_of_commands(initialization_cmds)

        #reading coordinate file
        self.lmp_instance.command('read_data	${dat_file}')
        pair_coeffs = ['pair_coeff * * lj96/cut 0.14   4.00   8.50',
                       'pair_coeff * * coul/debye']
        self.lmp_execute_list_of_commands(pair_coeffs)

        #grouping
        groups = ['group protein type ${bb_id}',
                  'group ghost_protein type ${ghost_id}',
                  'group polymer type ${monomer_ids}',
                  'group contacts subtract all protein ghost_protein polymer',
                  'group interactions union polymer contacts',
                  'group solid subtract all ghost_protein',
                  'group solid type 1 2',
                  'group activesite id ${active_site_ids}',
                  'group distance_group union polymer activesite']
        self.lmp_execute_list_of_commands(groups)

        #modify neighbor list
        self.lmp_instance.command('neigh_modify exclude group ghost_protein all')

        #computes
        self.lmp_instance.command('compute 2 contacts group/group polymer pair yes')
        self.lmp_instance.command('compute SolidTemp solid temp')

        # minimization
        self.lmp_instance.command('minimize 0.0 1.0e-8 100 1000')

        #equilibrate
        equilibration_cmds = ['velocity 	solid create 100.0 1321 dist gaussian',
                              'fix 		3 solid  langevin 1.0 300.0 100.0 699483 gjf yes',
                              'fix 		4 solid nve',
                              'write_dump	solid xyz ${start_xyz}.xyz']
        self.lmp_execute_list_of_commands(equilibration_cmds)
        equilibration_timesteps = 170
        self.lmp_instance.run(equilibration_timesteps)

        #production
        production_MD_cmds = ['reset_timestep 0',
                              'fix 		3 solid  langevin 300.0 300.0 100.0 699483 gjf yes']
        self.lmp_execute_list_of_commands(production_MD_cmds)

        #output
        output_vars = {'group0': 'c_2',
                       'energy_': 'etotal',
                       'potential_e': 'pe',
                       'kinetic_e': 'ke',
                       'sim_temp': 'c_SolidTemp',
                       'time_step': 'step'}

        self.set_dictionary_as_lammps_variables(output_vars, 'equal')
        self.lmp_instance.command('fix 		5 all print 100 "${energy_} ${potential_e} ${kinetic_e} ${sim_temp} ${time_step}" file ${group_output} screen no')

    def lmp_execute_list_of_commands(self, list_of_commands):
        for line in list_of_commands:
            self.lmp_instance.command(line)

    def lmps_run(self, Id, parameters, paths, fifos={}, run_script=False):
        self._start_fifo_capture(fifos, Id)
        # submitting parameters
        self.set_dictionary_as_lammps_variables(parameters, 'string')
        # submitting paths
        self.set_dictionary_as_lammps_variables(paths, 'string')
        # submitting run Id
        self.lmp_variable('num', 'string','%05d' % Id)
        # use default settings
        self.default_settings()
        # starting script if desired
        if run_script:
            self.lmp_instance.file(paths['script'])
        # specify fifo dumps
        self.set_fifos(fifos)
        self.lmp_instance.command('run ${time_steps}')
        # write snapshot of end-comformation
        self.lmp_instance.command('write_dump solid xyz ${end_xyz}.xyz')
        self.lmp_instance.close()
