from lammps import PyLammps


class lmp_controler(object):
    """Initialize and control a lammps instance"""

    def set_dictionary_as_lammps_variables(self, dictionary, var_style, lmp_instance):
        for name, val in dictionary.items():
            # lists are converted to strings
            if isinstance(val, list):
                val = lmp_controler.convert_python_list_to_lammps_list(val)
            self.lmp_variable(lmp_instance, name, var_style, val)

    @staticmethod
    def lmp_variable(lmp_instance, name, style, val):
        lmp_instance.variable('%s %s %s' % (name, style, val))

    @staticmethod
    def convert_python_list_to_lammps_list(val):
        """convert from [a,b,c] to "a b c" which is required in lammps"""
        val = ' '.join(map(str, val))
        val = '"' + val + '"'
        return val

    @staticmethod
    def set_fifos(fifos, lmp_instance):
        for name,fifo in fifos.items():
            lmp_instance.command(fifo.lammps_string())

    @staticmethod
    def _start_fifo_capture(fifos, index):
        for name, fifo in fifos.items():
            fifo.activate(index)

    def default_settings(self, lmp_instance):

        #variables
        file_paths = {'dat_file': '${input}/${num}',
                      'start_xyz': '${output}/trajectoryS${num}',
                      'end_xyz': '${output}/trajectoryE${num}',
                      'group_output': '${output}/Energy${num}'}
        self.set_dictionary_as_lammps_variables(file_paths, 'string', lmp_instance)

        #initialization
        initialization_cmds = ['units		real',
                               'timestep 	${timestep}',
                               'boundary	p p p',
                               'atom_style	hybrid angle charge',
                               'bond_style  harmonic',
                               'angle_style harmonic',
                               'pair_style hybrid/overlay lj96/cut 10.5 coul/debye ${debye_kappa} 25.0',
                               'dielectric ${dielectric_par}']
        self.lmp_execute_list_of_commands(initialization_cmds, lmp_instance)

        #reading coordinate file
        lmp_instance.command('read_data	${dat_file}')
        pair_coeffs = ['pair_coeff * * lj96/cut 0.14   4.00   8.50',
                       'pair_coeff * * coul/debye']
        self.lmp_execute_list_of_commands(pair_coeffs, lmp_instance)

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
        self.lmp_execute_list_of_commands(groups, lmp_instance)

        #modify neighbor list
        lmp_instance.command('neigh_modify exclude group ghost_protein all')

        #computes
        lmp_instance.command('compute 2 contacts group/group polymer pair yes')
        lmp_instance.command('compute SolidTemp solid temp')

        # minimization
        lmp_instance.command('minimize 0.0 1.0e-8 100 1000')

        #equilibrate
        equilibration_cmds = ['velocity 	solid create 100.0 1321 dist gaussian',
                              'fix 		3 solid  langevin 1.0 300.0 100.0 699483 gjf yes',
                              'fix 		4 solid nve',
                              'write_dump	solid xyz ${start_xyz}.xyz']
        self.lmp_execute_list_of_commands(equilibration_cmds,lmp_instance)
        equilibration_timesteps = 170
        lmp_instance.run(equilibration_timesteps)

        #production
        production_MD_cmds = ['reset_timestep 0',
                              'fix 		3 solid  langevin 300.0 300.0 100.0 699483 gjf yes']
        self.lmp_execute_list_of_commands(production_MD_cmds, lmp_instance)

        #output
        output_vars = {'group0': 'c_2',
                       'energy_': 'etotal',
                       'potential_e': 'pe',
                       'kinetic_e': 'ke',
                       'sim_temp': 'c_SolidTemp',
                       'time_step': 'step'}

        self.set_dictionary_as_lammps_variables(output_vars, 'equal', lmp_instance)
        lmp_instance.command('fix 		5 all print 100 "${energy_} ${potential_e} ${kinetic_e} ${sim_temp} ${time_step}" file ${group_output} screen no')

    def lmp_execute_list_of_commands(self, list_of_commands, lmp_instance):
        for line in list_of_commands:
            lmp_instance.command(line)

    def lmps_run(self, Id, parameters, paths, fifos={}, run_script=False):
        self._start_fifo_capture(fifos, Id)
        # lammps_sim = self.lmp_instance
        lammps_sim = PyLammps()
        # submitting parameters
        self.set_dictionary_as_lammps_variables(parameters, 'string', lammps_sim)
        # submitting paths
        self.set_dictionary_as_lammps_variables(paths, 'string', lammps_sim)
        # submitting run Id
        self.lmp_variable(lammps_sim, 'num', 'string','%05d' % Id)
        # use default settings
        self.default_settings(lammps_sim)
        # starting script if desired
        if run_script:
            lammps_sim.file(paths['script'])
        # specify fifo dumps
        self.set_fifos(fifos, lammps_sim)
        lammps_sim.command('run ${time_steps}')
        # write snapshot of end-comformation
        lammps_sim.command('write_dump solid xyz ${end_xyz}.xyz')
        lammps_sim.close()