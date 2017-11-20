from lammps import PyLammps


class LmpController(object):
    """Initialize and control a lammps instance"""

    def __init__(self, Id, parameters, paths, parameterisation, fifos={}, previous_instance='', stoichiometry=[1,1]):
        if isinstance(previous_instance, PyLammps):
            self.instance = previous_instance
        elif previous_instance != '':
            raise ValueError('The previous_instance passed to LmpController is not a PyLammps instance!')
        else:
            self.instance = PyLammps()

            self.Id = Id
            self.parameters = parameters
            self.paths = paths
            self.parameterisation = parameterisation
            self.fifos = fifos
            self.stoichiometry=stoichiometry

            self.is_protein_present = (self.stoichiometry[0] > 0)
            self.is_polymer_present = (self.stoichiometry[1] > 0)

    def get_lmp_styles(self):
        styles = OrderedDict()
        globals_ = self.parameterisation['globals']
        available_styles = filter(lambda x: x in globals_, 
            ['atom_style', 'bond_style', 'angle_style', 'dihedral_style'])
        for type_ in available_styles:
            styles[type_] = self.add_style(type_, globals_)
        styles['pair_style'] = self.add_pair_style(globals_, self.parameterisation['Pairs'])
        return styles

    @staticmethod
    def add_style(style_type, globals_):
        style = []
        if 'hybrid' in globals_[style_type]:
            sub_name = style_type.replace('_', '_sub')
            substyles = {key: val for key,val in globals_.items() if sub_name in key}
            if len(substyles) == 0:
                raise AttributeError("%s is hybrid but has no substyles." % style_type)
            style = [globals_[style_type]] + [v for k,v in sorted(substyles.items())]
        else:
            style = [globals_[style_type]]
        return style

    @staticmethod
    def add_pair_style(globals_, pairs):
        style_type = globals_['pair_style']
        pair_list = []
        if style_type in ['hybrid', 'hybrid/overlay']:
            pair_list += [style_type]
            sub_name = 'pair_substyle'
            substyles = {key: val for key,val in globals_.items() if sub_name in key}
            if len(substyles) == 0:
                raise AttributeError("pair_style is hybrid but has no substyles.")
            sub_styles = [v for k,v in sorted(substyles.items())]
            for sub_style in sub_styles:
                if pairs[sub_style]['kind'] in ['lj/cut', 'lj96/cut', 'coul/cut', 'coul/diel', 'soft']:
                    pair_list += [pairs[sub_style]['kind'], pairs[sub_style]['cutoff']]
                elif pairs[sub_style]['kind'] in ['coul/debye']:
                    pair_list += [pairs[sub_style]['kind'], pairs[sub_style]['kappa'], pairs[sub_style]['cutoff']]
                else:
                    raise KeyError('pair_style %s is currently not implemented.' % sub_style)
        else:
            pair_list += [pairs[style_type]['kind'], pairs[style_type]['cutoff']]
        return pair_list

    def set_dictionary_as_lammps_variables(self, dictionary, var_style):
        """dictionary keys will be variable names, dictionary values will be variable values"""
        for name, val in dictionary.items():
            # lists are converted to strings
            if isinstance(val, list):
                val = self.convert_python_list_to_lammps_list(val)
            self.variable(name, var_style, val)

    def variable(self, name, style, val):
        """
        >>> test_controller = controler()
        LAMMPS output is captured by PyLammps wrapper
        >>> test_controller.variable('a', 'string', '2')
        >>> test_controller.instance.variables['a'].value == 2.0
        True"""
        self.instance.variable('%s %s %s' % (name, style, val))

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
            self.instance.command(fifo.lammps_string())

    @staticmethod
    def _start_fifo_capture(fifos, index):
        for name, fifo in fifos.items():
            fifo.activate(index)

    def settings(self):
        '''Set styles and read data file into Lammps.'''
        initialization_cmds = ['units		real',
                               'boundary	p p p']
        for name,style in self.get_lmp_styles().items():
            initialization_cmds += ['%s	%s' % (name, ' '.join(map(str, style)))]
        self.execute_list_of_commands(initialization_cmds)

        #reading coordinate file
        self.instance.command('read_data	%s/%05d' % (self.paths['input'], self.Id))

    def minimization(self):
        # minimization
        self.instance.command('minimize 0.0 1.0e-8 10000 100000')

    def configure_output(self):
        # output & computes

        print_interval = 100 # every n timesteps...
        self.lmp_instance.command('thermo %s' % print_interval) # ...determine thermodynamics

        self.lmp_instance.command('compute SolidTemp solid temp')

        list_of_output_var_tuples = [('energy_', 'etotal'),
                                     ('potential_e', 'pe'),
                                     ('kinetic_e', 'ke'),
                                     ('sim_temp', 'c_SolidTemp'),
                                     ('time_step', 'step')]
        if self.is_protein_present and self.is_polymer_present:
            self.lmp_instance.command('compute 2 contacts group/group polymer pair yes')
            list_of_output_var_tuples = [('group0', 'c_2')] + list_of_output_var_tuples

        output_vars = OrderedDict(list_of_output_var_tuples)
        self.set_dictionary_as_lammps_variables(output_vars, 'equal')
        self.instance.command('fix 		5 all print 100 %s file %s/Energy%05d screen no' % (
        self.convert_python_list_to_lammps_list([self.variable_sigil_for_lammps(x) for x in output_vars.keys()]), self.paths['output'], self.Id))

    def production_MD(self, temperature_production):
        # production
        self.instance.command('reset_timestep 0')
        self.langevin_dynamics(T_0=temperature_production, T_end=temperature_production)
        self.instance.run(self.parameters['time_steps'])
        # write snapshot of end-comformation
        output_file_final_conformation = self.paths['output'] + '/trajectoryE' + '%05d' % self.Id
        self.instance.command('write_dump solid xyz %s.xyz' % output_file_final_conformation)

    def equilibration_MD(self, temperature_production, temperature_start):
        output_file_start_conformation = self.paths['output'] + '/trajectoryS' + '%05d' % self.Id

        equilibration_cmds = ['velocity 	solid create %s 1321 dist gaussian' % temperature_start,
                              'write_dump	solid xyz %s.xyz' % output_file_start_conformation]
        self.execute_list_of_commands(equilibration_cmds)
        self.langevin_dynamics(T_0=temperature_start, T_end=temperature_production)
        equilibration_timesteps = 170
        self.instance.run(equilibration_timesteps)

    def group_declaration(self):
        """define groups as needed, """
        groups = ['group solid union all']

        if self.is_protein_present:
            protein_exclusive_groups = ['group protein type %s' % self.parameters['bb_id'],
                                        'group ghost_protein type %s' % self.parameters['ghost_id'],
                                        'group activesite id %s' % self.convert_python_list_to_lammps_list(self.parameters['active_site_ids'], with_quotes=False),
                                        'group solid subtract all ghost_protein'] # overwrite default
            groups += protein_exclusive_groups

        if self.is_polymer_present:
            polymer_exclusive_groups = ['group polymer type %s' % self.convert_python_list_to_lammps_list(self.parameters['monomer_ids'], with_quotes=False)]
            groups += polymer_exclusive_groups

        if self.is_protein_present and self.is_polymer_present:
            protein_and_polymer_groups = ['group contacts subtract all protein ghost_protein polymer',
                                          'group interactions union polymer contacts',
                                          'group distance_group union polymer activesite']
            groups += protein_and_polymer_groups

        self.execute_list_of_commands(groups)

    def langevin_dynamics(self, T_0, T_end):
        cmds = ['fix 		3 solid  langevin %s %s 100.0 699483 gjf yes' % (T_0, T_end),
                'fix 		4 solid nve']
        self.execute_list_of_commands(cmds)

    def model_specifications(self):
        # specify the model
        pair_coeffs = ['timestep 	%s' % self.parameters['timestep'],
                       'dielectric %s' % self.parameters['dielectric_par']]
        for pair in [pair for key,pair in self.parameterisation['Pairs'].items() if not pair['single_type_parametrised']]:
            pair_parameter_string = ' '.join(map(str, sorted([p for k,p in 
                                                 sorted(pair.items()) if 'coef' in k])))
            if 'coul' in pair['kind']:
                pair_coeffs.append('pair_coeff * * %s' % (pair['kind']))
            else:
                pair_coeffs.append('pair_coeff * * %s %s' % (pair['kind'], 
                                                             pair_parameter_string))
        self.execute_list_of_commands(pair_coeffs)

    def execute_list_of_commands(self, list_of_commands):
        for line in list_of_commands:
            self.instance.command(line)

    def lmps_run(self, run_script=False):
        self._start_fifo_capture(self.fifos, self.Id)
        # submitting run Id
        self.variable('num', 'string','%05d' % self.Id)
        # use default settings
        self.settings()
        # starting script if desired
        if run_script:
            self.instance.file(self.paths['script'])
        self.model_specifications()
        self.group_declaration()
        #modify neighbor list
        if self.is_protein_present:
            self.instance.command('neigh_modify exclude group ghost_protein all')
        self.minimization()
        #equilibrate
        temperature_start = 100.0
        temperature_production = 300.0
        self.equilibration_MD(temperature_production, temperature_start)
        self.configure_output()
        # specify fifo dumps
        self.set_fifos(self.fifos)
        self.production_MD(temperature_production)
        self.instance.close()


if __name__ == "__main__":
    import doctest

    doctest.testmod()
