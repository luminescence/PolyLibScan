import numpy as np

from collections import OrderedDict
from lammps import PyLammps
from PolyLibScan.Analysis.sim_run import AtomFilter


class LmpController(object):
    """Initialize and control a lammps instance"""

    def __init__(self, Id, parameters, paths, parameterisation, fifos, previous_instance=''):
        '''
        Control LAMMPS simulation.

        Input:
            Id: number of simulation - used for input-, output-name. [int]
            parameters: parameters from config.yml [dict]
            paths: paths found in config.yml [dict]
            parameterisation: parameters of beads and force fields, found in parameters.yml [dict]
            fifos: fifo information from config.yml [dict]
            previous_instance: Pylammps object
        '''
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
        self.stoichiometry = parameters['stoichiometry']
        self.pairs = PairList(self.parameterisation['globals'], self.parameterisation['Pairs'])

        self.is_protein_present = (self.stoichiometry[0] > 0)
        self.is_polymer_present = (self.stoichiometry[1] > 0)

    def get_lmp_styles(self):
        styles = OrderedDict()
        globals_ = self.parameterisation['globals']
        available_styles = filter(lambda x: x in globals_, 
            ['atom_style', 'bond_style', 'angle_style', 'dihedral_style'])
        for type_ in available_styles:
            styles[type_] = self.add_style(type_, globals_)
        styles['pair_style'] = [self.pairs.pair_style_str]
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

    def set_dictionary_as_lammps_variables(self, dictionary, var_style):
        """dictionary keys will be variable names, dictionary values will be variable values"""
        for name, val in dictionary.items():
            # lists are converted to strings
            if isinstance(val, list):
                val = self.convert_python_list_to_lammps_list(val)
            self.variable(name, var_style, val)

    def variable(self, name, style, val):
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

    def configure_output(self, print_interval=100):
        # output & computes

        self.instance.command('thermo %s' % print_interval) # ...determine thermodynamics

        self.instance.command('compute SolidTemp solid temp')

        list_of_output_var_tuples = [('energy_', 'etotal'),
                                     ('potential_e', 'pe'),
                                     ('kinetic_e', 'ke'),
                                     ('sim_temp', 'c_SolidTemp'),
                                     ('time_step', 'step')]
        if self.is_protein_present and self.is_polymer_present:
            self.instance.command('compute 2 contacts group/group polymer pair yes')
            list_of_output_var_tuples = [('group0', 'c_2')] + list_of_output_var_tuples

        output_vars = OrderedDict(list_of_output_var_tuples)
        self.set_dictionary_as_lammps_variables(output_vars, 'equal')
        
        output_file_path = '%s/Energy%05d' % (self.paths['output'], self.Id)
        csv_header = self.convert_python_list_to_lammps_list(output_vars.keys())

        formatted_output_vars = self.convert_python_list_to_lammps_list([self.variable_sigil_for_lammps(x)
                                                                         for x in output_vars.keys()])

        self.instance.command('fix 5 all print %s %s file %s screen no title %s' % (print_interval,
                                                                                    formatted_output_vars,
                                                                                    output_file_path,
                                                                                    csv_header))

    def production_MD(self, temperature_production, time_bug_fix=True):
        # production
        self.instance.command('reset_timestep 0')
        self.langevin_dynamics(T_0=temperature_production, T_end=temperature_production)
        # on Niklas's local machine, LAMMPS stops after 86200 steps
        # one can run 80k steps n-times, though!
        if time_bug_fix:
            max_no_time_steps = 80000
            full_runs, timesteps_left = divmod(self.parameters['time_steps'], max_no_time_steps)

            for x in range(full_runs):
                self.instance.run(max_no_time_steps)
            self.instance.run(timesteps_left)
        else:
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
        equilibration_timesteps = 1700
        self.instance.run(equilibration_timesteps)

    def group_declaration(self):
        """define groups as needed, """
        groups = ['group solid union all']

        if self.is_protein_present:

            active_site_ids_lammps_list = self.convert_python_list_to_lammps_list(self.parameters['active_site_ids'],
                                                                                  with_quotes=False)
            protein_exclusive_groups = ['group protein type %s' % self.parameters['bb_id'],
                                        'group ghost_protein type %s' % self.parameters['ghost_id'],
                                        'group activesite id %s' % active_site_ids_lammps_list,
                                        'group solid subtract all ghost_protein']   # overwrite default
            groups += protein_exclusive_groups

        if self.is_polymer_present:
            type_filter = AtomFilter(self.parameters['particle_ids'],
                                     self.parameters['poly_sequence'],
                                     monomer_id=self.parameters['monomer_ids'],
                                     molecule='polymer',
                                     filter_specification='type')
            polymer_bead_ids = np.where(type_filter.mask == True)[0].tolist()
            polymer_bead_ids = [x+1 for x in polymer_bead_ids]  # increment by one: LAMMPS starts counting at 1,
                                                                # python at 0

            polymer_ids_lammps_list = self.convert_python_list_to_lammps_list(polymer_bead_ids, with_quotes=False)
            polymer_exclusive_groups = ['group polymer id %s' % polymer_ids_lammps_list]
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
        for style in self.pairs.single_parameterized():
            pair_coeffs.append(self.pairs.pair_coef_str(style))
        self.execute_list_of_commands(pair_coeffs)

    def execute_list_of_commands(self, list_of_commands):
        for line in list_of_commands:
            self.instance.command(line)

    def lmps_run(self, run_script=False):
        self._start_fifo_capture(self.fifos, self.Id)
        # submitting run Id
        self.variable('num', 'string', '%05d' % self.Id)
        # use default settings
        self.settings()
        # starting script if desired
        if run_script:
            self.instance.file(self.paths['script'])
        self.model_specifications()
        self.group_declaration()
        # modify neighbor list
        if self.is_protein_present:
            self.instance.command('neigh_modify exclude group ghost_protein all')
        self.minimization()
        # equilibrate
        temperature_start = 100.0
        temperature_production = 300.0
        self.equilibration_MD(temperature_production, temperature_start)
        if 'output_interval' in self.parameters:
            self.configure_output(self.parameters['output_interval'])
        else:
            self.configure_output()
        # specify fifo dumps
        self.set_fifos(self.fifos)
        self.production_MD(temperature_production)
        self.instance.close()

class PairList(object):
    
    def __init__(self, globals_, pairs):
        self.has_substyles = False
        self.pairs = pairs
        self.styles = self.set_styles(globals_)
        
    def set_styles(self, styles):
        style_type = {'main': styles['pair_style']}
        if style_type['main'] in ['hybrid', 'hybrid/overlay']:
            style_type['sub'] = self.get_substyles(styles)
        else:
            style_type['sub'] = []
        return style_type
        
    def get_substyles(self, styles):
        sub_name = 'pair_substyle'
        substyles = {key: val for key,val in styles.items() if sub_name in key}
        if len(substyles) == 0:
            raise AttributeError("pair_style is hybrid but has no substyles.")
        self.has_substyles = True
        return [v for k,v in sorted(substyles.items())]
    
    @property
    def pair_style_str(self):
        output = self.pair_string(self.styles['main'])
        if self.has_substyles:
            for style in self.styles['sub']:
                output += self.pair_string(style)
        return ' '.join(map(str, output))
        
    def single_parameterized(self):
        if self.has_substyles:
            return [style for style in self.styles['sub'] if self.pairs[style]['single_type_parametrised']]
        elif self.pairs[self.styles['main']]['single_type_parametrised']:
            return [self.styles['main']]
        else:
            return []
        
    def pair_string(self, style):
        if style in ['hybrid', 'hybrid/overlay']:
            return [style]
        else:
            if self.pairs[style]['kind'] in ['lj/cut', 'lj96/cut', 'coul/cut', 'coul/diel', 'soft']:
                return [self.pairs[style]['kind'], self.pairs[style]['cutoff']]
            elif self.pairs[style]['kind'] in ['coul/debye']:
                return [self.pairs[style]['kind'], self.pairs[style]['kappa'], self.pairs[style]['cutoff']]
            else:
                raise KeyError('pair_style %s is currently not implemented.' % style)
    
    def pair_coef_str(self, style, atom_type1='*', atom_type2='*'):
        data = {}
        data['a1'] = str(atom_type1)
        data['a2'] = str(atom_type2)
        data['style'] = self.pairs[style]['kind']
        data['cutoff'] = self.pairs[style]['cutoff']

        coefficients = sorted([p for k,p in self.pairs[style].items() if 'coef' in k])
        if len(coefficients) > 0:
            data['parameters'] = ' '.join(map(str, coefficients))
        else:
            data['parameters'] = ''
        if self.has_substyles:
            format_str = 'pair_coeff {a1} {a2} {style} {parameters} {cutoff}'
        else:
            format_str = 'pair_coeff {a1} {a2} {parameters} {cutoff}'
        
        return format_str.format(**data)
