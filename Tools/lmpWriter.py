import jinja2 as jga
import os
import itertools as it

class LmpWriter(object):
    '''Writes the pdb data in a lammps compatible data format.
    Additional parameters are specified in the config file.
    '''
    def __init__(self, environment):
        template_path = os.path.dirname(os.path.abspath(__file__))+'/lmpsTemplate.tpl'
        self.template = jga.Template(open(template_path).read())
        self.env = environment
        self.env.coefficient_string = self.transfer_coeff_str(self.env)
        self.env.particle_str = self.transfer_particle_str(self.env)
        self.env.calc_box(add_margin=True)
        self.entryList = {'bond': ['bond', 'angle', 'full', 'molecular'],
                          'angle': ['angle', 'full', 'molecular'],
                          'dihedral': ['full', 'molecular']}

    def chain_data(self, key):
        molecules = sorted(self.env.molecules.values(), key=lambda x: x.Id)
        return it.chain(*[molecule.data[key] for molecule in molecules])

    def count_number(self, kind):
        return len(list(self.chain_data(kind)))

    def pair_list(self, pair_obj):
        sorted_atom_types = sorted(self.env.atom_type.values(), key=lambda x:x.Id)
        return [pair_obj[(atomType1, atomType2)] for i, atomType1 in enumerate(sorted_atom_types)
            for atomType2 in sorted_atom_types[i:] ]

    def write(self, out_path):
        with open(out_path, 'w') as f:
            f.write(self.template.render(env=self.env, writer=self))

    def entry_is_needed(self, entry):
        if self.env.globals['atom_style'] in self.entryList[entry]:
            return True
        elif self.env.globals['atom_style'] == 'hybrid':
            if self.env.globals['atom_substyle1'] in self.entryList[entry]:
                return True
            elif self.env.globals['atom_substyle2'] in self.entryList[entry]:
                return True
        return False

    def pair_string(self, pair):
        if self.env.globals['pair_style'] == 'hybrid':
            if pair.pair_type.kind == 'lj/cut':
                bond_template = '{:> 3d} {:> 3d} {} {:> 6.2f} {:> 6.2f} {:> 6.2f}'
                return bond_template.format(
                            pair.atom_type2.Id, pair.atom_type1.Id, 
                            pair.pair_type.kind, pair.epsilon, pair.sigma, 
                            pair.cutoff)
            elif pair.pair_type.kind == 'lj96/cut':
                bond_template = '{:> 3d} {:> 3d} {} {:> 6.2f} {:> 6.2f} {:> 6.2f}'
                return bond_template.format(
                            pair.atom_type2.Id, pair.atom_type1.Id, 
                            pair.pair_type.kind, pair.epsilon, pair.sigma, 
                            pair.cutoff)
            elif pair.pair_type.kind == 'morse':
                bond_template = '{:> 3d} {:> 3d} {} {:> 6.2f} {:> 6.2f} {:> 6.2f} {:> 6.2f}'
                return bond_template.format(
                            pair.atom_type2.Id, pair.atom_type1.Id, 
                            pair.pair_type.kind, pair.epsilon, pair.alpha, 
                            pair.sigma, pair.cutoff)
            elif pair.pair_type.kind == 'coul/diel':
                bond_template = '{:> 3d} {:> 3d} {} {:> 6.2f} {:> 6.2f} {:> 6.2f}'
                return bond_template.format(
                            pair.atom_type2.Id, pair.atom_type1.Id, 
                            pair.pair_type.kind, pair.pair_type.parameters['coeffs'][0], 
                            pair.pair_type.parameters['coeffs'][1], pair.pair_type.parameters['coeffs'][2])
            elif pair.pair_type.kind == 'coul/cut':
                bond_template = '{:> 3d} {:> 3d} {}'
                return bond_template.format(
                            pair.atom_type2.Id, pair.atom_type1.Id, 
                            pair.pair_type.kind)
            elif pair.pair_type.kind == 'soft':
                bond_template = '{:> 3d} {:> 3d} {} {:> 6.2f}'
                #The following coefficients must be defined for each pair of atom types 
                # via the pair_coeff command as in the examples above, or in the data file 
                # or restart files read by the read_data or read_restart commands, or by 
                # mixing as described below:
                #
                #     A (energy units)
                #     cutoff (distance units)
                #
                # The last coefficient is optional. If not specified, the global soft cutoff is used.
                return bond_template.format(
                            pair.atom_type2.Id, pair.atom_type1.Id, 
                            pair.pair_type.kind, pair.epsilon)
            else:
                raise Exception('Kind of pair style not implemented/known.')
        else:
            if pair.pair_type.kind == 'lj/cut':
                bond_template = '{:> 3d} {:> 3d} {:> 6.2f} {:> 6.2f} {:> 6.2f}'
                return bond_template.format(
                            pair.atom_type2.Id, pair.atom_type1.Id, 
                            pair.epsilon, pair.sigma, 
                            pair.cutoff)
            elif pair.pair_type.kind == 'lj96/cut':
                bond_template = '{:> 3d} {:> 3d} {} {:> 6.2f} {:> 6.2f} {:> 6.2f}'
                return bond_template.format(
                            pair.atom_type2.Id, pair.atom_type1.Id, 
                            pair.pair_type.kind, pair.epsilon, pair.sigma, 
                            pair.cutoff)

            elif pair.pair_type.kind == 'morse':
                bond_template = '{:> 3d} {:> 3d} {:> 6.2f} {:> 6.2f} {:> 6.2f} {:> 6.2f}'
                return bond_template.format(
                            pair.atom_type2.Id, pair.atom_type1.Id, 
                            pair.epsilon, pair.alpha, 
                            pair.sigma, pair.cutoff)
            elif pair.pair_type.kind == 'coul/diel':
                bond_template = '{:> 3d} {:> 3d} {:> 6.2f} {:> 6.2f} {:> 6.2f}'
                return bond_template.format(
                            pair.atom_type2.Id, pair.atom_type1.Id, 
                            pair.pair_type.parameters['coeffs'][0], 
                            pair.pair_type.parameters['coeffs'][1], pair.pair_type.parameters['coeffs'][2])
            elif pair.pair_type.kind == 'coul/cut':
                bond_template = '{:> 3d} {:> 3d}'
                return bond_template.format(
                            pair.atom_type2.Id, pair.atom_type1.Id)
            elif pair.pair_type.kind == 'soft':
                bond_template = '{:> 3d} {:> 3d} {:> 6.2f}'
                return bond_template.format(
                            pair.atom_type2.Id, pair.atom_type1.Id, 
                            pair.epsilon)
            else:
                raise Exception('Kind of pair style not implemented/known.')
        return pair_string

        
    def transfer_coeff_str(self, molecule):
        def coefficient_string(style_type):
            if molecule.globals['bond_style'] == 'hybrid':
                bond_template = '{:> 3d} {:10s}'+ '{:> 6.2f}'*len(style_type.parameters['coeffs'])
                return bond_template.format(style_type.Id, style_type.kind, *style_type.parameters['coeffs'])
            else:
                bond_template = '{:> 3d} '+ '{:> 7.2f}'*len(style_type.parameters['coeffs'])
            return bond_template.format(style_type.Id, *style_type.parameters['coeffs'])
        return coefficient_string


    def transfer_particle_str(self, molecule):
        sub_styles = filter(lambda x:'substyle' in x, molecule.globals.keys())
        sub_styles.sort(key=lambda x: [x[-1]])
        def particle_str(particle):
            if molecule.globals['atom_style'] == 'atomic':
                return '{:> 7d}{:> 4d}{:> 10.3f}{:> 10.3f}{:> 10.3f}'.format(
                    particle.Id, particle.type_.Id, *particle.position)
            elif molecule.globals['atom_style'] == 'bond':
                return '{:> 7d}{:> 4d}{:> 4d}{:> 10.3f}{:> 10.3f}{:> 10.3f}'.format(
                    particle.Id, particle.mol_id, particle.type_.Id, *particle.position)
            elif molecule.globals['atom_style'] == 'angle':
                return '{:> 7d}{:> 4d}{:> 4d}{:> 10.3f}{:> 10.3f}{:> 10.3f}'.format(
                    particle.Id, particle.mol_id, particle.type_.Id, *particle.position)
            elif molecule.globals['atom_style'] == 'hybrid':
                if molecule.globals['atom_substyle1'] == 'angle' and molecule.globals['atom_substyle2'] == 'charge':
                    return '{:> 7d}{:> 4d}{:> 10.3f}{:> 10.3f}{:> 10.3f}{:> 4d}{:> 6.2f}'.format(
                        particle.Id, particle.type_.Id, 
                        particle.position[0], particle.position[1], particle.position[2], 
                        particle.mol_id, particle.type_.charge)
                elif molecule.globals['atom_substyle1'] == 'charge' and molecule.globals['atom_substyle2'] == 'angle':
                    return '{:> 7d}{:> 4d}{:> 10.3f}{:> 10.3f}{:> 10.3f}{:> 6.2f}{:> 4d}'.format(
                        particle.Id, particle.type_.Id, 
                        particle.position[0], particle.position[1], particle.position[2], 
                        particle.type_.charge, particle.mol_id)
                else:
                    raise Exception('Something is wrong here... style: %s, substyle1: %s, substyle2: %s'%(
                        molecule.globals['atom_style'], molecule.globals['atom_substyle1'], molecule.globals['atom_substyle2']))
            else:
                raise Exception('Only the "atomic", "bond" and "angle" styles is supported for now. Edit your config file accordingly.')
        return particle_str