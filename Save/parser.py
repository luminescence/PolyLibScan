import json
import numpy as np
import gzip
import itertools as it

class Parser(object):

    def energy(self, file_name):
        with open(file_name) as f:
            data = f.read().split('\n')
        try:
            return np.array([map(float, line.split()) for line in data[1:-1]])
        except ValueError:
            raise ValueError('Error in File %s. Could not parse string to float' % file_name)

    def meta(self, config):
        meta = {}
        weights = config.sim_parameter['weights']
        monomers = config.sim_parameter['monomers']
        meta['weights'] = np.array([(mono, weight) for mono,weight in zip(monomers, weights)], 
                                                     dtype=[('monomers', '|S10'), ('weight', '<f4')])
        
        meta['misc'] = np.array([(key, str(value)) for key,value in config.sim_parameter.items()
                                 if (not isinstance(value, list)) and (not isinstance(value, dict))], 
                                 dtype=[('name', '|S15'), ('value', '|S22')])
        meta['parameter'] = np.array([(key, value) for key,value in config.lmp_parameter.items()
                                      if not isinstance(value, list)],
                                      dtype=[('name', '|S20'), ('value', '<i4')])
        meta['sequence'] = np.array([(id_, name) for name,id_ in zip(config.sim_parameter['named_sequence'], 
                                                                     config.sim_parameter['id_sequence']   )], 
                                                     dtype=[('ID', '>i2'), ('monomer', '|S10')])
        as_val = config.sim_parameter['active_site']
        meta['active_site'] = np.array([ (as_val['xyz'][i], as_val['chain'][i], as_val['pdb_id'][i], as_val['iCode'][i]) 
                                                                                    for i in xrange(len(as_val['xyz']))], 
                                        dtype=[('xyz', '<i2'), ('chain', '|S1'), ('pdb_id', '<i2'), ('iCode', '|S1')])
        return meta
        
    def xyz(self, file_path):
        with open(file_path) as f:
            data = f.read().split('\n')[:-1]
        try:
            atoms_count = int(data[0])
            timestep = int(data[1].split()[-1])
        except ValueError:
            raise ValueError('Faulty format found in file %s' % file_path)

        dtype = {'names':['atom_type', 'x', 'y', 'z'],
                 'formats': [np.int] + 3*[np.float]}
        atoms = np.zeros(atoms_count, dtype=dtype)
        
        try:
            for i, line in enumerate(data[2:]):
                items = line.split()
                atoms[i][0] = int(items[0])
                atoms[i][1] = float(items[1])
                atoms[i][2] = float(items[2])
                atoms[i][3] = float(items[3])
        except ValueError:
            raise ValueError('Faulty format found in file %s' % file_path)        
        
        return atoms

    def distance(self, file_path):
        with open(file_path) as f:
            data = f.read().split('\n')[:-1]
        dtype = [('time_step', np.int),('distance', np.float)]
        dist_data = np.zeros(len(data), dtype=dtype)
        for i, line in enumerate(data):
            items = line.split()
            dist_data[i] = (int(items[0]), float(items[1]))
        return dist_data

    def trajectory(self, file_path):
        with gzip.open(file_path) as f:
            try:
                lines = liner(f)
                chunks = chunk_data(lines)
                traj = np.array(list(it.chain.from_iterable([el[1] for el in chunks])))
            except IndexError:
                raise IndexError('Faulty format found in file %s' % file_path)
        return traj

    def version(self, version_info):
        np_format = [('program', '|S11'), ('version', '|S40')]
        return np.array(version_info.items(), dtype=np_format)


    def particle_list(self, path):
        '''Extract particles from list
        '''
        dtype = [('xyz', np.int), ('p_id', np.int), ('name', 'S6'), 
                 ('chain', 'S1'), ('atom', 'S6'), ('res_id', np.int), ('iCode', 'S1')]
        return np.fromfile(path, dtype=dtype)

    def trajectory_meta(self, file_path):
        '''Extract the list which specifies the 
        '''
        with gzip.open(file_path) as f:
            try:
                lines = liner(f)
                chunks = chunk_data(lines)
                first = chunks.next()
                second = chunks.next()
            except (IndexError, ValueError):
                raise ValueError('Faulty format found in file %s' % file_path)
            else:
                type_list = first[1]['type_']
                step_size = second[0] - first[0]
        return type_list, step_size

    def traj_info(self, file_path):
        type_list, step_size = self.trajectory_meta(file_path)
        data = np.array([('particle_number', len(type_list)),
                         ('step_size', step_size)], dtype=[('key', '|S15'),
                                                           ('value', np.int)])
        return data

def liner(handle):
    for raw_line in handle:
        line = raw_line[:-1]
        yield line

def chunk_data(data_stream):
    for line in data_stream:
        line_num = int(line)
        xyz = np.zeros(line_num, dtype=[('type_', np.int), ('xyz', np.float,3)])
        timestep_string = data_stream.next()
        time_step = int(timestep_string.split()[2])
        for i,d in enumerate(it.islice(data_stream, line_num)):
            val = d.split()
            xyz[i] = (int(val[0]), (float(val[1]), float(val[2]), float(val[3])))
        yield time_step, xyz
