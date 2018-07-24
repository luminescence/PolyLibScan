import sys
import numpy as np
import itertools as it
from PolyLibScan.Analysis.sim_run import AtomFilter

# Functions

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

def group(chunks, polymer_ids):
    for time_step,xyz in chunks:
        type_filter = AtomFilter(xyz['type_'], polymer_sequence, polymer_ids, molecule='polymer', filter_specification='type')
        group1 = xyz[type_filter.mask]
        group2 = xyz[~type_filter.mask]
        yield time_step, group1['xyz'], group2['xyz']

def dist(groups):
    '''calculates the distance.
    '''
    for time_step, group1, group2 in groups:
        distance = 1000
        for coords1 in group1:
            for coords2 in group2:
                d = np.linalg.norm(coords1-coords2)
                if d < distance:
                    distance = d
        yield time_step, distance

def write_gen(handle, generator):
	for g in generator:
		line = ' '.join(map(str, g))
		handle.write(line+'\n')

def main(fifo_file, out_path, polymer_ids):
    with open(fifo_file) as f_in:
    	with open(out_path, 'w') as f_out:
    		lines = liner(f_in)
    		chunks =  chunk_data(lines)
    		groups = group(chunks, polymer_ids)
    		results = dist(groups)
    		write_gen(f_out, results)


if __name__ == '__main__':
    fifo_file = sys.argv[1]
    out_path = sys.argv[2]
    polymer_ids = map(int, sys.argv[3].split('-'))
    polymer_sequence = map(int, sys.argv[4].split('-'))

    main(fifo_file, out_path, polymer_ids)