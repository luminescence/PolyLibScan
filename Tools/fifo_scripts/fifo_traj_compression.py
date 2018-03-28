import sys
import gzip
import shutil


def compress_stream(infile, outfile):
    with open(infile, 'r') as f_in, gzip.open(outfile, 'w') as f_out:
        f_out.write(f_in.read())

if __name__ == '__main__':
    traj_file = sys.argv[1]
    out_file = sys.argv[2]
    compress_stream(traj_file, out_file)