import re
import numpy as np


def read_data(path):
    with open(path) as f:
        return f.read().split('\n')


def find_line(data, string):
    for i,d in enumerate(data):
        if string in d:
            return i
    else:
        raise Exception(string + ' not in data')

def get_summary(data):
    summary_line = find_line(data, 'SUMMARY OF THIS PREDICTION')
    end_of_summary = find_line(data[summary_line:], '------------')
    return data[summary_line+2:summary_line+end_of_summary]

def parse_summary(summary):
    resn = '(?P<resn>[A-Z-+]+)'
    resi = '(?P<resi>\d{1,3})'
    chain = '(?P<chain>[A-Z])'
    pka = '(?P<pka>\d+\.\d+)'
    #std_pka = '(?P<std_pka>\d+\.\d+)'
    sum_line = re.compile(r'.*?%s.*?%s %s.*?%s.*' % (
    resn, resi, chain, pka))
    data = []
    for line in summary:
        try:
            data.append(re.match(sum_line, line).groupdict())
        except:
            raise Exception('Error in line: %s' % line)
    return data

def to_array(parsed_summary):
    dtype = [('resn', 'S3'), ('chain', 'S1'), ('resi', np.int), ('iCode', 'S1'), ('pka', np.float32)]
    ar = np.zeros(len(parsed_summary), dtype=dtype)
    for field in parsed_summary[0].keys():
        ar[field] = [r[field] for r in parsed_summary]
    ar['iCode'] = ' '
    return ar

def create_pipeline(pipeline):
    '''partial of through_pipeline.
    '''
    return lambda data: through_pipeline(data, pipeline)

def through_pipeline(data, pipeline):
    '''Pipeline implementation.
    '''
    return reduce(lambda x,f: f(x), pipeline, data)

def retrieve_propka_result(path):
    pipeline = create_pipeline([read_data, 
                                get_summary, 
                                parse_summary, 
                                to_array])
    return pipeline(path)

