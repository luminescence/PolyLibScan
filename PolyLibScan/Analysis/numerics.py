import numpy as np
import pandas as pd
import numba as nb
import PolyLibScan.helpers.bootstrap as bstrap

nb.jit()
def binning(data, energy, bins):
    '''Generate a histogramm for the frequency and
    summed energy. The results are stored in the bins
    array.

    input:
        data: distance_data
        energy: energy of the distance data
        bins: numpy array [m, 2]
    '''
    # resolution of bins are 0.1 Angstrom
    discrete_data = (10*data).astype(np.int)
    for i, d in enumerate(discrete_data[1:]):
        bins[d][0] += 1
        bins[d][1] += energy[i]

def cdf(array, dist):
    '''Returns the cumulative distance 
    probability in angstrom (binsize is 0.1 angstrom).
    Note that the array range is reduced by one.
    '''
    prob = 0.0
    # dist in angstrom conversion to bins
    max_bin = int(dist*10) - 1
    if max_bin < 0:
        return 0.0
    return np.sum(array[:max_bin])

nb.jit()
def cumulative_bins(bins):
    c_bins = np.zeros_like(bins)
    c_bins[0] = bins[0]
    for i, prob in enumerate(bins[1:], 1):
        c_bins[i] = c_bins[i-1] + prob
    return c_bins

def discard_tailing_zeros(vector):
    '''discards all zero-element at the end 
    of an array.

    Input:
        vector: numpy array (type integer)

    output:
        reduced vector: numpy array
    '''
    if vector.dtype.kind != 'i':
        raise ValueError("Vector must be of type integer!")
    for i,val in enumerate(vector[::-1]):
            if val != 0:
                max_idx = vector.shape[0] - i
                break
    return vector[:max_idx]

def Binding_ratio(distance_vector):
    """Calculate the fraction of binding over non-binding simulations.
    The distance vector should mark non-binding sims as np.nan (with the help of pandas).

    input:
    distance_vector -- numpy array

    output:
    float
    """
    return float(np.sum(~np.isnan(distance_vector)))/distance_vector.shape[0]


def mean_and_error(values):
    mean = values.mean()
    min_error = mean - values.min()
    max_error = values.max() - mean
    return mean, min_error, max_error


def Four_Fractions(vector):
    # guard against too little data
    if vector.shape[0] == 0:
        return np.array(4*[np.nan])
    elif vector.shape[0]<4:
        divisor = vector.shape[0]
    else:
        divisor = int(np.floor(vector.shape[0]/4))
    
    v = np.random.permutation(vector)
    
    binding_fractions = np.zeros(4)
    for i, n_min, n_max in zip(xrange(4), xrange(0, v.shape[0], divisor), xrange(divisor, v.shape[0]+1, divisor)):
        sub_sample = v[n_min:n_max]
        binding_fractions[i] = Binding_ratio(sub_sample)
    return binding_fractions


def distance_with_error(distance_matrix):
    results = np.zeros(distance_matrix.columns.shape[0], dtype=[('poly_name', '|S20'),
                                                                ('dist_mean', np.float),
                                                                ('dist_min_error', np.float),
                                                                ('dist_max_error', np.float)])
    for i, poly_name in enumerate(distance_matrix.columns):
        sub_results = Four_Fractions(distance_matrix[poly_name])
        stats = mean_and_error(sub_results)
        results[i] = tuple([poly_name]) + stats
    return pd.DataFrame(index=results['poly_name'], data=results[['dist_mean',
                                                                 'dist_min_error',
                                                                 'dist_max_error']])


def energy_with_error(energy_matrix, confidence_level=0.95):
    results = np.zeros(energy_matrix.columns.shape[0], dtype=[('poly_name', '|S20'),
                                                                ('energy_mean', np.float),
                                                                ('energy_min_error', np.float),
                                                                ('energy_max_error', np.float)])
    for i, poly_name in enumerate(energy_matrix.columns):
        mean_ = energy_matrix[poly_name].mean()
        mean_bs, lower_boundary, upper_boundary = bstrap.confidence_interval(bstrap.Bootstrap(energy_matrix[poly_name].dropna()), confidence_level)
        error_min = mean_ - lower_boundary
        error_max = upper_boundary - mean_
        results[i] = tuple([poly_name, mean_, error_min, error_max])
    return pd.DataFrame(index=results['poly_name'], data=results[['energy_mean',
                                                                 'energy_min_error',
                                                                 'energy_max_error']])

@nb.jit
def calc_box(coords, margin=20):
    """calculate the minimal box of an array of 3 dimensional coordinates and 
    add a margin to that minimal box.
    """
    box = np.zeros((2,3))
    for i in xrange(3):
        box[0,i] = coords[:,i].min()-margin
        box[1,i] = coords[:,i].max()+margin
    return box

@nb.jit
def in_box(coord, box):
    """check if coordinates are inside box.
    """
    if not np.any((coord - box[0])<0.0):
        if not np.any((coord - box[1])>0.0):
            return True
    return False