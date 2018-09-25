import numpy as np

def confidence_interval(bootstrap_data, confidence_level):
    # make sure confidence_level has sane value
    if not 0.0 < confidence_level < 1.0:
        raise ValueError('confidence_level must be in [0,1]')
    non_zero_data = bootstrap_data[~np.isnan(bootstrap_data)]
    data_length = non_zero_data.shape[0]
    if data_length == 0:
        return (np.nan, np.nan, np.nan)
    elif data_length == 1:
        return (data_length[0], np.nan, np.nan)
    alpha = (1.0-confidence_level)/2.0
    sorted_data = np.sort(non_zero_data)
    lower_limit = int(       alpha  * data_length)
    upper_limit = int((1.0 - alpha) * data_length)
    
    #results
    mean_ = np.nanmean(non_zero_data)
    lower_boundary = sorted_data[lower_limit+1]
    upper_boundary = sorted_data[upper_limit]
    return (mean_, lower_boundary, upper_boundary)

def Bootstrap(values, iterations=100):
    '''Takes a data array and build up subsamples of size
    sampling_size. The mean of each subsample is calculated. 
    The number of subsamples is set by the argument resampling_size.
    In Case an array with no values is given, an array is returned 
    with the sampling size
    '''
    if values.shape[0] == 0:
        return np.array(2*[np.nan])
    bootstraps = np.zeros((iterations, values.shape[0]))
    for i in xrange(iterations):
        bootstraps[i] = np.random.choice(values, size=values.shape[0])
    return bootstraps