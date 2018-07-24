import numpy as np

def confidence_interval(bootstrap_data, confidence_level):
    # make sure confidence_level has sane value
    if confidence_level<0.0 or confidence_level>1.0:
        raise ValueError('Confidence Level must be in [0,1]')

    data_length = bootstrap_data.shape[0]
    alpha = (1.0-confidence_level)/2.0
    sorted_data = np.sort(bootstrap_data)
    lower_limit = int(       alpha  * data_length)
    upper_limit = int((1.0 - alpha) * data_length)
    
    #results
    mean_ = np.mean(bootstrap_data)
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
    means = np.zeros(iterations)
    for i in xrange(iterations):
        means[i] = np.random.choice(values, size=values.shape[0]).mean()
    return means
