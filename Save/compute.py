import numpy as np

def distance_to_active_site(xyz_coords, polymer_ids, active_site_no):
    '''calculates the minimal distance between the polymer and the 
    active site.
    '''
    mask = np.in1d(xyz_coords['atom_type'], polymer_ids)
    poly_coords = xyz_coords[mask]
    # lammps atom ids start at 1 while numpy arrays
    # ids start at 0, thus decreasing the index 
    # yields the right index
    active_site = xyz_coords[active_site_no-1]
    min_dist = 1000
    for resi in active_site:
        for mono in poly_coords:
            md = _dist(resi, mono)
            if md < min_dist:
                min_dist = md
    return min_dist

def _dist(a, b):
    an = np.array([a['x'], a['y'], a['z']])
    bn = np.array([b['x'], b['y'], b['z']])
    c = an - bn
    return np.sqrt(np.sum(c**2))
