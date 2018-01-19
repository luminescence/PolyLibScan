import numpy as np
import numba as nb

from PolyLibScan.Analysis.sim_run import AtomFilter


def distance_to_active_site(xyz_coords, db, polymer_ids, active_site_no):
    '''calculates the minimal distance between the polymer and the 
    active site.
    '''

    type_filter = AtomFilter(db.traj_type_order, db, polymer_ids, molecule='polymer')

    poly_coords = xyz_coords[type_filter.mask]
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

@nb.jit
def _dist(a, b):
    an = np.array([a['x'], a['y'], a['z']])
    bn = np.array([b['x'], b['y'], b['z']])
    c = an - bn
    return np.sqrt(np.sum(c**2))
