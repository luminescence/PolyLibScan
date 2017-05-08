import numpy as np
import math

def Spherical_np(xyz):
    '''Returns all given coordinates into 
    spherical coordinates. 
    
    Input:
        xyz: numpy array [x, y, z]
    Output:
        rpt: numpy array [radius, phi, theta]

    In order to take advantage of serialization, an zeroed
    numpy array is added to the input in order to form a 
    ndarray that supports serialization operations.
    The zeroed vector is stripped on output.    '''
    xyz = np.vstack((xyz, np.ones((1,3))))
    rpt = np.zeros_like(xyz)
    xy = xyz[:,0]**2 + xyz[:,1]**2
    # radius
    rpt[:,0] = np.sqrt(xy + xyz[:,2]**2)
    # phi
    rpt[:,1] = np.arctan2(xyz[:,1], xyz[:,0])
    # theta
    rpt[:,2] = np.arctan2(np.sqrt(xy), xyz[:,2])
    return rpt[:-1]

def Cartesian_np(rpt):
    '''Transforms spherical coordinates to 
    cartesian coordinates.
    Input:
        rpt: numpy array [radius, phi, theta]
    Output:
        xyz: numpy array [x, y, z]

    In order to take advantage of serialization, an zeroed
    numpy array is added to the input in order to form a 
    ndarray that supports serialization operations.
    The zeroed vector is stripped on output.
    '''
    rpt = np.vstack((rpt, np.ones((1,3))))
    xyz = np.zeros_like(rpt)
    xyz[:,0] = rpt[:,0] * np.cos(rpt[:, 1]) * np.sin(rpt[:, 2])
    xyz[:,1] = rpt[:,0] * np.sin(rpt[:, 1]) * np.sin(rpt[:, 2])
    xyz[:,2] = rpt[:,0] * np.cos(rpt[:, 2])
    return xyz[:-1]


def rotation_matrix(axis, theta):
    """
    Return the rotation matrix associated with counterclockwise rotation about
    the given axis by theta radians.
    """
    axis = np.asarray(axis)
    theta = np.asarray(theta)
    axis = axis/math.sqrt(np.dot(axis, axis))
    a = math.cos(theta/2)
    b, c, d = -axis*math.sin(theta/2)
    aa, bb, cc, dd = a*a, b*b, c*c, d*d
    bc, ad, ac, ab, bd, cd = b*c, a*d, a*c, a*b, b*d, c*d
    return np.array([[aa+bb-cc-dd, 2*(bc+ad), 2*(bd-ac)],
                     [2*(bc-ad), aa+cc-bb-dd, 2*(cd+ab)],
                     [2*(bd+ac), 2*(cd-ab), aa+dd-bb-cc]])


def create_default_conf():
    config = ConfigParser.RawConfigParser(allow_no_value=True)

    # When adding sections or items, add them in the reverse order of
    # how you want them to be displayed in the actual file.
    # In addition, please note that using RawConfigParser's and the raw
    # mode of ConfigParser's respective set functions, you can assign
    # non-string values to keys internally, but will receive an error
    # when attempting to write to a file or when you get it in non-raw
    # mode. SafeConfigParser does not allow such assignments to take place.
    config.add_section('globals')
    config.set('globals', 'input_path', "")
    config.set('globals', 'output_path', "")
    config.set('globals', 'atom_style', 'atomic')


    config.add_section('Atom%d'% 1)
    config.set('Atom%d'% 1, 'name', 'CA')
    config.set('Atom%d'% 1, 'id', 1)
    config.set('Atom%d'% 1, 'mass', 12.0)
    # Writing our configuration file to 'example.cfg'
    with open('example.cfg', 'wb') as configfile:
        config.write(configfile)


def normalize(v):
    ar = np.array(v).astype(float)
    norm=ar.sum()
    if norm==0: 
        return ar
    return ar/norm

def normal(v1, v2):
    n = np.zeros(3)
    n[0] = v1[1] * v2[2] - v2[1] * v1[2]
    n[1] = v1[2] * v2[0] - v2[2] * v1[0]
    n[2] = v1[0] * v2[1] - v2[0] * v1[1]
    return n/np.linalg.norm(n)


def unit_vector(vector):
    """ Returns the unit vector of the vector.  """
    return vector / np.linalg.norm(vector)

def angle_between(v1, v2):
    """ Returns the angle in radians between vectors 'v1' and 'v2'::

            >>> angle_between((1, 0, 0), (0, 1, 0))
            1.5707963267948966
            >>> angle_between((1, 0, 0), (1, 0, 0))
            0.0
            >>> angle_between((1, 0, 0), (-1, 0, 0))
            3.141592653589793
    """
    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    return np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))

def align(obj, rot_setup):
    '''assumes that the reference and the mobule 
    points are the same object from different timesteps
    '''
    rot_matrix, point_source, point_destination = rot_setup
    return np.array([np.dot(rot_matrix,(p-point_source))+point_destination for p in obj])

def rotationSetup(ref_pts, mobile_pts):
    ref_vec = ref_pts[1] - ref_pts[0]
    mobile_vec = mobile_pts[1] - mobile_pts[0]

    # Calculating rotation matrix
    n1 = normal(mobile_vec, ref_vec)
    rot_matrix = rotation_matrix(n1, angle_between(mobile_vec, ref_vec))
    return rot_matrix, mobile_pts[0], ref_pts[0]

