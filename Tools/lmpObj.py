import numpy as np
from lmp_helpers import *

class LmpObj(object):
    def __init__(self, environment):
        self.env = environment
        self.mol_type = 'Unknown'
        self.Id = self.env.new_id['molecule']
        self.data  = {'particles': [],
                      'bonds':     [],
                      'angles':    [],
                      'dihedrals': []}
        self.box = np.zeros([3,2])

    def __repr__(self):
        return 'Molecule Type: %s %d | Length: %d' % (self.mol_type, self.Id, len(self.data['particles']))

    def name():
        doc = "The name property."
        def fget(self):
            return self.mol_type+' '+str(self.Id)
        def fset(self, value):
            print 'Not possible.'
        return locals()
    name = property(**name())

    def update(self):
        pass

    def calc_box(self, add_margin=False):
        margin = 0
        if add_margin:
            margin = float(self.globals['box_margin'])
        total_coords = np.zeros([len(self.data['particles']),3])
        box = np.zeros([3,2])
        for i, atom in enumerate(self.data['particles']):
            total_coords[i] = atom.position
        box[0,0] = total_coords.T[0].min()-margin
        box[0,1] = total_coords.T[0].max()+margin
        box[1,0] = total_coords.T[1].min()-margin
        box[1,1] = total_coords.T[1].max()+margin
        box[2,0] = total_coords.T[2].min()-margin
        box[2,1] = total_coords.T[2].max()+margin
        self.box = box

    def mol_info(self):
        return self.mol_type + str(self.Id)

    def move(self, move_coord):
        '''Moves all particles by the vector specified in move_coord.
        move_coord must be in Cartesian coordinates.

        input:
            move_coord: np.array([]) shape: (3)
        '''
        for particle in self.data['particles']:
            particle.position += move_coord
        self.calc_box()

    def rotate(self, axis, angle):
        '''Rotates all particles by the degree of p:angle around
        the p:axis.
        The pivot point of the object is at the geometrical center.

        input:
            axis: np.array([]); shape: (3)
            angle: float;       [radiant]
        '''
        rot_matrix = rotation_matrix(axis, angle)

        # selecting minor box corner as pivot
        pivot = (self.box[:,0]+self.box[:,1])/2.0
        # calculate vector from pivot to particle in spherical coordinates
        vectors = np.array([p.position for p in self.data['particles']]) - pivot
        rotated_vectors = np.array([np.dot(rot_matrix, v) for v in vectors])
        end_coords = pivot + rotated_vectors
        for i, particle in enumerate(self.data['particles']):
            particle.position = end_coords[i]
        self.calc_box()

    def center(self):
        return (self.box[:,0] + self.box[:,1])/2.0
