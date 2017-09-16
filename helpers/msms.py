import pathlib2 as pl
import os
import subprocess as sp
import tempfile
import numpy as np
import pandas as pd
import collections as col
import Bio.PDB as PDB
import PolyLibScan.Database.db as DB

class HydrophobicParameterisation(object):
    
    def __init__(self, pdb_path, 
                 pdb_to_xyzrn=None, atmtypenumbers=None, msms_bin=None, verbose=False):
        self.own_path = pl.Path(__file__).parent
        self.path = self.set_paths(pdb_to_xyzrn, atmtypenumbers, msms_bin)
        self.pdb_path = pdb_path
        self.residues = {}
        self.phob_resis = {}
        self.phil_resis = {}
        self.correlation_parameter = 0.02
        self.parser = PDB.PDBParser()
        self.struc = self.parser.get_structure('protein', self.pdb_path.as_posix())
        
    @property
    def pdb_path(self):
        return self._pdb_path

    @pdb_path.setter
    def pdb_path(self, value):
        '''Setter makes sure that the path is a PosixPath object
        '''
        if isinstance(value, basestring):
            self._pdb_path = pl2.Path(value)
        else:
            self._pdb_path = value
        self.check_pdb_path()

    def check_pdb_path(self):
        if not self.pdb_path.exists():
            raise IOError('pdb file could not be found.')

    def set_paths(self, pdb_to_xyzrn, atmtypenumbers, msms_bin):
        paths = {}
        if atmtypenumbers and pl.Path(atmtypenumbers).exists():
            paths['cmd_dir'] = pl.Path(atmtypenumbers).parent.absolute().resolve()
        else:
            paths['cmd_dir'] = self.own_path.joinpath('../external_bins/msms').absolute().resolve()
            
        if pdb_to_xyzrn and pl.Path(pdb_to_xyzrn).exists():
            paths['pdb_to_xyzrn'] = pl.Path(pdb_to_xyzrn)
        else:
            paths['pdb_to_xyzrn'] = self.own_path.joinpath('../external_bins/msms/pdb_to_xyzrn').absolute().resolve()
        
        if msms_bin and pl.Path(msms_bin).exists():
            paths['msms_bin'] = pl.Path(msms_bin)
        else:
            paths['msms_bin'] = self.own_path.joinpath('../external_bins/msms/msms.x86_64Linux2.2.6.1.staticgcc').absolute().resolve()
        return paths
    
    def msms_surface(self, pdb_file=None, verbose=False):
        if pdb_file:
            path = pdb_file
        else:
            path = self.pdb_path

        input_file = pl.Path(path).absolute().resolve().as_posix()
        msms_folder = tempfile.mkdtemp()
        xyz_file = tempfile.mktemp(dir=msms_folder)
        with open(xyz_file, 'w') as f:
            xyz_run = sp.Popen('%s %s' % (self.path['pdb_to_xyzrn'].as_posix(), input_file), 
                               cwd=self.path['cmd_dir'].as_posix(),
                               shell=True, stdin=None, stdout=f, stderr=sp.PIPE)
            xyz_out = xyz_run.communicate()
        if xyz_out[1] != '':
            raise RuntimeError('pdb_to_xyz:\n %s' % xyz_out[1])
            
        file_path = os.path.join(msms_folder,'msms')

        msms = sp.Popen('%s  -probe_radius 1.5 -if %s -of %s -af %s' % (self.path['msms_bin'], 
                                                                        xyz_file, file_path, file_path),
                        shell=True, stdin=None, stdout=sp.PIPE, stderr=sp.PIPE)
        msms_out = msms.communicate()
        if verbose:
            print out[0]
        if msms_out[1] != '':
            raise RuntimeError('msms:\n %s' % xyz_out[1])
        vertices_file = pl.Path(file_path + '.vert')
        if not vertices_file.exists():
            raise IOError('%s does not exist.' % vertices_file)
        face_file = pl.Path(file_path + '.face')
        if not face_file.exists():
            raise IOError('%s does not exist.' % face_file)
        area_file = pl.Path(file_path + '.area')
        if not area_file.exists():
            raise IOError('%s does not exist.' % area_file)
        self.path['msms_vertex'] = vertices_file
        self.path['msms_face'] = face_file
        self.path['msms_area'] = area_file
        return vertices_file, face_file, area_file
    
    def area(self, path=None):
        '''Parse msms-area file and return a numpy array.
        '''
        if path:
            area_file = pl.Path(path)
        else:
            area_file = self.path['msms_area']
        
        with open(area_file.as_posix()) as f:
            content = f.read().split('\n')[1:-1]
        atom_count = len(content)
        data = np.zeros(atom_count, dtype=[('count', np.int), ('atom_id', 'S4'), ('resn', 'S3'), 
                                           ('chain', 'S1'), ('resi', np.int), ('iCode', 'S1'),
                                           ('ses', np.float), ('sas', np.float)])
        for i, line in enumerate(content):
            cnt, ses, sas, id_string = line.split()
            a_id, resn, chain, resi = id_string.split('_')
            try:
                res_id = int(resi)
                icode = ' '
            except ValueError:
                res_id = int(resi[:-1])
                icode = resi[-1]
            data[i] = (int(cnt), a_id, resn, chain, res_id, icode, float(ses), float(sas))
        return data
    
    def faces(self, path=None):
        if path:
            face_file = pl.Path(path)
        else:
            face_file = self.path['msms_face']
        with open(face_file.as_posix()) as f:
            # discard first two lines
            f.next()
            f.next()
            n_faces, n_spheres, dens, probe_r = f.next().split()
            data = np.zeros(int(n_faces), dtype=[('vertices', np.int,3), ('type', np.int8), 
                                                 ('f_no', np.int)])
            for i,face in enumerate(f):
                face_data = map(int, face.strip().split())
                data[i] = (tuple(face_data[:3]), face_data[3], face_data[4])
        return data
    
    def vertices(self, path=None):
        if path:
            vertex_file = pl.Path(path)
        else:
            vertex_file = self.path['msms_vertex']
        with open(vertex_file.as_posix()) as f:
            # discard first two lines
            f.next()
            f.next()
            n_vertex, n_spheres, dens, probe_r = f.next().split()
            data = np.zeros(int(n_vertex), dtype=[('xyz', np.float, 3), ('normal', np.float, 3), 
                                                 ('face_id', np.int), ('sph_idx', np.int), ('type', np.int8),
                                                 ('atom', 'S4'), ('resn', 'S3'), ('chain', 'S1'), ('resi', np.int), 
                                                 ('icode', 'S1')])
            for i,line in enumerate(f):
                vertex = line.split()
                try:
                    coords = tuple(map(float, vertex[:3]))
                except:
                    print vertex[:3]
                    print i
                    break
                normal = tuple(map(float, vertex[3:6]))
                info = tuple(map(int, vertex[6:9]))
                res_str = vertex[9].split('_')
                try:
                    res_id = int(res_str[3])
                    icode = ' '
                except ValueError:
                    res_id = int(res_str[3][:-1])
                    icode = res_str[3][-1]
                res_info = res_str[0], res_str[1], res_str[2], res_id, icode
                data[i] = (coords, normal) + info + res_info
        return data
    
    def get_residue_hydro_levels(self, parameters_file=None):
        '''Read in the clogP values of all residues from the 
        parameter file and sort all members into sets hydrophobic 
        and hydrophilic particles.
        '''
        if parameters_file and pl.Path(parameters_file).exists():
            path = pl.Path(parameters_file)
        else:
            path = self.own_path.joinpath('../parameters/clogP.h5')
        store = pd.HDFStore(path.as_posix())
        hy = store['clogP']
        store.close()
        hydrophobic = hy[hy>0]
        hydrophilic = hy[hy<0]

        self.hydrophobic = set(hydrophobic.index)
        self.hydrophilic = set(hydrophilic.index)

    def create_residues(self, parent, surface_resis):
        for res_id in surface_resis:
            self.residues[res_id] = Residue(parent, res_id)
            if res_id[0] in self.hydrophilic:
                self.phil_resis[res_id] = self.residues[res_id]
            else:
                self.phob_resis[res_id] = self.residues[res_id]

    def add_area(self, areas):
        '''add area information to residue and 
        delete the residue, if it has an area of 0.0.
        '''
        for res_id, area in areas.items():
            if area > 0.0:
                self.residues[res_id].area = area
            else:
                del self.residues[res_id]
                if res_id in self.phob_resis:
                    del self.phob_resis[res_id]
                else:
                    del self.phil_resis[res_id]

    def add_vertices(self, vertices):
        for vert in vertices:
            res_id = (vert[-4], vert[-3], vert[-2], vert[-1])
            # there are corner cases, where there is no area 
            # but vertices. In case the residue has no area,
            # we discard the vertices.
            if res_id in self.residues:
                self.residues[res_id].vertices.append(vert[0]) 

    def resi_surface(self, area_data):
        res_area = col.defaultdict(float)
        for atom in area_data:
            res_id = (atom[2], atom[3], atom[4], atom[5])
            res_area[res_id] += atom[7]
        return res_area

    def calculate_distances(self):
        for residue in self.phob_resis.values():
            residue.distance_to(self.phil_resis)

    def add_neighbors(self):
        for residue in self.phob_resis.values():
            residue.neighboring_resis(self.phob_resis)

    def additive_area(self):
        for residue in self.phob_resis.values():
            residue.additive_area()

    def complete(self):
        self.get_residue_hydro_levels()
        self.msms_surface()
        area = self.resi_surface(self.area())
        self.create_residues(self, area.keys())
        self.add_area(area)
        vertices = self.vertices()
        self.add_vertices(vertices)
        self.calculate_distances()
        self.add_neighbors()
        self.additive_area()

    def to_numpy_array(self, threshold=0.0):
        arr_type = [('resname', 'S3'), ('chain', 'S1'), ('id', np.int16), 
                    ('iCode', 'S1'), ('singleParameter', np.float16), ('areaParameter', np.float16)]
        hydrophobic_array = np.zeros(len(self.phob_resis), dtype=arr_type)
        for i,resi in enumerate(self.phob_resis.values()):
            hydrophobic_array[i] = (resi.name[0], resi.pdb_id[0], resi.pdb_id[1], resi.pdb_id[2], 
                                    resi.hydrophobic_energy_single(), resi.hydrophobic_energy_area())
        mask = hydrophobic_array['singleParameter'] > threshold
        return hydrophobic_array[mask]

    def to_hdf5(self, db_path, table_name, threshold=0.1):
        '''Save residue info and parameters to HDF5 Database.
        '''
        hydrophobic_array = self.to_numpy_array(threshold)
        dBase = DB.Database(db_path)
        dBase._save_table(hydrophobic_array, '/', table_name)
        dBase.close()

    def to_dataframe(self, threshold=0.0):
        hydrophobic_array = self.to_numpy_array(threshold)
        df = pd.Dataframe(data=hydrophobic_array)
        return df

class Residue(object):

    def __init__(self, parent, name):
        self.name = name
        self.parent = parent
        self.clogP = None
        self.min_radius = -1
        self.neighbors = set([])
        self.area = -1
        self.surrounding_area = -1
        self.vertices = []
        self.distance = -1
        self.pdb_id = None
        self.set_pdb_id()

    def center(self):
        return np.array(self.vertices).mean(axis=0)

    def distance_to(self, group):
        self.distance = self.closest_of(group)
        return self.distance

    def closest_of(self, group):
        min_dist = 1000.0
        res_center = self.center()
        for residue in group.values():
            for vertex in residue.vertices:
                dist = np.linalg.norm(res_center - vertex)
                if dist < min_dist:
                    min_dist = dist
        return min_dist

    def neighboring_resis(self, hphob_resis):
        if self.distance < 0:
            raise ValueError('calculate distance first!')
        center = self.center()
        for residue_id, hphob_resi in filter(lambda x:x[0]!=self.name, hphob_resis.items()):
            for vert in hphob_resi.vertices:
                if np.linalg.norm(vert-center) < self.distance:
                    self.neighbors.add(hphob_resi)
                    break

    def additive_area(self):
        if self.area < 0:
            raise ValueError('Set area first.')
        self.surrounding_area = self.area
        for neighbor in self.neighbors:
            self.surrounding_area += neighbor.area

    def hydrophobic_energy_single(self):
        '''
        See: Reynolds1974
        '''
        return self.parent.correlation_parameter * self.area 

    def hydrophobic_energy_area(self):
        '''
        See: Reynolds1974
        '''
        return self.parent.correlation_parameter * self.surrounding_area

    def set_pdb_id(self):
        matches = []
        for pdb_resi in self.parent.struc[0].get_residues():
            pdb_id = self.reduced_id(pdb_resi.get_full_id())
            if self.name[1:] == pdb_id:
                matches.append(self.reduced_id(pdb_resi.get_full_id()))
        if len(matches) !=1:
            raise ValueError('Amino Acid %s was not matched to a single pdb residue but %s' % (self.name, matches))
        self.pdb_id = matches[0]

    @staticmethod
    def reduced_id(full_id):
        return (full_id[2], full_id[3][1], full_id[3][2])

