import collections as col
'''creates a new pdb models where chosen residues are mutated to glycines 
'''

ResidueId = col.namedtuple('ResidueId', 'chain, id, iCode, name')
AtomId = col.namedtuple('AtomId', 'name, element, id')

class Residue(object):
    def __init__(self):
        self.atoms = []
        self.atoms_dict = {}
        self.res_id = None

    def set_res_id(self, res_id):
        if self.res_id:
            raise AttributeError('Residue Id is already set')
        self.res_id = res_id
    
    def add(self, atom):
        if not self.res_id:
            self.set_res_id(atom.res_id)
        elif self.res_id != atom.res_id:
            raise ValueError('Atom of %s does not belong to Residue with %s' % (self.res_id, atom.res_id))
        self.atoms.append(atom)
        self.atoms_dict[atom.type] = atom

    def to_glycine(self):
        self.res_id = ResidueId(self.res_id.chain, self.res_id.id, self.res_id.iCode, 'GLY')
        old_atom_list = self.atoms[:]
        self.atoms = []
        self.atoms_dict = {}
        for atom in old_atom_list:
            if atom.id.name in ['N', 'CA', 'C', 'O']:
                atom.res_id = self.res_id
                atom.update()
                self.add(atom)
    
    def __repr__(self):
        out = []
        for atom in self.atoms:
            out.append(str(atom))
        return '\n'.join(out)
        
class Atom(object):
    
    def __init__(self, line):
        self.line = line
        self.type = line[:6].strip()
        self.res_id = self.extract_res_id(self.line)
        self.name = line[12:16].strip()
        self.id = AtomId(line[12:16].strip(), line[76:78].strip(), int(line[6:11]))
    
    def __repr__(self):
        return self.line
    
    def update(self):
        self.line = '{:6}{:>5}  {:3} {:3} {}{:>4}{}   {}{:>2}{}'.format(
            self.type, self.id.id, self.id.name, 
            self.res_id.name, self.res_id.chain, self.res_id.id,
            self.res_id.iCode, self.line[30:76], self.id.element, self.line[78:])
    
    @staticmethod
    def extract_res_id(line):
        chain = line[21]
        id_ = int(line[22:26])
        iCode = line[26]
        name = line[17:20]
        return ResidueId(chain, id_, iCode, name)

class Protein(object):
    def __init__(self):
        self.residues = col.defaultdict(Residue)
        
    def __getitem__(self, key):
        return self.residues[key]

    def __setitem__(self, key, value):
        self.residues[key] = value

    def add(self, atom):
        self.residues[atom.res_id].add(atom)      
    
    def reindex(self):
        i = 0
        for res_id in sorted(self.residues):
            for atom in self.residues[res_id].atoms:
                i += 1
                if atom.id.id == i:
                    continue
                atom.id = AtomId(atom.id.name, atom.id.element, i)
                atom.update()
    
    def __repr__(self):
        out = []
        for res_id in sorted(self.residues):
            out.append(str(self.residues[res_id]))
        return '\n'.join(out)
    
class Glycinator(object):
    '''Class for creating mutated
    receptors, by changing residues
    to a Glycine.
    '''
    def __init__(self, receptor_path):
        self.receptor_path = receptor_path
        self.protein = Protein()
        self.read()

    def read(self):
        with open(self.receptor_path) as f:
            for line in f.read().split('\n'):
                if line[:6].strip() in ['ATOM', 'HETATM']:
                    atom = Atom(line)
                    self.protein.add(atom)
    
    def glycinate(self, res_ids):
        if isinstance(res_ids, tuple):
            resi_ids = [res_ids]
        elif isinstance(res_ids, list):
            pass
        else:
            raise ValueError('Intput must be tuple or list')
        for res_id in res_ids:
            self.protein.residues[res_id].to_glycine()
        self.protein.reindex()
    
    def write(self, path):
        with open(path, 'w') as f:
            f.write(repr(self.protein))
