#!/usr/bin/env python
import networkx as nx
import numpy as np
import os
from copy import deepcopy

"""
  Simple molecule class, holds atom types, positions and bonds.
"""
class molecule:
    def __init__(self):
        '''
        All atom information is stored in a fixed order:
        atom-id(0),  belong-to-chain, coordinate, bonded-atoms, etc. 
        atom-id(1),  belong-to-chain, coordinate, bonded-atoms, etc. 
        ...
        '''
        self.n_atom      = 0   # Number of beads in system.
        self.atom_coords = []  # Contains the wrapped coordinate of the current atom.
        self.crystal_atom = [] # Contains the belonging-phase of the current atom.
                               # [0 ~ amorphous, 1 ~ crystalline]
        self.n_chain     = 0   # Number of chains in system.
        self.chainids    = []  # Contains the ids of chains in system. 0-indexed
        self.chain       = []  # Contains which chain the current atom belongs to.
        self.chains      = []  # Contains nested lists of atoms belonging to the same chains.
        self.sorted_chains = []  # Nexted list, clusters the bonded atoms into sublist.
                                 # The atom-ids are sorted according to connectivity.
        self.ends        = []  # Contains the ids of all end beads in cg system.
                               # [chain-id][atom-id].
        self.bonds       = []  # Contains which atoms the current atom is connected with.
                               # atom-id is 0-indexed.
        self.box         = []  # Contains the bounds of box in each direction.
        self.box_length  = []  # Contains the lenght of box in each direction.
        self.time        = 0   # The simulation time.

    # Rerurns the end-atoms of a chain with chain-id c.
    def chain_ends(self, c):
        atoms = [i for i, x in enumerate(self.chain) if x == c]
        # Determine end atoms in chain c.
        ends  = []
        for b in atoms:
            if len(self.bonds[b]) == 1:
                ends.append(b)
        return ends

    # Rerurns all the atoms belonging to a chain with chain-id c, sorted according to connectivity.
    def atoms_in_chain(self, c):
        # Determine end atoms in chain c.
        ends = self.chain_ends(c)
        # Sort atoms according to connectivity.
        atoms = [ends[0]]
        queue = self.bonds[ends[0]][:]  # Use deepcopy to avoid modifying self.bonds
        while len(queue) > 0:
            next = queue[-1]
            queue.pop(-1)               # We modify queue here
            atoms.append(next)
            for neighbor in self.bonds[next]:
                if neighbor not in atoms:
                    queue.append(neighbor)
        return atoms


"""
  Read coordinate and bond information from LAMMPS data file.
"""
def read_system(lmp_data):
    # Initialize molecule class.
    mol = molecule()
    fid   = open(lmp_data,'r') 
    # Read # of atoms in system 
    while True:
        line = fid.readline()
        if len(line.split()) == 4:
            break
        if 'atoms' in line:
            mol.n_atom = int(line.split()[0])
    # Read box scopes.
    line = line.split()
    mol.box.append([ float(j) for j in line[0:2] ])
    for i in range(2):
        line = fid.readline().split()
        mol.box.append([ float(j) for j in line[0:2] ])
    for i in range(3):
        mol.box_length.append( abs(mol.box[i][1]-mol.box[i][0]) )

    # Read bead coordinates.
    # [atom-id,mol-id,atom-type,x,y,z,image-x,image-y,image-z]
    while True:
        line = fid.readline()
        if line.startswith('Atoms'):
            fid.readline()
            break
    mol.atom_coords = np.zeros((mol.n_atom, 3))
    mol.chain       = [-1]*mol.n_atom
    # Read coordinates.
    atomids = range(mol.n_atom)
    for i in atomids:
        ll = fid.readline().split()
        atom_id = int(ll[0])-1  # Convert from 1-indexed to zero-indexed.
        mol.chain[atom_id] = int(ll[1])
        x = [float(s) for s in ll[3:6]]
        x = put_in_box(x, mol.box, mol.box_length)
        mol.atom_coords[atom_id] = x
    mol.chainids = list(set(mol.chain))
    mol.chainids.sort()
    # Read bonds.
    # [bond-id,bond-type,atom-id1,atom-id2]
    # For now mol.bonds contains the bonding list, not yet those atoms the current atom is connected with.
    while True:
        ll = fid.readline()
        if ll.startswith('Bonds'):
            fid.readline()
            break
    while True:
        ll = fid.readline().split()
        if len(ll) < 4 or not ll: break
        # Convert connectivity from 1-indexed to 0-indexed.
        mol.bonds.append( [int(s)-1 for s in ll[2:4]] )

    # Sort mol.bonds according to index.
    mol.bonds.sort(key=lambda x: x[0])
    mol.atom_coords = np.array(mol.atom_coords)

    # Build graph to generate the neighbor list for all beads.
    # https://networkx.github.io/documentation/latest/reference/classes.graph.html
    G = nx.Graph()
    nodes = atomids
    edges = [tuple(b) for b in mol.bonds]
    G.add_nodes_from(nodes)
    G.add_edges_from(edges)
    mol.bonds = [ list(G.neighbors(i)) for i in nodes ]
    # Determine number of chains in system.
    mol.chains  = list(nx.connected_components(G))
    mol.n_chain = len(mol.chains)

    # Determine all end atoms in system.
    for c in mol.chainids:
        mol.ends.append(mol.chain_ends(c))
        mol.sorted_chains.append(mol.atoms_in_chain(c))
    return mol


def read_frame_positions(lmp_trj):
    '''
    Read positions in trajectory file corresponding to
    time-step and atom-data.
    '''
    ts_pos   = []
    data_pos = []
    fid = open(lmp_trj,'r')
    while True:
        line = fid.readline()
        if not line or len(line.split())==0:
            break
        if line.startswith('ITEM: TIMESTEP'):
            ts_pos.append(fid.tell())
        elif line.startswith('ITEM: ATOMS id'):
            data_pos.append(fid.tell())
    fid.close()
    return ts_pos, data_pos


def update_mol(scaled, mol_ref, lmp_trj, tsk, dsk, dt):
    '''
    Read coordinate from trajectory file at a specified time frame.
    If the coordinates in lmp_trj is unwrapped, then the updated coordinates are also unwrapped.
    '''
    mol = deepcopy(mol_ref)
    fid = open(lmp_trj,'r')
    # Identify the indices of 'x', 'y', 'z'.
    while True:
        line = fid.readline()
        if line.startswith('ITEM: ATOMS id'):
            line = line.split()
            xid = line.index('x') - 2
            yid = line.index('y') - 2
            zid = line.index('z') - 2
            break
    # Update simulation time.
    fid.seek(tsk)
    mol.time = float(fid.readline().split()[0]) * dt
    # Update box information.
    for i in range(3):
        fid.readline()
    for i in range(3):
        line = fid.readline().split()
        mol.box[i] = [float(j) for j in line[0:2]]
        mol.box_length[i] = abs(mol.box[i][1]-mol.box[i][0])
    # Modify box.
    midz       = np.average(mol.box[2])
    mol.box[2] = [midz, midz+mol.box_length[2]]
    # Update coordinate information.
    fid.seek(dsk)
    while True:        
        line = fid.readline()
        if not line or line.startswith('ITEM: TIMESTEP'):
            break
        atom_id = int(line.split()[0])-1
        xs = [float(n) for n in line.split()[xid:zid+1]]
        if scaled:
            [[xl,xh],[yl,yh],[zl,zh]] = mol.box
            x = xl + xs[0]*(xh-xl)
            y = yl + xs[1]*(yh-yl)
            z = zl + xs[2]*(zh-zl)
            xs = [x,y,z]
        xs = put_in_box(xs, mol.box, mol.box_length)
        mol.atom_coords[atom_id] = xs
    return mol


''' Put atom x in box using periodic boundary conditions. '''
def put_in_box(x, box, box_length):
    for i in range(3):
        while x[i] < min(box[i]): x[i] += box_length[i]
        while x[i] > max(box[i]): x[i] -= box_length[i]
    return x


'''
  Take an atom x0 and another atom x in box,
  modify the coordinate of x according to periodic boundary conditions,
  so that x is in the vicinity of x0.
'''
def move_to_vicinity(x, x0, box_length):
    dl = x - x0
    for i in range(3):
        while dl[i] < -0.5*box_length[i]: dl[i] += box_length[i]
        while dl[i] > 0.5*box_length[i]:  dl[i] -= box_length[i]
    return x0 + dl


''' Write lammps trajectory file for semicrystalline system. '''
def write_trj_shifted(sys, mol, target):
    wid = open(target,'a')
    wid.write('ITEM: TIMESTEP\n')
    wid.write('{}\n'.format(mol.time))
    wid.write('ITEM: NUMBER OF ATOMS\n')
    wid.write('{}\n'.format(mol.n_atom))
    wid.write('ITEM: BOX BOUNDS pp pp pp\n')
    for i in range(3):
        [lo,hi] = mol.box[i]
        wid.write('{} {}\n'.format(lo,hi))
    wid.write('ITEM: ATOMS id mol type x y z\n')
    for i in range(mol.n_atom):
        achain  = mol.chain[i]
        atype   = mol.type[i]
        [x,y,z] = list(mol.atom_coords[i])
        row = (i+1, achain, atype, x, y, z)
        wid.write(' %d %d %d %.4f %.4f %.4f\n' %row)

