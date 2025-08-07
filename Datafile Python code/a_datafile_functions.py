import collections 
from collections import Counter
import numpy as np
import copy
import operator 


''' Functions for final datafile, no matter what functional group is there. Impropers are neglected though'''

def generate_datafile(atoms_list, bonds_list, angles_list, dihedrals_list, dimensions, mass_list):
        
    natoms  = atoms_list[-1][0]
    nbonds  = bonds_list[-1][0]
    nangles = angles_list[-1][0]
    ndihedr = dihedrals_list[-1][0]
    nimprop = 0
    
    nattypes = sorted(atoms_list, key=operator.itemgetter(2))[-1][2]
    nbondtypes = sorted(bonds_list, key=operator.itemgetter(1))[-1][1]
    nangletypes = sorted(angles_list, key=operator.itemgetter(1))[-1][1]
    ndihedtypes = sorted(dihedrals_list, key=operator.itemgetter(1))[-1][1]
    #ndihedtypes = 'something'
    nimproptypes = 0
    
    atoms = 'atoms'
    angles = 'angles' 
    bonds = 'bonds'
    dihedrals = 'dihedrals'
    impropers = 'impropers'
    
    atom_types = 'atom types'
    bond_types = 'bond types'
    angle_types = 'angle types'
    dihedral_types = 'dihedral types'
    improper_types = 'improper types'
        
    the_header = [
    ['This is a LAMMPS datafile \n'],
    [natoms , atoms],
    [nbonds , bonds],
    [nangles , angles],
    [ndihedr , dihedrals],
    [nimprop , impropers, "\n"],
    
    [nattypes, atom_types],
    [nbondtypes, bond_types],
    [nangletypes, angle_types],
    [ndihedtypes, dihedral_types, "\n"]]
    #[nimproptypes, improper_types, ],
    
    body_data_file = [the_header, dimensions,' ', [['Masses \n']], mass_list,' ', [['Atoms \n']], atoms_list, ' ', [['Bonds \n']], 
    bonds_list, ' ', [['Angles \n']], angles_list, ' ', [['Dihedrals \n']], dihedrals_list] 
    
    return body_data_file

