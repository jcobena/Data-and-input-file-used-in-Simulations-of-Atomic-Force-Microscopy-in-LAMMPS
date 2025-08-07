'''======================  Functions for CONNECT section  ===================================== '''

from d_translations_and_rotations import *

from miscellaneous import *
#from b_hydrocarbon_chain import cc_recursivef_x, cc_recursivef_z
#from f_functional_termination_CH3 import create_last_H_CH3_x, create_last_H_CH3_y

'''The enumeration is as follows: First C, then H, then atoms in functional group'''


def cconections(lenchain):
    #For Carbons in chain (except if they belong to the terminal function group)
    cconected = []
    conect = 'CONECT'
    i=1
    
    #carbon connected
    while i != lenchain:
        conection1 = [conect, i, i+1]
        cconected.append(conection1)
        i = i + 1
    #H attached to C, notice that the number of H are twice as C (except functional group)
    j = 1
    k = lenchain+1
    while j != lenchain+1 :
        conection2 = [conect, j, k, k+1]
        conection3 = [conect, k, j]
        conection4 = [conect, k+1, j]
        conection5 = [conect, j, k]
        conection6 = [conect, j, k+1]
        cconected.append(conection2)
        cconected.append(conection3)
        cconected.append(conection4)
        cconected.append(conection5)
        cconected.append(conection6)
        j = j + 1
        k = k + 2
        
    for x in cconected:
        print(x)
        
    return cconected
    
'''PDB file lines'''
def pdblines(cchain, hydrogens, fgroup):
    pdblineslist = []
    atoms = 'ATOM'
    xxx = 'XXX'
    anumb = 'MOL'
    one = '1.00'
    zero = '0.00'
    
    structurecoords = cchain + hydrogens + fgroup
    
    lchain = len(cchain)
    lhydro = len(hydrogens)
    lfgroup = len(fgroup)
    
    n = 1
    for x in structurecoords:
        print(x)
        final_line = [atoms, n, anumb, xxx, n, format(x[0], '.3f'),\
                    round(x[1],3), round(x[2],3),  one,  zero,  x[3]]
        n = n+1
        pdblineslist.append(final_line)
        
    return pdblineslist
    
   
    
        



