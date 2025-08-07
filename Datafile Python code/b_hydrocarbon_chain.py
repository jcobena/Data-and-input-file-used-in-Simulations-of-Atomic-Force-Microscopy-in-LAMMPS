'''Create a hydrocarbon chain to be used as part of SAMs structures '''


'''Hydrocarbon chain functions '''

import numpy as np
from miscellaneous import *
#from h_global_var import *


def get_hydrocarbon_xyz(length_chain, bondl_sc, bondl_cc, angle_sc, angle_cc, xz_even, xz_even_3):
    
    mylist_x = []
    mylist_y = []
    mylist_z = []
    n = 1

    x_even = xz_even[0]
    z_even = xz_even[1]
    
    x_even_3 = xz_even_3[0]
    z_even_3 = xz_even_3[1]
 
    while n <= length_chain:
        
        my_x = cc_recursivef_x(n, bondl_sc, bondl_cc, angle_sc, angle_cc, x_even, x_even_3)
        my_y = 0.0
        my_z = cc_recursivef_z(n, bondl_sc, bondl_cc, angle_sc, angle_cc, z_even, z_even_3)
        
        mylist_x.append(my_x)
        mylist_y.append(my_y)
        mylist_z.append(my_z)
        
        n = n+1    
    xyzcoords = cc_get_xyz_coordinates(mylist_x, mylist_y, mylist_z, length_chain)
    return xyzcoords

def cc_get_xyz_coordinates(xlist, ylist, zlist, length_chain):
    
    mergedcoords = []
    n = 1
    for x,y,z in zip(xlist, ylist, zlist):
        if n == 1:
            atom_label = 'S'
        elif n == 2:
            atom_label = 'CH2'
        elif n == length_chain-1:
            atom_label = 'CH2ne'
        elif n == length_chain:
            atom_label = 'CH3'
        coords = [x,y,z, atom_label]
        mergedcoords.append(coords)
        n = n + 1
    return mergedcoords

def cc_recursivef_x(n, bondl_sc, bondl_cc, angle_sc, angle_cc, x_even, x_even_3):
    coords = []
    #initial conditions----
    if n == 1:
        x_n = 0.0
        coords.append(0.0)
        return x_n
    
    if n == 2:
        x_n = bondl_sc
        coords.append(x_n)
        return x_n 
    
    if n == 3:
        x_n = cc_recursivef_x(n-1, bondl_sc, bondl_cc, angle_sc, angle_cc, x_even, x_even_3) + x_even_3
        coords.append(x_n)
        return x_n 
        
    else:
        if n%2 == 0: #here, since n starts at 1 instead of 0, instead of using the ex_even I use the bond_length
            x_n = cc_recursivef_x(n-1, bondl_sc, bondl_cc, angle_sc, angle_cc, x_even, x_even_3) + bondl_cc
            coords.append(x_n)
        
        else:
            x_n = cc_recursivef_x(n-1, bondl_sc, bondl_cc, angle_sc, angle_cc, x_even, x_even_3) + x_even
            coords.append(x_n)

        return x_n 

def cc_recursivef_z(n, bondl_sc, bondl_cc, angle_sc, angle_cc, z_even, z_even_3):
    coords = []
    #initial conditions----
    if n == 1 or n == 2:
        z_n = 0.0
        coords.append(z_n)
        return z_n
    
    if n == 3:
        z_n = z_even_3
        coords.append(z_n)
        #print(y_n, bondl_cc)
        return z_n
    #-------------------------
    else:
        if n%2 != 0: #means it's an odd number
            myn = n + 1
            mylevel = myn/2
            z_n = (mylevel-1)*z_even + z_even_3 - z_even
            coords.append(z_n)
            return z_n
        
        else:
            z_n = cc_recursivef_z(n-1, bondl_sc, bondl_cc, angle_sc, angle_cc, z_even, z_even_3)
            coords.append(z_n)
            return z_n
    
# def cc_get_xyz_coordinates_molecules(nmol, xlist, ylist, zlist, length_chain):
#     
#     s_type = 1
#     ch2_type = 2
#     ch3_type = 3
#     
#     mergedcoords = []
#     id = 1
#     for x,y,z in zip(xlist, ylist, zlist):
#         if n == 1:
#             atype = s_type
#             q = 0
#         elif n == 2:
#             atype = ch2_type
#             q = 0
#         elif n == length_chain:
#             atype = ch3_type
#             q = 0
#         coords = [id, nmol, atype, q, x, y, z]
#         mergedcoords.append(coords)
#         id = id + 1
#     return mergedcoords

