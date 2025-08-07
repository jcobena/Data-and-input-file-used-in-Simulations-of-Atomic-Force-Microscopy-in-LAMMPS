'''Create a hydrocarbon chain to be used as part of SAMs structures. 
    The hydrocarbon chain ends with a CN terminal group'''


#import fileinput
import numpy as np
import math as mt
import random
import glob
import os
import copy
import scipy
#from b_hydrocarbon_chain import *
from d_translations_and_rotations import *
from miscellaneous import *
from i_connections import *

#length of chain:
length_chain = 14

#bond length in \AA
blSC = 1.815
blCC = 1.54
blCNH2 = 1.454
blNh = 1.002

#hydrocarbon angle (degrees)
cc_angle = 109.5
sc_angle = 114.4
cnh_angle = 120.0
nhh_angle = 120.0
H_angle = 60 #careful with this one


#x and zeven for C in CN... let's see
#xz_even_CN   = x_even_z_even(cc_angle, blC11C12) #for the C in CN
xz_even_NH2   = x_even_z_even(cnh_angle, blCNH2) #for the C in CN
xz_even_HinNH2 = x_even_z_even(nhh_angle, blNh)

#for the last missing H:
xz_even_H = x_even_z_even(H_angle, blNh)


def driver_NH2():
    
    #for the first 3 atoms: S-C-C
    xz_even_n = x_even_z_even(cc_angle, blCC)
    xz_even_3 = x_even_z_even(sc_angle, blCC)
    
    #create C of hydrocarbon chain
    hydrocarbon_chain = get_hydrocarbon_xyz(length_chain, blSC, blCC, sc_angle, cc_angle, xz_even_n, xz_even_3)

    #rotate the hydrocabron chain:
    for i in hydrocarbon_chain:#[1:]:
        #print(i)
        j = rotate_counterclockwise_about_y_axis(i, 35)
        i[0] = j[0]
        i[1] = j[1]
        i[2] = j[2]
    #the rotation is done, now we sent this to the main program so that it can be translated
    return hydrocarbon_chain

            
def get_hydrocarbon_xyz(length_chain, bondl_sc, bondl_cc, angle_sc, angle_cc, xz_even, xz_even_3):
    
    mylist_x = []
    mylist_y = []
    mylist_z = []

    x_even = xz_even[0]
    z_even = xz_even[1]
    
    x_even_3 = xz_even_3[0]
    z_even_3 = xz_even_3[1]
    
    n = 1
    
    while n <= length_chain:
        
        my_x = cc_recursivef_x(n, bondl_sc, bondl_cc, angle_sc, angle_cc, x_even, x_even_3)
        my_y = 0.0
        my_z = cc_recursivef_z(n, bondl_sc, bondl_cc, angle_sc, angle_cc, z_even, z_even_3)
        
        mylist_x.append(my_x)
        mylist_y.append(my_y)
        mylist_z.append(my_z)
        
        n = n + 1    
        
    xyzcoords = cc_get_xyz_coordinates(mylist_x, mylist_y, mylist_z, length_chain)
    
    #=================================
    
    n_in_nh2 = xyzcoords[-2]
        
    #p1 is C11 and o2 is C in CN
    
    x_n_in_nh2 = n_in_nh2[0]
    y_n_in_nh2 = n_in_nh2[1]
    z_n_in_nh2 = n_in_nh2[2]
    
    Hx = x_n_in_nh2 + xz_even_H[0] 
    Hy = y_n_in_nh2
    Hz = z_n_in_nh2 + xz_even_H[1]
    atom_label = 'HinNH2'
    Ncoords = [Hx, Hy, Hz, atom_label]
        
    xyzcoords.append(Ncoords)
    
    #=================================
    
    #add the last H:
    
    

    return xyzcoords
        

def cc_get_xyz_coordinates(xlist, ylist, zlist, length_chain):
    
    mergedcoords = []
    n = 1
    for x,y,z in zip(xlist, ylist, zlist):
        if n == 1:
            atom_label = 'S'
        # elif n == 2:
        #     atom_label = 'CH2'
        elif n == length_chain-3:
            atom_label = 'C10'
        elif n == length_chain-2:
            atom_label = 'C11'
        elif n == length_chain-1:
            atom_label = 'NinNH2'
        elif n == length_chain:
            atom_label = 'HinNH2'
        else:
            atom_label = 'CH2'
            
            
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
            x_n = cc_recursivef_x(n-1, bondl_sc, bondl_cc, angle_sc, angle_cc, x_even, x_even_3) + bondl_cc #x_even
            coords.append(x_n)
            
            if n == 14:
                x_n = cc_recursivef_x(n-1, bondl_sc, bondl_cc, angle_sc, angle_cc, x_even, x_even_3) + blNh #xz_even_HinNH2[0]
                coords.append(x_n)
            
        else:
            if n == 13:
                x_n = cc_recursivef_x(n-1, bondl_sc, bondl_cc, angle_sc, angle_cc, x_even, x_even_3) + xz_even_NH2[0]
                coords.append(x_n)
                
            else:
                x_n = cc_recursivef_x(n-1, bondl_sc, bondl_cc, angle_sc, angle_cc, x_even, x_even_3) + x_even #bondl_cc 
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
            
            if n == 13:
                myn = n + 1
                mylevel = myn/2
                z_n = (mylevel-3)*z_even + z_even_3  + xz_even_NH2[1]
                coords.append(z_n)
                return z_n
                          
            else:
                myn = n + 1
                mylevel = myn/2
                z_n = (mylevel-1)*z_even + z_even_3 - z_even
                coords.append(z_n)
                return z_n
        
        else:
            z_n = cc_recursivef_z(n-1, bondl_sc, bondl_cc, angle_sc, angle_cc, z_even, z_even_3)
            coords.append(z_n)
            return z_n
        
            if n == 14:
                z_n = cc_recursivef_z(n-1, bondl_sc, bondl_cc, angle_sc, angle_cc, z_even, z_even_3) + xz_even_HinNH2[1]
                coords.append(z_n)
                return z_n
    