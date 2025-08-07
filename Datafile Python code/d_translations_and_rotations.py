'''====================== Geometry functions: midpoints, translations and rotations ============================ '''



''' Preparation Functions to get D (diagonal) and L (length of the cube) '''

import math as mt
import numpy as np
from miscellaneous import *

def call_D_L_midpoint(f1, f2, f3, c0, c2):
    #This function calls get_d, get_L_of_cube and get_midpoint
    D = f1(c0, c2)
    L = f2(D)
    mdp = f3(c0, c2)
    
    return [L, mdp]


def get_D(c0, c2):
    #find long diagonal. The diagonal is formed by joining the 2 carbons that are attached to the even c.
    #In other words, distance between the 2 points (odd/even Cs)
    p0x = c0[0]
    p0y = c0[1]
    p0z = c0[2]
    
    #odd  C 2 (although is carbon number 3):
    p2x = c2[0]
    p2y = c2[1]
    p2z = c2[2]
    
    D = mt.sqrt ((p2x - p0x)**2 + (p2y - p0y)**2 + (p2z - p0z)**2)
    
    return D


def get_midpoint(c0, c2):
    #get the midpoint of D or any line
    xm = (c0[0] + c2[0])/2
    ym = (c0[1] + c2[1])/2
    zm = (c0[2] + c2[2])/2
    #print(zm)
    midpoint = [xm, ym, zm]
        
    return midpoint

def get_L_of_cube(D):
    #length of the side of the cube:
    L = D/(mt.sqrt(2))
    
    return L


''' General translation functions '''

def translate_a_point_to_0(p1, mdpoint):
    #move the system such that the midpoint is 0.0, 0.0, 0.0.
    #assumes that the point is located in the +x, +y and +z area
    p1x = p1[0]
    p1y = p1[1]
    p1z = p1[2]
    #1st, odd carbon (1)
    p1x_t = p1x - mdpoint[0]
    p1y_t = p1y - mdpoint[1]
    p1z_t = p1z - mdpoint[2]
    #print("p1z_t =",p1z_t)
    
    return [p1x_t, p1y_t, p1z_t]

''' Translation functions for even C'''
def translate_a_point_L_and_beta_even(c0, L, hc_angle):
    
    angle_to_move = ((180-hc_angle)/2)*ang_to_rad
    
    p1x = c0[0]
    p1y = c0[1]
    p1z = c0[2]
    
    Lx = L*np.sin(angle_to_move)
    Ly = L*np.cos(angle_to_move)
    Lz = p1z #this one does not change
    
    p1x_t = p1x - Lx 
    p1y_t = p1y + Ly 
    p1z_t = Lz

    return [p1x_t, p1y_t, p1z_t]


''' Translation functions for odd C'''
def reverse_translation_from_0(p1, mdpoint):
    p1x = p1[0]
    p1y = p1[1]
    p1z = p1[2]
    #1st, odd carbon (1)
    p1x_t = p1x + mdpoint[0]
    p1y_t = p1y + mdpoint[1]
    p1z_t = p1z + mdpoint[2]

    return [p1x_t, p1y_t, p1z_t]



def translate_a_point_L_and_beta_odd(c0, L, hc_angle):
    
    angle_to_move = ((180-hc_angle)/2)*ang_to_rad
    
    p1x = c0[0]
    p1y = c0[1]
    p1z = c0[2]
    
    Lx = L*np.sin(angle_to_move)
    Ly = L*np.cos(angle_to_move)
    Lz = p1z #this one does not change
    
    p1x_t = p1x + Lx 
    p1y_t = p1y - Ly 
    p1z_t = Lz

    return [p1x_t, p1y_t, p1z_t]

''' Counterclockwise rotation functions'''
def rotate_counterclockwise_about_z_axis(mypoint, hc_angle):
    
    rad_angle = ((180 - hc_angle)/2)*ang_to_rad
    
    x = mypoint[0]
    y = mypoint[1]
    z = mypoint[2]
    
    #use rotation matrix for Rz(gamma)
    newx = (x*np.cos(rad_angle)) - (y*np.sin(rad_angle))
    newy = (x*np.sin(rad_angle)) + (y*np.cos(rad_angle))
    newz = z
    
    return [newx, newy, newz]

def rotate_counterclockwise_about_y_axis(mypoint, roty_degree):
    
    rad_angle = roty_degree*ang_to_rad
    #print(rad_angle)
    
    x = mypoint[0]
    y = mypoint[1]
    z = mypoint[2]
    
    newx = (x*np.cos(rad_angle)) - (z*np.sin(rad_angle))
    newy = y
    newz = (+x*np.sin(rad_angle)) + (z*np.cos(rad_angle))
    
    return [newx, newy, newz]

''' Clockwise rotation functions'''
def rotate_clockwise_about_z_axis(mypoint, hc_angle):
    
    rad_angle = ((180 - hc_angle)/2)*ang_to_rad

        
    x = mypoint[0]
    y = mypoint[1]
    z = mypoint[2]
    
    #use rotation matrix for Rz
    newx = (x*np.cos(rad_angle)) + (y*np.sin(rad_angle))  # + -
    newy = (-x*np.sin(rad_angle)) + (y*np.cos(rad_angle))  # + +
    newz = z
    
    return [newx, newy, newz]

def rotate_clockwise_about_y_axis(mypoint, roty_degree): #angle = 90
    
    rad_angle = roty_degree*ang_to_rad
    
    x = mypoint[0]
    y = mypoint[1]
    z = mypoint[2]
    
    #use rotation matrix for Ry
    newx = (x*np.cos(rad_angle)) - (z*np.sin(rad_angle))
    newy = y
    newz = (+x*np.sin(rad_angle)) + (z*np.cos(rad_angle))
    
    return [newx, newy, newz]

# def rotation_matrix(axis, theta):
#     """
#     Return the rotation matrix associated with counterclockwise rotation about
#     the given axis by theta radians.
#     """
#     axis = np.asarray(axis)
#     axis = axis / mt.sqrt(np.dot(axis, axis))
#     a = mt.cos(theta / 2.0)
#     b, c, d = -axis * mt.sin(theta / 2.0)
#     aa, bb, cc, dd = a * a, b * b, c * c, d * d
#     bc, ad, ac, ab, bd, cd = b * c, a * d, a * c, a * b, b * d, c * d
#     return np.array([[aa + bb - cc - dd, 2 * (bc + ad), 2 * (bd - ac)],
#                      [2 * (bc - ad), aa + cc - bb - dd, 2 * (cd + ab)],
#                      [2 * (bd + ac), 2 * (cd - ab), aa + dd - bb - cc]])

################   END    
#if __name__ == "__main__": main()
