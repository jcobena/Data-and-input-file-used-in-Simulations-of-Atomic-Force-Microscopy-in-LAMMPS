'''Create a hydrocarbon chain to be used as part of SAMs structures '''


'''Miscellaneous functions and constants'''

import numpy as np

ang_to_rad = 2*np.pi/360

    

def x_even_z_even(angle_bond, bond_length):
    cosangle_x = np.cos((180-angle_bond)*ang_to_rad)
    sinagle_z  = np.sin((180-angle_bond)*ang_to_rad)
    x_even = cosangle_x*bond_length
    z_even = sinagle_z*bond_length

    #x_even is the projection of the hypotenusa (bl) on the x-axis
    #y_even is the projection of the hypotenusa (bl) on the y-axis

    return [x_even, z_even]
    
