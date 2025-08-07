import collections 
from collections import Counter
import numpy as np
import copy
import operator 
from a_datafile_functions import *
from a_driver_main import *


fungroup = 'NH2'
tip = 'diamond'

x = 21.5
y = 27.5

dimensions = [
    [0, 49.9482, 'xlo', 'xhi'],
    [0, 51.9076, 'ylo', 'yhi'],
    [1.45, 88.1, 'zlo', 'zhi'],
    ]


#================================ new types ========================================

#first group:
new_vat_type      = 1 #virtual atom
new_tip_type      = 2 #C in diamond tip  
new_au_type       = 3 #Au 
new_h_watr_type   = 4 #H in H2O
new_o_watr_type   = 5 #O in H2O

#second group   
new_S_type        = 6  #S 
new_CH2_type      = 7  #CH2     
new_CH2_10        = 8  #third CH2 counting from the N in NH2 
new_CH2_11        = 9  #CH2 next to the NH2. C11 
new_N_inNH2       = 10 #N in NH2
new_H_inNH2       = 11 #Hs in NH2

#charges
q_S = 0
q_CH2 = 0
q_C10 = 0
q_C11 = 0.18 #pending
q_N_inNH2  = -0.9 #-0.837 #pending
q_H_inNH2  = 0.36 #pending
q_Au = 0
q_Ow = -0.8476
q_Hw = 0.4238

#bond types
S_C_type = 1
C_C_type = 2
C10_C11_type = 3
C11_N_in_NH2_type = 4
N_H_type = 5
O_H_water_type = 6

#angle types:
s_ch2_ch2_type = 1
ch2_ch2_ch2_type = 2
c10_c11_n_type = 3
c11_n_h_type = 4
n_h_h_type  = 5
h_o_h_water_type = 6

#diheral types:
s_c_c_c_type = 1


#masses:
masses_list = [
    [1, 12.0],      # virtual atom
    [2, 12.0],      #C in tip
    [3, 196.96657], #Au
    [4, 1.0],       #H in water
    [5, 16.0],      #O in water
    [6, 32.065],    #S        
    [7, 14.0],    #CH2        
    [8, 14.0],    #C10
    [9, 14.0],    #C11
    [10, 14.0],   #N in NH2
    [11, 1.0]]    #H in NH2
#number of sams chains 
nchains = 120

#for water positions:
n_x = 16
n_y = 16
n_z = 18 #19


filetestx = open('mytest', 'w') 

def main():
    
    for i in range(1, 2, 1):
        driver(i)
        

def driver(numberdatafile):
    
    with open('vat.txt', 'r') as vatfile, open('the_tip.txt', 'r') as tipfile, \
         open('gold_bottom.txt', 'r') as goldbottomfile, \
         open('datafile-{}-{}-{}-{}-{}.data'.format(tip, fungroup, x, y, numberdatafile), 'w') as outfile:
    
        '''============ Get vat, tip and bottom layer of gold from the files ========================='''
        vat_list, tip_list, gold_a_list = ([] for i in range(3))
        
        for i in vatfile:
            j = i.split()
            k = [int(j[0]), int(j[1]), int(j[2]), float(j[3]), float(j[4]), float(j[5]), float(j[6])]
            vat_list.append(k)
            
        for i in tipfile:
            j = i.split()
            k = [int(j[0]), int(j[1]), int(j[2]), float(j[3]), float(j[4]), float(j[5]), float(j[6])]
            tip_list.append(k)
            
        for i in goldbottomfile:
            j = i.split()
            k = [int(j[0]), int(j[1]), int(j[2]), float(j[3]), float(j[4]), float(j[5]), float(j[6])]
            gold_a_list.append(k)
       
        ''' ====================== Get the rest of the gold ==========================================='''
        
        mylayer_b_gold = layer_b_scratch(gold_a_list)
            
        mylayer_c_gold = layer_c_scratch(gold_a_list, mylayer_b_gold)
        
        gold_atoms_first_layer_ABC = gold_a_list + mylayer_b_gold + mylayer_c_gold
        
        gold_atoms_second_layer_ABC = get_second_ABC_gold_layer(gold_atoms_first_layer_ABC)
        
        mygold_all = gold_atoms_first_layer_ABC + gold_atoms_second_layer_ABC
        
        ''' This section calls the SAMS with the CH3 termination function. The SAMs will be created, rotated if
        needed and placed on the gold plate'''
        
        ''' ========================= Sulfur coords ============================================='''
        #we need the S coordinates, so we need some gold atoms coords 
        sulfur_coords = get_coords_sulfur()
        
        start_id_sams = mygold_all[-1][0]+1 #my_updated_water_mols[-1][0]+1
        start_mol_tag_sams = 1 #my_updated_water_mols[-1][1]+1
        
        print(len(sulfur_coords))
        
        ''' =========================== SAMs coords ============================================='''
        sams_finished = really_get_sams(sulfur_coords, start_id_sams, start_mol_tag_sams)
        
        '''=============================== SAMs Bonds, angles and dihedrals================================= '''
            
        my_sams_bonds = sams_bonds(sams_finished, 1)#my_water_bonds[-1][0]+1)
        my_sams_angles = sams_angle(sams_finished, 1) #my_water_angles[-1][0]+1)
        my_sams_dihedrals = sams_dihedral(sams_finished, 1)

        ''' ==================================== water positions =========================================== '''
        #oxygen positions
        my_better_oxygens = water_oxygen_positions() 
        #print(my_better_water)
        
        #hydrogen creation
        noxygens = n_x*n_y*n_z
        print('noxygens = ', noxygens)
        my_better_h1s = adding_H1(noxygens, 109.47, 1.0) #4445 
        my_better_h2s = adding_H2(my_better_h1s, 109.47, 1.0)
        
        #hydrogen positions:
        h1_final = final_position_for_H(my_better_oxygens, my_better_h1s) #move the water atoms from the origin to the random positions generated in oxygens
        h2_final = final_position_for_H(my_better_oxygens, my_better_h2s)
        
        start_water_id = sams_finished[-1][0]+1
        start_mol_tag = sams_finished[-1][1] + 1
        all_water = final_water(my_better_oxygens, h1_final, h2_final, start_water_id, start_mol_tag)
        
        # for i in all_water:
        #     print(*i, file = filetestx)
        
        
        ''' ==================================== all atoms =================================================''' 
        #first, delete the water that overlaps the tip, 355 must be destroyed
        all_atoms =  vat_list + tip_list + mygold_all + sams_finished + all_water
        
        for i in all_atoms:
            print(i, file=filetestx)
        
        print('water before deleting = ', len(all_water))
        
        water_list_updated = get_tip_shield_and_delete_water(all_atoms, all_water)
        
        all_atoms_updated = vat_list + tip_list + mygold_all + sams_finished + water_list_updated
        

        ''' ============================== All bonds, angles and dihedrals =================================='''
        #create water bonds:
        bond_id_start = my_sams_bonds[-1][0]+1
        my_water_bonds = generate_water_bonds(water_list_updated, bond_id_start)
        
        #create water angles:
        angle_id_start = my_sams_angles[-1][0]+1
        my_water_angles =  generate_water_angles(water_list_updated, angle_id_start)
        
        all_bonds = my_sams_bonds + my_water_bonds
        all_angles = my_sams_angles + my_water_angles
        all_dihedrals = my_sams_dihedrals 
        
        ''' ========================== Generate the data file ==============================================''' 
        
#         mydatafile = generate_datafile(all_atoms, 
#             all_bonds, all_angles, all_dihedrals, dimensions, masses_list)
        
        mydatafile = generate_datafile(all_atoms_updated, 
            all_bonds, all_angles, all_dihedrals, dimensions, masses_list)
        
        for i in mydatafile:
            for j in i:
                print(*j,  file = outfile)
        
        '''======================================  END =========================================================='''

''' The functions below prepare the data file and modify, add, change stuff... '''
                
def get_tip_shield_and_delete_water(atomlist, water_list):
    
    start_id = water_list[0][0]
    
    safety_clearance_x = 3.2
    safety_clearance_y = 4
    safety_clearance_z = 2.2
    
    only_tip = []
    for i in atomlist:
        if i[2] == new_tip_type:
            #print(i)
            only_tip.append(i)
            
    # for j in only_tip:
    #     print(j, file=filetestx)
    
    #for i in only_tip:
    min_x = sorted(only_tip, key=operator.itemgetter(4), reverse=False)[0][4] - safety_clearance_x
    max_x = sorted(only_tip, key=operator.itemgetter(4), reverse=True)[0][4] + safety_clearance_x
    
    min_y = sorted(only_tip, key=operator.itemgetter(5), reverse=False)[0][5] - safety_clearance_y
    max_y = sorted(only_tip, key=operator.itemgetter(5), reverse=True)[0][5] + safety_clearance_y
    
    min_z = sorted(only_tip, key=operator.itemgetter(6), reverse=False)[0][6] - safety_clearance_z
    max_z = sorted(only_tip, key=operator.itemgetter(6), reverse=True)[0][6] + safety_clearance_z
    
    print('min_x = ', min_x, 'max_x = ', max_x)
    print('min_y = ', min_y, 'max_y = ', max_y)
    print('min_z = ', min_z, 'max_z = ', max_z)
    
    water_molecules_to_delete = []
    for i in water_list:
        if (i[4] > min_x and i[4] < max_x) and (i[5] > min_y and i[5] < max_y) and (i[6] > min_z and i[6] < max_z):
            if i[1] not in water_molecules_to_delete:
                water_molecules_to_delete.append(i[1])
            else:
                continue
        
    print(len(water_molecules_to_delete))
    
    water_list_updated = []
    for i in water_list:
        if i[1] in water_molecules_to_delete:
            pass
        else:
            water_list_updated.append(i)
            
    #print('len updated = ', len(water_list_updated))
    
    print('initial number of water molecules = ', len(water_list)/3)
    print('new number of water after deletion = ', len(water_list_updated)/3) #water that was on the tip
    
    #when I ran this it showed 4480, I need to bring down the number to 4445 from 4480
    # so just chose randomnly the molecues to be deleted.
    goal_val = 4000 #4300, 4445
    
    currentvalue = len(water_list_updated)/3 # number of current water molecules
    diff_value = int(currentvalue - goal_val)*3 #number of atoms to be deleted
    print('diff value = ',diff_value)
    
    #print(diff_value)
    
    sliced_list = water_list_updated[:-diff_value]
    
    print('len of water in atoms final = ', len(sliced_list))
    
    
    final_water_l = []
    for i in sliced_list:
        coords = [start_id, i[1], i[2], i[3], round(i[4], 4), round(i[5], 4), round(i[6],4)]
        final_water_l.append(coords)
        start_id = start_id + 1
        
    print('final number of water molecules = ', len(final_water_l)/3)
        
    return final_water_l
        
                
def final_water(Olist, H1list, H2list, startid, nmol):
    
    Oready, h1ready, h2ready  = ([] for i in range(3))
    
    my_nmol = nmol
    id0 = startid
    for o in Olist:
        #print(o)
        coords = [id0, my_nmol, new_o_watr_type, q_Ow, o[0], o[1], o[2]]
        Oready.append(coords)
        id0 = id0 + 3
        my_nmol = my_nmol + 1
    
    my_nmol = nmol
    id0 = startid + 1
    for h1 in H1list:
        coords = [id0, my_nmol, new_h_watr_type, q_Hw, h1[0], h1[1], h1[2]]
        h1ready.append(coords)
        id0 = id0 + 3
        my_nmol = my_nmol + 1
    
    my_nmol = nmol
    id0 = startid + 2
    for h2 in H2list:
        coords = [id0, my_nmol, new_h_watr_type, q_Hw, h2[0], h2[1], h2[2]]
        h2ready.append(coords)
        id0 = id0 + 3
        my_nmol = my_nmol + 1
        
    all_water = Oready + h1ready + h2ready
        
    water_sorted = sorted(all_water, key=operator.itemgetter(0))
    
    for i in water_sorted:
        x = round(i[4], 4)
        y = round(i[5], 4)
        z = round(i[6], 4)
        
        i[4] = x
        i[5] = y
        i[6] = z
    return water_sorted
                
def water_oxygen_positions():
    
#     dimensions = [
#     [0, 49.9482, 'xlo', 'xhi'],
#     [0, 51.9076, 'ylo', 'yhi'],
#     [1.45, 88.1, 'zlo', 'zhi'],
#     ]
    
    xstart = 2.66171
    xend = 47.28649
    
    ystart = 2.66171 #33.4385
    yend = 51.9076
    
    zstart = 36.3308 #33.4385
    zend = 85.5
    
    bond_length = 0.9572
    angle = 104.52
    safety_clearance = 1.0
    
    #this produces 4800 positions (some will be deleted because of the tip)
    x_pos = np.linspace(xstart, xend, n_x)
    y_pos = np.linspace(ystart, yend, n_y)
    z_pos = np.linspace(zstart, zend, n_z)
    
    new_oxygens = []
    for x in x_pos:
        for y in y_pos:
            for z in z_pos:
                coords = [x, y, z]
                new_oxygens.append(coords)
    
    return new_oxygens

    #create the H1
def adding_H1(number_of_atoms, angle, bond_length):
    bond_length = 1.0
    angle = 109.47
    rand_rot = np.random.randint(1, 3)
    values = [0,0,0]
    totalh1 = []
    #sprint(values)
    for _numbers_atoms in range((1*number_of_atoms)+1):
        rand_ini = np.random.randint(1, 7)
        
        rand_ini = 1
        rand_rot = 2
        
        if rand_ini == 1 and rand_rot == 1:
            basic_vector = [0.0, 0.0, bond_length]
        
        elif rand_ini == 1 and rand_rot == 2:
            basic_vector = [0.0, 0.0, bond_length]
        #===========================================================================
        elif rand_ini == 2 and rand_rot == 1:
            basic_vector = [0.0, bond_length, 0.0]
            
        elif rand_ini == 2 and rand_rot == 2:
            basic_vector = [0.0, bond_length, 0.0]

        #===============================================================================
        elif rand_ini == 3 and rand_rot == 1:
            basic_vector = [bond_length, 0.0, 0.0]
            
        elif rand_ini == 3 and rand_rot == 2:
            basic_vector = [bond_length, 0.0, 0.0]
                
        ##########################=====================================================
        elif rand_ini == 4 and rand_rot == 1:
            basic_vector = [0.0, 0.0, -bond_length]
            
        elif rand_ini == 4 and rand_rot == 2:
            basic_vector = [0.0, 0.0, -bond_length]
            
        #===========================================================================
        elif rand_ini == 5 and rand_rot == 1:
            basic_vector = [0.0, -bond_length, 0.0]
            
        elif rand_ini == 5 and rand_rot == 2:
            basic_vector = [0.0, -bond_length, 0.0]
            
        #===============================================================================
        elif rand_ini == 6 and rand_rot == 1:
            basic_vector = [-bond_length, 0.0, 0.0]
            
        elif rand_ini == 6 and rand_rot == 2:
            basic_vector = [-bond_length, 0.0, 0.0]
             
        totalh1.append(basic_vector)
    return totalh1

def adding_H2(h1_list, angle, bond_length):
    r_angle = mt.radians(angle)
    cos_ang = np.cos(r_angle)
    total_h2 = []

    for basic_vector in h1_list:
        #print(basic_vector)
        
        rand_rot = np.random.randint(1, 3)
        
        rand_rot = 1

        if basic_vector[0] == 0.0 and basic_vector[1] == 0.0:
            #print(basic_vector)
            if rand_rot == 1:
                xh2 = basic_vector[0]
                yh2 = (basic_vector[1]*np.cos(r_angle) - basic_vector[2]*np.sin(r_angle))
                zh2 = (basic_vector[1]*np.sin(r_angle) + basic_vector[2]*np.cos(r_angle))
            elif rand_rot == 2:
                xh2 = basic_vector[0]*np.cos(r_angle) + basic_vector[2]*np.sin(r_angle)
                yh2 = basic_vector[1]
                zh2 = -basic_vector[0]*np.sin(r_angle) + basic_vector[2]*np.cos(r_angle) 
                
        elif basic_vector[0] == 0.0 and basic_vector[2] == 0.0:
            if rand_rot == 1:
                xh2 = basic_vector[0]
                yh2 = (basic_vector[1]*np.cos(r_angle) - basic_vector[2]*np.sin(r_angle))
                zh2 = (basic_vector[1]*np.sin(r_angle) + basic_vector[2]*np.cos(r_angle))
            elif rand_rot == 2:
                xh2 = basic_vector[0]*np.cos(r_angle) - basic_vector[1]*np.sin(r_angle)
                yh2 = basic_vector[0]*np.sin(r_angle) + basic_vector[1]*np.cos(r_angle)
                zh2 = basic_vector[2]
                
        elif basic_vector[1] == 0.0 and basic_vector[2] == 0.0:
            if rand_rot == 1:
                xh2 = basic_vector[0]*np.cos(r_angle) + basic_vector[2]*np.sin(r_angle)
                yh2 = basic_vector[1]
                zh2 = -basic_vector[0]*np.sin(r_angle) + basic_vector[2]*np.cos(r_angle) 
                
            elif rand_rot == 2:
                xh2 = basic_vector[0]*np.cos(r_angle) - basic_vector[1]*np.sin(r_angle)
                yh2 = basic_vector[0]*np.sin(r_angle) + basic_vector[1]*np.cos(r_angle)
                zh2 = basic_vector[2]
                
        vector_h2 = [xh2, yh2, zh2]
        total_h2.append(vector_h2)
        
    return total_h2

def final_position_for_H(oxy_list, a_list):
    final_position_of_a = []
    for a1, a0 in zip(oxy_list, a_list):
        x = a1[0] + a0[0]
        y = a1[1] + a0[1]
        z = a1[2] + a0[2]
        coords = [x, y, z]
        final_position_of_a.append(coords)
        
    return final_position_of_a
                
def really_get_sams(Scoords, id0, mol_tag):
    
    #    #first group:
    # new_vat_type      = 1 #virtual atom
    # new_tip_type      = 2 #C in diamond tip  
    # new_au_type       = 3 #Au 
    # new_h_watr_type   = 4 #H in H2O
    # new_o_watr_type   = 5 #O in H2O
    #
    # #second group   
    # new_S_type        = 6  #S 
    # new_CH2_type      = 7  #CH2     
    # new_CH2_10        = 8  #third CH2 counting from the N in NH2 
    # new_CH2_11        = 9  #CH2 next to the NH2. C11 
    # new_N_inNH2       = 10 #N in NH2
    # new_H_inNH2       = 11 #Hs in NH2
    #
    # #charges
    # q_S = 0
    # q_CH2 = 0
    # q_C10 = 0
    # q_C11 = 'pending'
    # q_N_inNH2  = 'pending'
    # q_H_inNH2  = 'pending'
    # q_Au = 0
    # q_Ow = -0.8476
    # q_Hw = 0.4238
        
    sams_chains = []
    
    for i in Scoords:
        a_sams_chain = driver_NH2()
        xn = i[0]
        yn = i[1]
        zn = i[2]
        
        for j in a_sams_chain:
            
            my_x = round(j[0] + xn, 4)
            my_y = round(j[1] + yn, 4)
            my_z = round(j[2] + zn, 4)
            
            if j[3] == 'S':
                atype = new_S_type
                qtype = q_S
                
            elif j[3] == 'CH2':
                atype = new_CH2_type
                qtype = q_CH2
                
            elif j[3] == 'C10':
                atype = new_CH2_10
                qtype = q_C10
                
            elif j[3] == 'C11':
                atype = new_CH2_11
                qtype = q_C11
                
            elif j[3] == 'NinNH2':
                atype = new_N_inNH2
                qtype = q_N_inNH2
                
            elif j[3] == "HinNH2":
                atype = new_H_inNH2
                qtype = q_H_inNH2
                
            coords = [id0, mol_tag, atype, qtype, my_x, my_y, my_z]
                
            id0 =  id0 + 1
            sams_chains.append(coords)
            
        mol_tag = mol_tag + 1
        
    return sams_chains
                
            
def sams_dihedral(sams_atoms, dihedral_id):
    #apparently there is only 1 dihedral type here
    start_sams = sams_atoms[0][0]
    end_sams = sams_atoms[-1][0] 
    
    sams_dihedral_list = []
    atom_counter = 1
    
    dihedral_type = 1
    
    print('start and end sams =', start_sams, end_sams)
    
    dihe_counter = 0
    i = start_sams
    while i <= end_sams-4:
    
    #for i in range(start_sams, end_sams-2, 1):
         #atom 10 is the first ch2 in ch2-ch2-ch2-ch3, so we need to skip the other 2 ch2 and the ch3
            
        new_dihedral = [dihedral_id, dihedral_type, i, i+1, i+2, i+3]
        sams_dihedral_list.append(new_dihedral)
        
        i = i + 1
        dihe_counter = dihe_counter + 1
        dihedral_id = dihedral_id + 1
        
        if dihe_counter == 7:
            dihe_counter = 0
            i = i + 3
            
    return sams_dihedral_list

    
def sams_angle(sams_atoms, angle_id):
    
    #angle types:
    # s_ch2_ch2_type = 1
    # ch2_ch2_ch2_type = 2
    # c10_c11_n_type = 3
    # c11_n_h_type = 4
    # n_h_h_type  = 5
    # h_o_h_water_type = 6
    
    start_sams = sams_atoms[0][0]
    end_sams = sams_atoms[-1][0] 
    
    sams_angles_list = []
    
    atom_counter = 1
    
    for i in range(start_sams, end_sams, 1):
        
        if atom_counter == 1:
            atom_counter = atom_counter + 1
            continue
            #angtype = s_ch2_ch2_type
            
        elif atom_counter == 2:
            angtype = s_ch2_ch2_type
            
        elif atom_counter > 2 and atom_counter <= 11 :
            angtype = ch2_ch2_ch2_type
            
        elif atom_counter == 12:
            angtype = c10_c11_n_type
            
        elif atom_counter == 13:
            angtype1 = c11_n_h_type
            
            newangle = [angle_id, angtype, i-1, i, i+1]
            sams_angles_list.append(newangle)
            
            newangle = [angle_id, angtype, i-1, i, i+2]
            sams_angles_list.append(newangle)
            
            newangle = [angle_id, angtype, i+1, i, i+2]
            sams_angles_list.append(newangle)
            
            
            atom_counter = atom_counter + 3
            angle_id = angle_id + 3
            continue 
            
            
        else:
            atom_counter = 1
            continue
                
        newangle = [angle_id, angtype, i-1, i, i+1]
        sams_angles_list.append(newangle)
        
        atom_counter = atom_counter + 1
        angle_id = angle_id + 1

    return sams_angles_list
    
            
def sams_bonds(sams_atoms, bond_id):
    
    #bond types
    # S_C_type = 1
    # C_C_type = 2
    # C10_C11_type = 3
    # C11_N_in_NH2_type = 4
    # N_H_type = 5
    # O_H_water_type = 6

    start_sams = sams_atoms[0][0]
    end_sams = sams_atoms[-1][0] 
    sams_bonds_list = []
    
    atom_counter = 1
    for i in range(start_sams, end_sams, 1):
        
        if atom_counter == 1:
            btype = S_C_type
            
        if atom_counter >1 and atom_counter < 11:
                btype = C_C_type
            
        if atom_counter == 11: #atom=11, but C=10
            btype = C10_C11_type
            
        if atom_counter == 12: #atom 11 but C=12
            btype = C11_N_in_NH2_type
            
        if atom_counter == 13:
            btype = N_H_type
            
            newbond = [bond_id, btype, i, i+1]
            sams_bonds_list.append(newbond)
            
            bond_id = bond_id + 1
            
            newbond = [bond_id, btype, i, i+2]
            sams_bonds_list.append(newbond)
            
            
            bond_id = bond_id + 1
            atom_counter = atom_counter + 1
            continue
            
        else:
            if atom_counter == 14:
                atom_counter = atom_counter + 1
                continue
            
            if atom_counter == 15:
                atom_counter = 1
                continue
                
            # else:
            #     #resetting to 1, going back to the top as 1 (because of continue)
            #     atom_counter = 1
            #     continue
                
        newbond = [bond_id, btype, i, i+1]
        sams_bonds_list.append(newbond)
        
        atom_counter = atom_counter + 1
        bond_id = bond_id + 1

    return sams_bonds_list
    
    
def get_coords_sulfur():
    zn = 15.7967
    x0_odd = 2.4974
    x0_even = 4.99481 # I think...
    
    y0_odd = 1.4419
    y0_even = 5.7675
    
    x_positions_odd = []
    x_positions_even = []
    
    step_x = 4.99481 #same for odd and even
    step_y = 8.6513
    
    def get_coords_xy(my_k, nlines, mystep):
        k_positions = []
        n = 1
        while n <= nlines:
            k_positions.append(round(float(my_k), 4))
            my_k = my_k + mystep
            n = n + 1
        
        return k_positions

    xodds = get_coords_xy(x0_odd, 10, step_x)
    xevens = get_coords_xy(x0_even, 10, step_x)
    
    yodds= get_coords_xy(y0_odd, 6, step_y)
    yevens = get_coords_xy(y0_even, 6, step_y)
    
    def combine_coords_odd_even(xlist, ylist):
        some_coords = []
        for x in xlist:
            for y in ylist:
                acoord = [x, y, zn] 
                some_coords.append(acoord)
                #id0 = id0 + 1
        return some_coords
    
    odd_coords = combine_coords_odd_even(xodds, yodds)
    #next_id = odd_coords[-1][0] + 1
    even_coords = combine_coords_odd_even(xevens, yevens)
    
    layer_S = odd_coords + even_coords
    
    return layer_S
    

def get_second_ABC_gold_layer(gold_first_layer):
    
    z_add = 6.86017
    
    gold_sorted = sorted(gold_first_layer, key=operator.itemgetter(0))
    id_0 = gold_sorted[-1][0] + 1
    
    second_layer = []
    for i in gold_sorted:
        x = i[4]
        y = i[5]
        z = i[6]
        new_z = z + z_add
        second_layer_coords = [id_0, 0, new_au_type, q_Au, x, y, new_z]
        second_layer.append(second_layer_coords)
        id_0 = id_0 + 1
        
    return second_layer

def layer_c_scratch(layer_a, layer_b):
    z = 6.86017
    moltag = 0
    id_0 = layer_b[-1][0] + 1
    
    both_layers = layer_a + layer_b
    
    #get the single file where y=0--------------------------------------
    first_line_xs = []
    for i in both_layers:
        ys = i[5]
        if ys == 0:
            first_line_xs.append(i[4])

    first_line_xs_sorted = sorted(first_line_xs)
    #print(first_line_xs_sorted)
    
    #now that we all the x's at y=0 then  we can get the x's (midpoints) that correspond to the odds x's positions
    #so we need to form the pairs to get the midpoint
    at1 = first_line_xs_sorted[0::2]
    at2 = first_line_xs_sorted[1::2]
    
    midpoint_odd_xs_lst = []
    for i,j in zip(at1, at2):
        #print(i,j)
        midpoint_odd_xs = round(float(((j-i)/2.0) + i),4)
        midpoint_odd_xs_lst.append(midpoint_odd_xs)
        
    #print(midpoint_odd_xs_lst)#important result
    #-------------------------------------------------------------------------
    #now even file: move the midpoint_odd_xs_lst a little up in the y axis (1.44188) and 
    #a little to the right on the x-axis (4.16235)
    x_add = 4.16235 - 1.66495
    
    midpoint_even_xs_lst = []
    for i in midpoint_odd_xs_lst:
        xp = round(i + x_add,4)
        midpoint_even_xs_lst.append(xp)
    #-----------------------------------------------------------------------------------
    #now let's duplicate the odd and even along the y axis!
    #1.44188 is y0 for 
    
    #odd layers
    add_y = 2.88376
    
    def get_coords_k(my_k, nlines, mystep):
        k_positions = []
        n = 1
        while n <= nlines:
            k_positions.append(round(float(my_k), 4))
            my_k = my_k + mystep
            n = n + 1
        
        return k_positions
    
    odd_points_y = get_coords_k(0, 18, add_y)
    even_points_y = get_coords_k(1.44188, 18, add_y)
    #print(odd_points_y)
    #print(even_points_y)
    
    
    def combine_coords_odd_even_k(xlist, ylist, id0):
        some_coords = []
        for x in xlist:
            for y in ylist:
                acoord = [id0, moltag, new_au_type, q_Au, x, y, z] 
                some_coords.append(acoord)
                id0 = id0 + 1
        return some_coords
    
    odd_coords = combine_coords_odd_even_k(midpoint_odd_xs_lst, odd_points_y, id_0)
    
    next_id = odd_coords[-1][0] + 1
    even_coords = combine_coords_odd_even_k(midpoint_even_xs_lst, even_points_y, next_id)
    
    layer_c = odd_coords + even_coords
    
    return layer_c
        
        
def layer_b_scratch(layer_a):
#def layer_b(start_x, end_x, start_y, end_y, step_x):
    #print(layer_a[-1][0])
    id_0 = layer_a[-1][0] + 1
    
    moltag = 0
    z = 4.63208
    layer_c_list = []
    
    #---------------the x positions only--------------
    x0_odd = 0.83247
    x0_even = 3.32988
    
    x_positions_odd = []
    x_positions_even = []
    
    #step_x = 2.49741 #same for odd and even
    step_x = 4.99481 #same for odd and even
    
    def get_coords_xy(my_k, nlines, mystep):
        k_positions = []
        n = 1
        while n <= nlines:
            k_positions.append(round(float(my_k), 4))
            my_k = my_k + mystep
            n = n + 1
        
        return k_positions

    xodds = get_coords_xy(x0_odd, 10, step_x)
    xevens = get_coords_xy(x0_even, 10, step_x)
    #----------------------------------------------------
    #now the y's-----------------------------------------
    y0_odd = 1.44188
    y0_even = 0 #3.32988
    step_y = 2.88376
    
    yodds= get_coords_xy(y0_odd, 18, step_y)
    yevens = get_coords_xy(y0_even, 18, step_y)
    
    #----------------------------------------------------
    #combine the pairs:
    def combine_coords_odd_even(xlist, ylist, id0):
        some_coords = []
        for x in xlist:
            for y in ylist:
                acoord = [id0, moltag, new_au_type, q_Au, x, y, z] 
                some_coords.append(acoord)
                id0 = id0 + 1
        return some_coords
    
    odd_coords = combine_coords_odd_even(xodds, yodds, id_0)
    
    next_id = odd_coords[-1][0] + 1
    even_coords = combine_coords_odd_even(xevens, yevens, next_id)
    
    layer_b = odd_coords + even_coords
    
    return layer_b
        
def layer_a_scracth(aulist):
    
    autype = aulist[0][2]
    aumol = 0
    aucharge = aulist[0][3]
    
    au_sorted_z = sorted(aulist, key=operator.itemgetter(6))
    
    layer_a = au_sorted_z[0:360]
    
    #the line below is just to get the first id!!, then it's useless
    au_sorted_id = sorted(aulist, key=operator.itemgetter(0))
    
    #adjust enumeration of bottom_layer_a:
    n = au_sorted_id[0][0]
    for i in layer_a:
        i[0] = n
        n = n + 1    

    for i in layer_a:
        x = round(float(i[4]), 4)
        y = round(float(i[5]), 4)
        z = round(float(i[6]), 4)
        
        i[4] = x
        i[5] = y
        i[6] = z
    
    return layer_a


def generate_water_angles(water_list, start_angle_id):
    
    water_angles = []
    
    start_number = water_list[0][0]
    last_water_atom = water_list[-1][0]
    
    
    while start_number < last_water_atom:
#     for i in range(start_number, last_water_atom+1, 1):
        angleline = [start_angle_id, h_o_h_water_type, start_number+1, start_number, start_number+2]
        water_angles.append(angleline)
        start_angle_id = start_angle_id + 1
        start_number = start_number + 3

    return water_angles
    
    
def generate_water_bonds(water_list, bid_start):
    
    water_bonds = []
    
    start_number = water_list[0][0]
    last_water_atom = water_list[-1][0]
    
    
    while start_number <= last_water_atom:
#     for i in range(start_number, last_water_atom+1, 1):
        bondline = [bid_start, O_H_water_type, start_number, start_number+1]
        water_bonds.append(bondline)
        bid_start = bid_start + 1
        #print(bondline)
        bondline = [bid_start, O_H_water_type, start_number, start_number+2]
        #print(bondline)
        water_bonds.append(bondline)
        
        start_number = start_number + 3
        bid_start = bid_start + 1
        
    return water_bonds
        

    
if __name__ == "__main__": main()
