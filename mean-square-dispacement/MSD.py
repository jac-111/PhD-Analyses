#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 17 03:57:49 2021

This code collects the CoM coordinates from polar space and also a single selected atom. Originally,
this was intended for the analysis of the end to end distance and msd of POSS systems.
Requires HISTORY and FIELD. Is only Parameterised to run on pure systems

@author: jacob
"""
import os
import numpy as np

import time

start_time = time.time()


print("\nThe current Program running is: ", os.path.basename(__file__), "\n")

if not os.path.exists('./output_files/msd_data'):
    os.mkdir('./output_files/msd_data')

    # /* ------------------------------------------------------------------ */
    # /*               Collect information from FIELD file                  */
	# /* ------------------------------------------------------------------ */


field_file = open ('./input_files/FIELD', 'r')
print("\n>>>> FIELD FILE DATA <<<<\n")
for a, line in enumerate(field_file):
    if a > 3:
        if a == 4:
            num_mols_f = int(line.split()[1])
        if a ==5:
            num_atoms_f = int(line.split()[1])
            atom_type_f = [0 for i in range(num_atoms_f)]
            atom_mass_f = [0 for i in range(num_atoms_f)]
            atom_charge_f = [0 for i in range(num_atoms_f)]
            repeat_counter_f = [0 for i in range(num_atoms_f)]
            frzn_num_f = [0 for i in range(num_atoms_f)]
            element_list_f = [0 for i in range(num_atoms_f)]
            print("THE NUMBER OF ATOMS IN EACH MOLECULE IS:\n%d\n" % num_atoms_f)
    if a > 5:
        if a <= 5 + num_atoms_f :

            atom_type_f[a-6] = line.split()[0]
            atom_mass_f[a-6] = float(line.split()[1])
            atom_charge_f[a-6] = float(line.split()[2])
            repeat_counter_f[a-6] = int(line.split()[3])
            frzn_num_f[a-6] = int(line.split()[4])
        if a > 5 + num_atoms_f :
            break            
field_file.close()

masses = np.array(atom_mass_f).reshape(-1,1)
all_masses = np.tile(masses, (num_mols_f ,1))
mol_mass = np.sum(masses)
inv_mol = 1/mol_mass


    # /* ------------------------------------------------------------------ */
    # /*         Collect First line information from history file           */
	# /* ------------------------------------------------------------------ */


history_file = open('./input_files/HISTORY', 'r')
for i, line in enumerate(history_file):
    if i == 1:
         traj_key = int(line.split()[0])
         boundary_key = int(line.split()[1])
         num_atoms_h = int(line.split()[2])
         num_frames_h = int(line.split()[3])
         num_lines_h = int(line.split()[4])
         break


    # /* ------------------------------------------------------------------ */
    # /*       Account for periodicity in system + calculate the CoM        */
	# /* ------------------------------------------------------------------ */


np.set_printoptions(suppress=True)


def calc_CoM(mass_xyz, dimens):
    p_coords = np.array([[0 for l in range(3)] for k in range(len(mass_xyz))]).astype(float)
    p_coords[:,0] = ((mass_xyz[:,1] + dimens[0]/2) / dimens[0]) *(2*np.pi)
    p_coords[:,1] = ((mass_xyz[:,2] + dimens[1]/2) / dimens[1]) *(2*np.pi)
    p_coords[:,2] = ((mass_xyz[:,3] + dimens[2]/2) / dimens[2]) *(2*np.pi)
    

    eta = np.cos(p_coords)
    zeta = np.sin(p_coords)
    
    bar_eta = np.array([0 for i in range(3)]).astype(float)
    bar_eta[0] = inv_mol * np.sum(mass_xyz[:,0]*eta[:,0])
    bar_eta[1] = inv_mol * np.sum(mass_xyz[:,0]*eta[:,1])
    bar_eta[2] = inv_mol * np.sum(mass_xyz[:,0]*eta[:,2])
    
    bar_zeta = np.array([0 for i in range(3)]).astype(float)
    bar_zeta[0] = inv_mol * np.sum(mass_xyz[:,0]*zeta[:,0])
    bar_zeta[1] = inv_mol * np.sum(mass_xyz[:,0]*zeta[:,1])
    bar_zeta[2] = inv_mol * np.sum(mass_xyz[:,0]*zeta[:,2])
    
    bar_theta = np.array([0 for l in range(3)]).astype(float)
    bar_theta = np.arctan2(-bar_zeta, -bar_eta)  
    
    CoM = np.array([0 for i in range(3)]).astype(float)
    CoM[0] = dimens[0] * (bar_theta[0] / (2*np.pi)) 
    CoM[1] = dimens[1] * (bar_theta[1] / (2*np.pi)) 
    CoM[2] = dimens[2] * (bar_theta[2] / (2*np.pi)) 
    
    return(CoM)



    # /* ------------------------------------------------------------------ */
    # /*          Bring molecule to origin and make whole w/ pbc            */
	# /* ------------------------------------------------------------------ */
    
    

def make_whole(xyz, com, dimens):
    differ = xyz - com
    
    for j in range(3):    
        for i in range(num_atoms_f):
            
            if differ[i,j] <= -dimens[j]/2:
                differ[i,j] += dimens[j]
                
            if differ[i,j] >= dimens[j]/2:
                differ[i,j] -= dimens[j]
    
    
    
    return(differ)


    # /* ------------------------------------------------------------------ */
    # /*                   DECLARE MAIN LOOP VARIABLES                      */
	# /* ------------------------------------------------------------------ */
    

sys_coords = np.array([[0 for l in range(3)] for k in range(num_atoms_f*num_mols_f)]).astype(float)
mol_coords_temp = np.array([[0 for l in range(3)] for k in range(num_atoms_f)]).astype(float)

cell_dimensions = np.array([0 for l in range(3)]).astype(float)
all_dimensions = np.array([[0 for l in range(3)] for k in range(num_frames_h)]).astype(float)
all_coms = np.array([[0 for l in range(3)] for k in range(num_mols_f)]).astype(float)
my_hist_dimen = np.array([[0 for l in range(3)] for k in range(3)]).astype(float)

coms_3d = np.array([[[0 for l in range(3)] for k in range(num_mols_f)] for j in range(num_frames_h)]).astype(float)
chosen_3d = np.array([[[0 for l in range(3)] for k in range(num_mols_f)] for j in range(num_frames_h)]).astype(float)


print("Main Parsing Loop\n\nNumber of frames : ", num_frames_h)
print("Number of molecules: ", num_mols_f)

    # /* ------------------------------------------------------------------ */
    # /*                        MAIN PARSING LOOP                           */
	# /* ------------------------------------------------------------------ */
    
iter_num = (num_atoms_h*2) + 4
p =-1
dimens_sup = 4
mol_count = 1
counter = 0                     # COUNTS NUMBER OF MOLECULES(i.e loop 0 - num_mols)
c = 0                           # counts number of atoms in each molecule (i.e. loop 0 - num_atoms_f)
count_atoms_fr = 0
tick = 0                        # COUNTS FOR EACH ATOM IN EACH FRAME (i.e num_mol * num_atoms * num_frames)
frame_count = 0




for i, line in enumerate(history_file):
    
    
    if i % iter_num == 0 :                                                                  # FOR EACH FRAME, MAKE TOP LINE FOR
        m = 0                                                                               # COM MY_HISTORY
        p += 1
        top_line = line
        top_line_new = top_line.replace(' '+str(num_atoms_h)+' ', ' '+str(num_mols_f)+' ')
        
        
    if i == ((p*iter_num) + 1):                                                             # COLLECT X CELL DIMENSION
        x_dim = float(line.split()[(0)])
        cell_dimensions[0] = x_dim                                                          # X CELL FOR COM CALCULATION
        my_hist_dimen[0,0] = x_dim                                                          # X CELL FOR MY_HISTORY PRINTING
        all_dimensions[frame_count,0]  = x_dim
        
    if i == ((p*iter_num) + 2):
        y_dim = float(line.split()[(1)])
        cell_dimensions[1] = y_dim                                                          # Y CELL FOR COM CALCULATION
        my_hist_dimen[1,1] = y_dim                                                          # Y CELL FOR MY_HISTORY PRINTING
        all_dimensions[frame_count,1]  = y_dim
        
        
    if i == ((p*iter_num) + 3):
        z_dim = float(line.split()[(2)])
        cell_dimensions[2] = z_dim                                                          # Z CELL FOR COM CALCULATION
        my_hist_dimen[2,2] = z_dim                                                          # Z CELL FOR MY_HISTORY PRINTING
        all_dimensions[frame_count,2]  = z_dim
        
        
    
        
    if i >= ((p * iter_num) + 4) and i % 2 == 1:                                            # COLLECTING THE XYZ COORDINATES
                                                                                            # OF INDIVIDUAL MOLECULES
        mol_coords_temp[c,m] = (float(line.split()[0]))                                     # COLLECT X COORDINATES
        sys_coords[count_atoms_fr,m] = (float(line.split()[0]))
        m += 1                                                                              # m COUNTS 0-1-2 for X-Y-Z
        
        mol_coords_temp[c,m] = (float(line.split()[1]))                                     # COLLECT Y COORDINATES
        sys_coords[count_atoms_fr,m] = (float(line.split()[1]))
        m += 1
        
        mol_coords_temp[c,m] = (float(line.split()[2]))                                     # COLLECT Z COORDINATES
        sys_coords[count_atoms_fr,m] = (float(line.split()[2]))
        m = 0
        c += 1
        tick += 1                                                                           # TICK COUNTS FOR EACH ATOM IN EACH
        count_atoms_fr += 1

        
    if i == (mol_count) * (num_atoms_f*2) + dimens_sup -1:                                  # IN CONTEXT OF 1 MOLECULE!
                                                                                            # STOPS WHEN END OF MOLECULE HAS BEEN
                                                                                            # MET
        m_xyz = np.concatenate((masses,mol_coords_temp), axis=1)                            # MAKE ARRAY WITH MASSES AND MOL 
                                                                                            # COORDS FOR 1 MOLECULE
        
        CoM = calc_CoM(m_xyz, cell_dimensions)                                              # CALCULATE COM AND PUT IN 1D ARRAY
        all_coms[counter] = CoM                                                             # COLLATE ALL COMS
        
        whole_xyz_dif = make_whole(m_xyz[:,1:4], CoM, cell_dimensions)                      # MOVE MOLECULE TO ORIGIN - INCORP
                                                                                            # PBC RULE
        mol_count += 1                                                                      # COUNTS FOR EVERY MOLECULE IN EVERY
                                                                                            # FRAME                     
        counter += 1                                                                        # COUNTS NUM MOLS IN EACH FRAME
        
        m=0                                                                                 # RESET XYZ COUNTER
        c=0                                                                                 # RESET XYZ COUNTER

        if tick % num_atoms_h == 0:                                                         # PRINTS THE COM_HISTORY FILE
            coms_3d[frame_count] = all_coms                                                 # ASSIGN FRAME OF COM TO 3d HOLDING ARRAY


            dimens_sup +=4
            counter = 0
            frame_count +=1
            count_atoms_fr = 0
            print(frame_count)

            


    # /* ------------------------------------------------------------------ */
    # /*         This is a fairly quick but reliable MSD code.              */
    # /*          MSD(tau) = < | x(tau)(i) - x(t-tau) |^2 >                 */
    # /* ------------------------------------------------------------------ */

def tau_3(xyz_1, xyz_0, dimens, fram):
    
    
    differ = xyz_1 - xyz_0                                  # This section is necessary for handling the periodic boundary condition
    for k in range(len(differ)):                            # all CoMs are referenced from position 0 CoM
        for j in range(3):                                  # All x y z are coordinates in the 'differ' array.
            for i in range(num_mols_f):
                    
                if differ[k,i,j] <= -dimens[k,j]/2:
                    differ[k,i,j] += dimens[k,j]
                    
                if differ[k,i,j] >= dimens[k,j]/2:
                    differ[k,i,j] -= dimens[k,j]
    
    
    
    s1 = np.array([[0 for l in range(3)] for k in range(num_mols_f)]).astype(float)     # declare 'step1' array
    msd2 = np.array([[0 for l in range(3)] for k in range(num_mols_f)]).astype(float)   #   

    
    MSD_a = np.array([0 for y in range(num_frames_h)]).astype(float)                    # final MSD array
    tau = np.array(range(1, (num_frames_h)))                                            # array for simplicity - all tau values
    
    for k in range(num_frames_h-1):
        print(k)
        for g in range(tau[-k-1]):
            
            s1 = (differ[tau[k] + g ] - differ[g])              #  x(tau)(i) - x(t-tau)                 
            mod1 = np.square(s1)                                #
            mod2 = np.sqrt(mod1)                                # modulus || section
            msd1 = np.square(mod2)                              # ^2
            msd2 = msd1.sum(axis=1)                             # sum 
            msd3 = np.mean(msd2)                                # average
        
        MSD_a[k+1] = msd3
            
        s1.fill(0)

    
    MSD_b = MSD_a.reshape(-1,1)
    index = np.array(range(num_frames_h)).astype(int).reshape(-1,1)
    ind_msd = np.append(index, MSD_b, axis=1)
    
    
    fee = open('./output_files/msd_data/msd_CoM_avg2.txt', 'w')
    np.savetxt(fee, ind_msd, delimiter='\t',fmt='%s')
    fee.close()
    return(MSD_a)



print("Managing 3d arrays - msds")

com_com_msd2 = tau_3(coms_3d, coms_3d[0], all_dimensions, num_frames_h)

print("Program finished")



