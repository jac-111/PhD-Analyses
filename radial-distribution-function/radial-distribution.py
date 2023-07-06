#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 11 18:50:20 2020

@author: jacob
"""




import numpy as np
import time
import json
import sys
#from tqdm import tqdm
import math
import os.path
#import matplotlib.pyplot as plt
print("\nThe current Program running is: ", os.path.basename(__file__), "\n")

dr = 0.2

m = 0
n = 0
my_history_file = open("./output_files/my_history_MOI4", 'r')
for i, line in enumerate(my_history_file):
    if i == 1:
         traj_key = int(line.split()[0])
         boundary_key = int(line.split()[1])
         num_atoms_h = int(line.split()[2])
         num_frames_h = int(line.split()[3])
         num_lines_h = int(line.split()[4])
         break

my_history_file.close()

iter_num_ = (num_atoms_h*2) + 4
all_dimensions = np.array([[0 for l in range(3)] for k in range(num_frames_h)]).astype(float)
    
my_history_2 = open("./output_files/my_history_MOI4", 'r')
for i, line in enumerate(my_history_2):
    if i >= 1:
        if i == (m*iter_num_) + 3:
            all_dimensions[m,n] = float(line.split()[(0)])
            n =n+1
        if i == (m*iter_num_) + 4:
            all_dimensions[m,n] = float(line.split()[(1)])
            n =n+1
        if i == (m*iter_num_) + 5:
            all_dimensions[m,n] = float(line.split()[(2)])
            m = m + 1
            n = 0
    
my_history_2.close()

xyz_min = np.amin(all_dimensions)
#xyz_min = np.floor(xyz_min)
print('smallest xyz box dimension = ', xyz_min)

dr       = dr / xyz_min
num_dr   = math.floor(0.5/dr)
r_max    = dr * num_dr

my_hist = np.zeros(num_dr,dtype=np.int_)

r = np.array([[0 for l in range(3)] for k in range(num_atoms_h)]).astype(float)


def read_hist(counter2):
    m=0
    c=0
    
    multiple1 = (counter2*iter_num_) + 3
    multiple2 = (counter2*iter_num_) + 4
    multiple3 = (counter2*iter_num_) + 5
    multiple4 = (counter2*iter_num_) + 6
    
    my_hist = open("./output_files/my_history_MOI4", 'r')
    for j, line in enumerate(my_hist):
        if j == multiple1:
            x_dimension = float(line.split()[(0)])
        if  j == multiple2:
            y_dimension = float(line.split()[(1)])
        if j == multiple3:
            z_dimension = float(line.split()[(2)])
            xyz_dimensions = np.array([x_dimension,y_dimension,z_dimension])
        
        if j >= multiple4 and j % 2 ==1:
            r[c,m] = float(line.split()[0])
            m = m+1
            r[c,m] = float(line.split()[1])
            m = m+1
            r[c,m] = float(line.split()[2])
            m = 0
            c =c+1
        
        if c == num_atoms_h:
            
            return r, xyz_dimensions
            break
    
start_time = time.time()

for k in range(num_frames_h):
#for k in tqdm(range(num_frames_h)):
    r, xyz_dimes = read_hist(k)                                 #Collect xyz postions from each frame of my_trajectory.
    
    rin1     = r / xyz_min                                      # Scale the xyz coordinates to be in units of box=1.
    
    
    #print(all_dimes, dr, num_dr)
    
    r_dv     = rin1[:,np.newaxis,:] - rin1[np.newaxis,:,:]      # Minus each xyz frome every other one and store each 
                                                                #itter in 3d array. - called Distance Vectors
    r_dv     = r_dv - np.rint(r_dv)                             # Apply periodic boundaries ~ the above operation can
                                                                #result in the distance vectors outside of the box.
                                                                #np.rint() rounds the values to the closest integer
    r_sep    = np.sqrt(np.sum(r_dv**2,axis=-1))                 # Separation distances as sqrt. sum of x**2 + y**2 + z**2
                                                                #in a 2d array of size num_atoms_h * num_atoms_h as each row 
                                                                #of r_dv is considered. Note, half of the values will be
                                                                #repeated because it is a magnitued calculation. ~ theory ~
                                                                #sqrt.(xyz_1**2 - xyz_2**2) == sqrt.(xyz_1**2 - xyz_1**2)
                                                                #But, this calculation is square - so xyz1 is minused from xyz1,
                                                                #Hence, triu_indices_from() is required to extract the values
    r_sep    = r_sep[np.triu_indices_from(r_sep,k=1)]           #                                                      
    #print(len(r_sep))
    hist,edges = np.histogram(r_sep,bins=num_dr,range=(0.0,r_max))
    my_hist = my_hist + 2*hist
    #print(all_dimes,k)



print("Total runtime:\n--- %s seconds ---\n" % (time.time() - start_time))

rho  = float(num_atoms_h)                                       # Our calculation is done in box=1 units so rho (density) is n/1 = n
h_id = ( 4.0 * np.pi * rho / 3.0) * ( edges[1:num_dr+1]**3 - edges[0:num_dr]**3 ) # Calculation to determine equivalent average number
                                                                          # of atoms in the same interval for an ideal gas.
#(edges[1:nk+1]**3- edges[0:nk]**3))                                      #This section, in code, makes an array where:
                                                                          #element = (n + 1) - n
g    = my_hist / h_id / (num_atoms_h*num_frames_h)                  # Average number

print('Output to pair_distrib_MOI4.out' )
edges = edges*xyz_min                                               # Convert bin edges back to sigma=1 units
r_mid = 0.5*(edges[0:num_dr]+edges[1:num_dr+1])                     # Mid points of bins
np.savetxt('./output_files/pair_distrib_MOI4.out', np.c_[r_mid,g],fmt="%15.8f")





