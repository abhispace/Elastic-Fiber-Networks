#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May  3 03:42:02 2025

@author: abhinav
"""
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 30 22:49:54 2025
this code reads in all the simulation data for srand 113....122.
does not work for srand = 112 yet

it reads all the num_center, srand and diff_network values of dfar and energies in this order.
Then it finds the mean, median and, std and sem and writes them to a file.

when we use model 1, it also write the mean radial displacement of bending dominated networks to a master file for a specific N_d.

@author: abhinav
"""

import numpy as np
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
from matplotlib import cm
import matplotlib.pyplot as plt
import math

def angle(x1,y1,x2,y2):
    dx = x2-x1
    dy = y2-y1
    return(math.atan2(dy,dx))

L = 16
L = 64
# L = 128
num_pts = L*L

ranseed = '667720,601210'    # -- network 1

num_center_list = [1, 2, 5, 10, 15, 20, 25, 30, 35, 40]
srand_list = srand_list = [113, 114, 116, 117, 118, 119, 120, 121, 122]                         # 115 is BAD!!!!!

# num_dip = num_center*6
num = 5+1                                                                        # the total number of force steps
num = 10+1

pbond = 1

mu = 1
mu_c = 1
tol = 2e-6
tol = 1.0e-6
tol = 1.0e-7

kappa = 1e-6

if kappa == 1e-6:
    kappa_fname = 'kappa2_e-6/'
elif kappa == 1e-5:
    kappa_fname = 'kappa2_e-5/'
elif kappa == 5e-6:
    kappa_fname = 'kappa2_5e-6/'
elif kappa == 5e-5:
    kappa_fname = 'kappa2_5e-5/'
elif kappa == 1e-4:
    kappa_fname = 'kappa2_e-4/'

rlen = 0.9
rlen_txt = "%.4f" % rlen
rlen_txt_ring = "%.9f" % 0.1

pbond_string = "%.2f" % pbond # to write the filename
tol_str = "%.2e" % tol # to write the filename
kappa_str = "%.2e" % kappa # to write the filename
if kappa == 0:
    kappa_str = "0.00e+00" # to write the filename
    
mu_str = "%.4f" % mu # to write the filename
mu_c_str = "%.4f" % mu_c # to write the filename

base = "/home/abhinav/david/"


bndry_force = []
loc_dip_moment = []
far_dip_moment = []
ratio_dip_moment = []
radial_disp = []
std_radial_disp = []

for i in range(0,len(num_center_list)):
    num_center = num_center_list[i]
    num_dip = num_center*6
    
    bndry_force.append([])
    loc_dip_moment.append([])
    far_dip_moment.append([])
    ratio_dip_moment.append([])
    
    
    radial_disp.append([])
    std_radial_disp.append([])
    
    for j in range(0,len(srand_list)):
        srand = srand_list[j]
        
        bndry_force[i].append([])
        loc_dip_moment[i].append([])
        far_dip_moment[i].append([])
        ratio_dip_moment[i].append([])
        
        radial_disp[i].append([])
        std_radial_disp[i].append([])
            
        folder = 'lattice_nopbd_bndry_all_clamp_restlength_circle_inner_fixed_radial_arp/'     # this is for hex bonds enforeced
                                
        # reading displacement values binned per ring in every simulation
        disp_fname = base+folder+"txt/displacement/"+"bndry_node_radial_disp_mean_"+str(srand)+"_"+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+".txt"
        # heading = 'bin mid       Mean Radial dist        Meadian Radial dist        Std Radial dist'
        disp_data = np.loadtxt(disp_fname)
        radial_disp[i][j].append(disp_data[:,1])    # this is the mean displacement of all the nodes in a given ring
        if i == 0:
            radial_bins = disp_data[:,0]
            
radial_disp = np.squeeze(np.array(radial_disp))

# radial displacement
mean_radial_disp = np.mean(radial_disp, axis = 1)
median_radial_disp = np.median(radial_disp, axis = 1)
std_radial_disp = np.std(radial_disp, axis = 1)
sem_radial_disp = std_radial_disp/np.sqrt(np.shape(radial_disp)[0]*np.shape(radial_disp)[1])
        
#============================================================================================
# writing the radial displacements to an output file
#============================================================================================
for i in range(0,len(num_center_list)):
    num_center = num_center_list[i]
    outfname = base + folder + "txt/displacement/radial_disp_srand_113_122_N_"+str(num_center)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+".txt"
    heading = 'radius         mean disp         median disp       std disp         sem disp'
    fmt = '%15.7e', '%15.7e', '%15.7e', '%15.7e', '%15.7e'
    np.savetxt(outfname, np.column_stack((radial_bins, mean_radial_disp[i], median_radial_disp[i], std_radial_disp[i], sem_radial_disp[i])), header = heading, fmt = fmt)
