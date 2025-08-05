#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 16 06:33:57 2025
reads the output of bndry_froces_circular.py

The output has pairwise distance for each pair of dipoles; stored in ...../means/

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

dont_throw_flag = 0                   # set to 1 if you do not want to throw the simulations where there maybe a completely disconnected node

# L = 16
L = 64
# L = 128

if L == 64:
    unconnected_const = 9999
if L == 128:
    unconnected_const = -1

num_pts = L*L

num_center = 5
num_dip = num_center * 6

# num_center_list = [5]
srand_list = ['113','114','116','117','118','119','120','121','122']

num = 10+1

pbond = 1
pbond = 0.55

mu = 1
mu_c = 1

tol = 1.0e-7

kappa = 1e-6

rlen = 0.9
rlen_txt = "%.4f" % rlen

radial_flag = 0    # for inner and outer boundaries different; dipoles are hexagons
radial_inner_flag = 0   # for dipoles with just radial forces and inner boundary different from outer
inner_radial_arp_flag = 0       # for hex bonds
cluster_flag = 0
cluster_new_flag = 0

hex_flag = 0  # these are for cluster 100 simulations where the dipoles have all 12 bonds present - including the outer hex bonds
hex_rand_flag = 1
cluster_radial_flag = 0

if hex_flag == 1 or hex_rand_flag == 1 or cluster_radial_flag == 1:
    cluster_new_flag = 1

radial_only_flag = 0
    
srand_flag = 1
if srand_flag == 1:
    srand = 113

if kappa == 1e-6:
    kappa_fname = 'kappa2_e-6/'

pbond_string = "%.2f" % pbond # to write the filename
tol_str = "%.2e" % tol # to write the filename
kappa_str = "%.2e" % kappa # to write the filename
if kappa == 0:
    kappa_str = "0.00e+00" # to write the filename
    
mu_str = "%.4f" % mu # to write the filename
mu_c_str = "%.4f" % mu_c # to write the filename

ranseed = '667720,601210'    # -- network 1

base = "/home/abhinav/david/"

############################################################################################################################
# reading all the pairwise distances for each N case and plotting
############################################################################################################################
mean_dist_pairs_all = []
std_dist_pairs_all = []
rg_all = []
convex_hull = []

for i in range(0,len(srand_list)):
    srand = srand_list[i]
    
    if cluster_new_flag == 1 and hex_flag == 0 and hex_rand_flag == 0 and cluster_radial_flag == 0:
        cent_dist_folder = 'cluster_download/lattice_nopbd_bndry_all_clamp_restlength_circle_inner_fixed_radial_arp_bash/'+ kappa_fname +ranseed+"/"+str(srand)+"/"
        plot_folder = 'cluster_download/lattice_nopbd_bndry_all_clamp_restlength_circle_inner_fixed_radial_arp_bash/'+ kappa_fname +ranseed+"/plots/"
    elif cluster_new_flag == 1 and hex_flag == 1 and hex_rand_flag == 0 and cluster_radial_flag == 0:
        cent_dist_folder = 'cluster_download/lattice_nopbd_bndry_all_clamp_restlength_circle_inner_fixed_radial_arp_bash_hex/'+ kappa_fname +ranseed+"/"+str(srand)+"/"
        plot_folder = 'cluster_download/lattice_nopbd_bndry_all_clamp_restlength_circle_inner_fixed_radial_arp_bash_hex/'+ kappa_fname +ranseed+"/plots/"
    elif cluster_new_flag == 1 and hex_flag == 0 and hex_rand_flag == 1 and cluster_radial_flag == 0:
        cent_dist_folder = 'cluster_download/lattice_nopbd_bndry_all_clamp_restlength_circle_inner_fixed_radial_arp_bash_random_dipole_hex_bonds/'+ kappa_fname +ranseed+"/"+str(srand)+"/"
        plot_folder = 'cluster_download/lattice_nopbd_bndry_all_clamp_restlength_circle_inner_fixed_radial_arp_bash_random_dipole_hex_bonds/'+ kappa_fname +ranseed+"/plots/"
    elif cluster_new_flag == 1 and hex_flag == 0 and hex_rand_flag == 0 and cluster_radial_flag == 1:
        cent_dist_folder = 'cluster_download/lattice_nopbd_bndry_all_clamp_restlength_circle_inner_fixed_radial_arp_bash_radial/'+ kappa_fname +ranseed+"/"+str(srand)+"/"    
        plot_folder = 'cluster_download/lattice_nopbd_bndry_all_clamp_restlength_circle_inner_fixed_radial_arp_bash_radial/'+ kappa_fname +ranseed+"/plots/"

    center_dist_fname = base+cent_dist_folder+"means/"+"center_dist_"+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+".txt"    
    # print('pairwise distance file:  ', center_dist_fname)
    dist_pairs = np.loadtxt(center_dist_fname)
    if num_center == 2:
        mean_dist_pairs = np.mean(dist_pairs[-1])
        std_dist_pairs = 0
    else:
        mean_dist_pairs = np.mean(dist_pairs[:,-1])
        std_dist_pairs = np.std(dist_pairs[:,-1])
    mean_dist_pairs_all.append(mean_dist_pairs)
    std_dist_pairs_all.append(std_dist_pairs)    
    
    
    center_quant_fname = base+cent_dist_folder+"means/"+"centers_quant_"+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+".txt"    
    quant_data = np.loadtxt(center_quant_fname)
    rg_all.append(quant_data[0])
    convex_hull.append(quant_data[1])    


#============================================================================================
# plotting the data
#============================================================================================
# plot of mean pairwise distances
fig = plt.figure(figsize=(8, 6), dpi=300)
plt.errorbar(np.arange(0,len(srand_list)), mean_dist_pairs_all, yerr= std_dist_pairs_all, ls = "None", marker = 'o', ms = 10, c = 'b')
plt.xlabel('Random Network', fontsize = 20)
plt.ylabel('Mean paiwise distance', fontsize = 20)
plt.xticks(fontsize = 20)
plt.yticks(fontsize = 20)
plt.subplots_adjust(top=0.96, left = 0.17, bottom = 0.15, right = 0.98)
# plt.legend(loc='best', fontsize = 15)
plt.savefig(base + plot_folder + 'pairwise_dist_centers_srand_113_122_Nd_'+str(num_center)+"_"+pbond_string+"_"+mu_str+"_"+mu_c_str+"_"+kappa_str+"_"+tol_str+"_"+rlen_txt+"_"+str(num-1)+".png")  # just plotting the last step data

plt.close()
plt.clf()

# plot of radius of gyration
fig = plt.figure(figsize=(8, 6), dpi=300)
plt.scatter(np.arange(0,len(srand_list)), rg_all, marker = 'o', s = 80, c = 'b')
plt.xlabel('Random Network', fontsize = 20)
plt.ylabel('$R_g$', fontsize = 20)
plt.xticks(fontsize = 20)
plt.yticks(fontsize = 20)
plt.subplots_adjust(top=0.96, left = 0.17, bottom = 0.15, right = 0.98)
# plt.legend(loc='best', fontsize = 15)
plt.savefig(base + plot_folder + 'radius_gyration_centers_srand_113_122_Nd_'+str(num_center)+"_"+pbond_string+"_"+mu_str+"_"+mu_c_str+"_"+kappa_str+"_"+tol_str+"_"+rlen_txt+"_"+str(num-1)+".png")  # just plotting the last step data

plt.close()
plt.clf()

# plot of convex hull area
fig = plt.figure(figsize=(8, 6), dpi=300)
plt.scatter(np.arange(0,len(srand_list)), convex_hull, marker = 'o', s = 80, c = 'b')
plt.xlabel('Random Network', fontsize = 20)
plt.ylabel('Convex Hull Area', fontsize = 20)
plt.xticks(fontsize = 20)
plt.yticks(fontsize = 20)
plt.subplots_adjust(top=0.96, left = 0.17, bottom = 0.15, right = 0.98)
# plt.legend(loc='best', fontsize = 15)
plt.savefig(base + plot_folder + 'convex_hull_srand_113_122_Nd_'+str(num_center)+"_"+pbond_string+"_"+mu_str+"_"+mu_c_str+"_"+kappa_str+"_"+tol_str+"_"+rlen_txt+"_"+str(num-1)+".png")  # just plotting the last step data

plt.close()
plt.clf()
