#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  6 11:24:05 2025
It reads stress on each concentric ring.
finds the mean and std over all data and plots it vs radius for different N.

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
num_pts = L*L

base = "/home/abhinav/david/"

cluster_new_flag = 1   # for dipoles with only radial bonds
hex_flag = 0           # for dipoles with hexagonal bonds forced to be present
hex_rand_flag = 0      # for dipoles with hexagonal bonds randomly removed as per p value
cluster_radial_flag = 1   # for dipoles with just radial bonds

# num_center_list = [5, 10, 15]
num_center_list = [1, 2, 5, 10, 15, 20, 25, 30, 35, 40]
# num_center_list = [1, 5, 10, 15, 20]    # for paper plot of dfar vs kappa

num = 10+1

pbond = 1
pbond = 0.5
pbond = 0.55
# pbond = 0.6

#pbond = 0.61
#pbond = 0.63
#pbond = 0.65
#pbond = 0.7

mu = 1
mu_c = 1
tol = 1.0e-7

kappa_list = [1e-6, 1e-5, 1e-4]
if num_center_list[-1] == 20:
    kappa_list = [1e-6, 5e-6, 1e-5, 5e-5, 1e-4]
if hex_flag == 1:
    kappa_list = [1e-6]

rlen = 0.9
rlen_txt = "%.4f" % rlen
rlen_txt_ring = "%.9f" % 0.1

srand = 112

pbond_string = "%.2f" % pbond # to write the filename
tol_str = "%.2e" % tol # to write the filename
    
mu_str = "%.4f" % mu # to write the filename
mu_c_str = "%.4f" % mu_c # to write the filename


data_all = []
i = 0

for kappa in kappa_list:
    data_all.append([])
    kappa_str = "%.2e" % kappa # to write the filename
    for num_center in num_center_list:
        num_dip = num_center*6
        
        ranseed = '667720,601210'    # -- network 1
    
        if kappa == 2e-7:
            kappa_fname = 'kappa2_2e-7/'
        elif kappa == 5e-7:
            kappa_fname = 'kappa2_5e-7/'
        elif kappa == 1e-6:
            kappa_fname = 'kappa2_e-6/'
        elif kappa == 1e-5:
            kappa_fname = 'kappa2_e-5/'
        elif kappa == 2e-6:
            kappa_fname = 'kappa2_2e-6/'
        elif kappa == 5e-6:
            kappa_fname = 'kappa2_5e-6/'
        elif kappa == 2e-5:
            kappa_fname = 'kappa2_2e-5/'
        elif kappa == 5e-5:
            kappa_fname = 'kappa2_5e-5/'
        elif kappa == 1e-4:
            kappa_fname = 'kappa2_e-4/'
        elif kappa == 2e-4:
            kappa_fname = 'kappa2_2e-4/'
        elif kappa == 5e-4:
            kappa_fname = 'kappa2_5e-4/'
        elif kappa == 1e-3:
            kappa_fname = 'kappa2_e-3/'
        elif kappa == 2e-3:
            kappa_fname = 'kappa2_2e-3/'
        elif kappa == 5e-3:
            kappa_fname = 'kappa2_5e-3/'
        elif kappa == 1e-2:
            kappa_fname = 'kappa2_e-2/'
        elif kappa == 0:
            kappa_fname = 'kappa2_0/'
        
        folder = 'lattice_nopbd_bndry_all_clamp_restlength_circle_inner_fixed_radial_arp_cluster/'+ kappa_fname +ranseed+"/"+str(srand)+"/means/"
        
        if cluster_new_flag == 1 and hex_flag == 0 and hex_rand_flag == 0 and cluster_radial_flag == 0:
            folder = 'cluster_download/lattice_nopbd_bndry_all_clamp_restlength_circle_inner_fixed_radial_arp_bash/'+ kappa_fname +ranseed+"/"+str(srand)+"/means/"
            plot_folder = 'cluster_download/lattice_nopbd_bndry_all_clamp_restlength_circle_inner_fixed_radial_arp_bash/'
        elif cluster_new_flag == 1 and hex_flag == 1 and hex_rand_flag == 0 and cluster_radial_flag == 0:
            folder = 'cluster_download/lattice_nopbd_bndry_all_clamp_restlength_circle_inner_fixed_radial_arp_bash_hex/'+ kappa_fname +ranseed+"/"+str(srand)+"/means/"
            plot_folder = 'cluster_download/lattice_nopbd_bndry_all_clamp_restlength_circle_inner_fixed_radial_arp_bash_hex/'
        elif cluster_new_flag == 1 and hex_flag == 0 and hex_rand_flag == 1 and cluster_radial_flag == 0:
            folder = 'cluster_download/lattice_nopbd_bndry_all_clamp_restlength_circle_inner_fixed_radial_arp_bash_random_dipole_hex_bonds/'+ kappa_fname +ranseed+"/"+str(srand)+"/means/"
            plot_folder = 'cluster_download/lattice_nopbd_bndry_all_clamp_restlength_circle_inner_fixed_radial_arp_bash_random_dipole_hex_bonds/'
        elif cluster_new_flag == 1 and hex_flag == 0 and hex_rand_flag == 0 and cluster_radial_flag == 1:
            folder = 'cluster_download/lattice_nopbd_bndry_all_clamp_restlength_circle_inner_fixed_radial_arp_bash_radial/'+ kappa_fname +ranseed+"/"+str(srand)+"/means/"
            plot_folder = 'cluster_download/lattice_nopbd_bndry_all_clamp_restlength_circle_inner_fixed_radial_arp_bash_radial/'
    
        fname = base+folder+"mean_ring_stress_"+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+".txt"    
        data = np.loadtxt(fname)        
        data_all[i].append(data)
        
    i += 1

data_all = np.array(data_all)

radial_dist = data_all[0,0,:,0]
cutoff_index = np.where(radial_dist >= 13)[0][0]

# plotting mean ring stress vs Radius for each N
for i in np.arange(0,len(num_center_list)):
    plt.figure(figsize=(8, 5), dpi=300)

    plt.errorbar(radial_dist, data_all[0,i,:,1], yerr = data_all[0,i,:,2], c = 'b', fmt = 'o', label = '$\widetilde \kappa = 10^{-6}$')
    plt.errorbar(radial_dist, data_all[1,i,:,1], yerr = data_all[1,i,:,2], c = 'g', fmt = 'd', label = '$\widetilde \kappa = 10^{-5}$')
    plt.errorbar(radial_dist, data_all[2,i,:,1], yerr = data_all[2,i,:,2], c = 'k', fmt = 's', label = '$\widetilde \kappa = 10^{-4}$')

    plt.xlabel("Radial Distance", fontsize = 15)
    plt.ylabel("Radial Stress", fontsize = 15)
    plt.xticks([5, 10, 15], fontsize = 15)
    plt.yticks(fontsize = 15)
    plt.yscale('log')
    plt.xscale('log')

    plt.subplots_adjust(right=0.96,left=0.18,top=0.94,bottom=0.12)

    plt.legend(loc = 'best')
    plt.savefig(base+plot_folder+"mean_ring_stress_vs_r_"+pbond_string+"_Nd_"+str(num_center_list[i])+"_"+tol_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+".png")

    plt.clf()
    plt.close()


    # only outer annulus now
    plt.figure(figsize=(8, 5), dpi=300)

    plt.errorbar(radial_dist[cutoff_index:], data_all[0,i,cutoff_index:,1], yerr = data_all[0,i,cutoff_index:,2], c = 'b', fmt = 'o', label = '$\widetilde \kappa = 10^{-6}$')
    plt.errorbar(radial_dist[cutoff_index:], data_all[1,i,cutoff_index:,1], yerr = data_all[1,i,cutoff_index:,2], c = 'g', fmt = 'd', label = '$\widetilde \kappa = 10^{-5}$')
    plt.errorbar(radial_dist[cutoff_index:], data_all[2,i,cutoff_index:,1], yerr = data_all[2,i,cutoff_index:,2], c = 'k', fmt = 's', label = '$\widetilde \kappa = 10^{-4}$')

    plt.xlabel("Radial Distance", fontsize = 15)
    plt.ylabel("Radial Stress", fontsize = 15)
    plt.xticks([5, 10, 15], fontsize = 15)
    plt.yticks(fontsize = 15)
    plt.yscale('log')
    plt.xscale('log')

    plt.subplots_adjust(right=0.96,left=0.18,top=0.94,bottom=0.12)

    plt.legend(loc = 'best')
    plt.savefig(base+plot_folder+"mean_outer_ring_stress_vs_r_"+pbond_string+"_Nd_"+str(num_center_list[i])+"_"+tol_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+".png")

    plt.clf()
    plt.close()



    # only outer annulus now and also normalized
    plt.figure(figsize=(8, 5), dpi=300)

    plt.errorbar(radial_dist[cutoff_index:], data_all[0,i,cutoff_index:,1]/data_all[0,i,cutoff_index,1], yerr = data_all[0,i,cutoff_index:,2]/data_all[0,i,cutoff_index:,1], c = 'b', fmt = 'o', label = '$\widetilde \kappa = 10^{-6}$')
    plt.errorbar(radial_dist[cutoff_index:], data_all[1,i,cutoff_index:,1]/data_all[1,i,cutoff_index,1], yerr = data_all[1,i,cutoff_index:,2]/data_all[1,i,cutoff_index:,1], c = 'g', fmt = 'd', label = '$\widetilde \kappa = 10^{-5}$')
    plt.errorbar(radial_dist[cutoff_index:], data_all[2,i,cutoff_index:,1]/data_all[2,i,cutoff_index,1], yerr = data_all[2,i,cutoff_index:,2]/data_all[2,i,cutoff_index:,1], c = 'k', fmt = 's', label = '$\widetilde \kappa = 10^{-4}$')

    plt.xlabel("Radial Distance", fontsize = 15)
    plt.ylabel("Normalized Radial Stress", fontsize = 15)
    plt.xticks(fontsize = 15)
    plt.yticks(fontsize = 15)
    plt.xlim([radial_dist[cutoff_index],radial_dist[-2]+0.5])
    # plt.yscale('log')
    # plt.xscale('log')

    plt.subplots_adjust(right=0.96,left=0.18,top=0.94,bottom=0.12)

    plt.legend(loc = 'best')
    plt.savefig(base+plot_folder+"mean_outer_ring_stress_vs_r_normalized_"+pbond_string+"_Nd_"+str(num_center_list[i])+"_"+tol_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+".png")

    plt.clf()
    plt.close()


marker_styles = ['o', 's', 'D', '^', 'v', '<', '>', 'p', '*', 'X']
colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k']  # blue, green, red, cyan, magenta, yellow, black
# Generate random selections
import random

# defining a line of slope 1
xaxis = np.logspace(-6, -4, 100)
# y_fit_linear40 = 5*xaxis   # works great for N =40
y_fit_linear = 0.1*xaxis   # works great for N = 10
    
random.seed(42)
# plotting mean ring stress vs kappa for each N
for i in np.arange(0,len(num_center_list)):
    plt.figure(figsize=(8, 5), dpi=300)

    for j in np.arange(0,len(radial_dist)-1):
        label = '$R = %.1f$' % (radial_dist[j])
        random_marker = random.choice(marker_styles)
        random_color = random.choice(colors)
        plt.errorbar(kappa_list, data_all[:,i,j,1], yerr = data_all[:,i,j,2], c = random_color, fmt = 'o', marker = random_marker, label = label)

    plt.xlabel("Bending modulus", fontsize = 15)
    plt.ylabel("Radial Stress", fontsize = 15)
    plt.xticks(fontsize = 15)
    plt.yticks(fontsize = 15)
    plt.xlim([8e-7,5e-4])
    plt.yscale('log')
    plt.xscale('log')

    plt.subplots_adjust(right=0.96,left=0.18,top=0.94,bottom=0.12)

    plt.legend(loc = 'best')
    plt.savefig(base+plot_folder+"mean_ring_stress_vs_kappa_"+pbond_string+"_Nd_"+str(num_center_list[i])+"_"+tol_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+".png")

    plt.clf()
    plt.close()


    # only outer annulus now
    plt.figure(figsize=(8, 5), dpi=300)

    if num_center_list[i] == 1:
        y_fit_linear = 0.01*xaxis   # works great for N = 10
        # yloc_text = 
    elif num_center_list[i] == 5:
        y_fit_linear = 0.1*xaxis   # works great for N = 10
    elif num_center_list[i] == 10:
        y_fit_linear = 0.1*xaxis   # works great for N = 10
    elif num_center_list[i] == 20:
        y_fit_linear = 0.5*xaxis   # works great for N = 10
    elif num_center_list[i] == 35:
        y_fit_linear = 2*xaxis   # works great for N = 10
    elif num_center_list[i] == 40:
        y_fit_linear = 5*xaxis   # works great for N = 10

    for j in np.arange(int(radial_dist[cutoff_index]),len(radial_dist)-1):
        label = '$R = %.1f$' % (radial_dist[j])
        random_marker = random.choice(marker_styles)
        random_color = random.choice(colors)
        plt.errorbar(kappa_list, data_all[:,i,j,1], yerr = data_all[:,i,j,2], c = random_color, fmt = 'o', marker = random_marker, label = label)
    plt.plot(xaxis, y_fit_linear)
    plt.text(3e-5, 5e-4, "$\widetilde{\\kappa}$", fontsize = 12)

    plt.xlabel("Bending Modulus", fontsize = 15)
    plt.ylabel("Radial Stress", fontsize = 15)
    plt.xticks(fontsize = 15)
    plt.yticks(fontsize = 15)
    plt.yscale('log')
    plt.xscale('log')
    plt.xlim([8e-7,5e-4])

    plt.subplots_adjust(right=0.96,left=0.18,top=0.94,bottom=0.12)

    plt.legend(loc = 'best')
    plt.savefig(base+plot_folder+"mean_outer_ring_stress_vs_kappa_"+pbond_string+"_Nd_"+str(num_center_list[i])+"_"+tol_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+".png")

    plt.clf()
    plt.close()

