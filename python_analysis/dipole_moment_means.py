#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 12 23:18:06 2024
it reads the mean values (that have been written to file using dipole_moment_plots.py)
for each kappa and N value and produces some master plots
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
    
        fname = base+folder+"mean_dip_mom_bndry_force_"+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+".txt"    
        data = np.loadtxt(fname)        
        data_all[i].append(data)
        
    i += 1

data_all = np.array(data_all)

# for k = 1e-5, N = 20, 25:
# mean_data_e5 = np.array([5.8117768e-03, 9.4174380e-03])
# std_data_e5 = np.array([1.8043998e-03, 3.9185328e-03])
# data_e5_all = np.concatenate([data_all[1,:,1], mean_data_e5])
# data_std_e5_all = np.concatenate([data_all[1,:,2], std_data_e5])
# num_center_e5 = [5,10,15,20,25]

# standard error of mean = standard deviation / sqrt(n)
yerr = data_all[:,:,2]/np.sqrt(data_all[:,:,0])


# plotting the dipole moments vs number of dipoles for all kappa
plt.figure(figsize=(8, 5), dpi=300)

# plt.errorbar(num_center_list, data_all[0,:,1], yerr = data_all[0,:,2], c = 'b', fmt = 'o', label = '$\widetilde \kappa = 10^{-6}$')
# # for e-5::
# # plt.errorbar(num_center_e5, data_e5_all, yerr = data_std_e5_all, c = 'g', fmt = 'd', label = '$\widetilde \kappa = 10^{-5}$')
# plt.errorbar(num_center_list, data_all[1,:,1], yerr = data_all[1,:,2], c = 'g', fmt = 'd', label = '$\widetilde \kappa = 10^{-5}$')
# plt.errorbar(num_center_list, data_all[2,:,1], yerr = data_all[2,:,2], c = 'k', fmt = 's', label = '$\widetilde \kappa = 10^{-4}$')

plt.errorbar(num_center_list, data_all[0,:,1], yerr = yerr[0,:], c = 'b', fmt = 'o', label = '$\widetilde \kappa = 10^{-6}$')


x_linspace = np.linspace(num_center_list[0],num_center_list[-1],100) 
# x_linspace_e5 = np.linspace(num_center_e5[0],num_center_e5[-1],100) 
line_fit1 = data_all[0,0,1]*x_linspace/num_center_list[0]
plt.plot(x_linspace, line_fit1, c = 'b') 


if cluster_radial_flag == 1:
    plt.errorbar(num_center_list, data_all[1,:,1], yerr = yerr[1,:], c = 'g', fmt = 'd', label = '$\widetilde \kappa = 10^{-5}$')
    plt.errorbar(num_center_list, data_all[2,:,1], yerr = yerr[2,:], c = 'k', fmt = 's', label = '$\widetilde \kappa = 10^{-4}$')

    # for e-5
    line_fit2 = data_all[1,0,1]*x_linspace/num_center_list[0]
    # line_fit2 = data_e5_all[0]*x_linspace_e5/num_center_e5[0]
    line_fit3 = data_all[2,0,1]*x_linspace/num_center_list[0]
    plt.plot(x_linspace, line_fit2, c = 'g') 
    # plt.plot(x_linspace_e5, line_fit2, c = 'g') 
    plt.plot(x_linspace, line_fit3, c = 'k') 

plt.xlabel("Number of Dipoles, $N_d$", fontsize = 20)
plt.ylabel("$<D_{far}>$", fontsize = 20)
plt.xticks(fontsize = 20)
plt.yticks(fontsize = 20)
ax = plt.gca()
# Increase tick width and length for major and minor ticks
ax.tick_params(axis='both', which='major', width=3, length=10)  # Major ticks
ax.tick_params(axis='both', which='minor', width=2, length=5)   # Minor ticks
plt.yscale('log')
plt.xscale('log')

plt.subplots_adjust(right=0.96,left=0.16,top=0.98,bottom=0.16)

plt.legend(loc = 'best', fontsize = 15)
plt.savefig(base+plot_folder+"dip_mom_vs_N_"+pbond_string+"_"+tol_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+".png")

plt.clf()
plt.close()

####################################################
# same as above: but xaxis is now packing fraction
####################################################
num_center_list = np.array(num_center_list)
r_in = 12
packing_frac2 = num_center_list*( 6*(np.sqrt(3)/4) ) / (np.pi * r_in * r_in)  # this treats the area of a dipole as a circluar region. not perfect.

# the following is rigrous for the geometry we have
# unit_traingle_area = (1/2)*(1)*(np.sqrt(3)/2)   # area of a triangle
# hex_area = 6*unit_traingle_area       # one hex is made of 6 traingles
# area_top_cut = np.sqrt(3)/16.
# area_bottom_cut = 3*np.sqrt(3)/16.
# unit_dipole_area = hex_area + 12*area_top_cut + 6*area_bottom_cut       # it includes some cuts of outside traingles as well
unit_dipole_area = 6 * (np.sqrt(3)/4) * (1.5*1.5)     # 6 * area of unit traingle that makes a unit hexagon (note the side lenght  = 1.5 and not 1)
area_of_dipoles = num_center_list*unit_dipole_area
area_of_inner_region = 509.22   # this is gotten by 6* area of a triangle of side = 14. it is approximate total inner area

packing_frac3 = area_of_dipoles/area_of_inner_region

plt.figure(figsize=(8, 5), dpi=300)

# plt.errorbar(num_center_list, data_all[0,:,1], yerr = data_all[0,:,2], c = 'b', fmt = 'o', label = '$\widetilde \kappa = 10^{-6}$')
# # for e-5::
# # plt.errorbar(num_center_e5, data_e5_all, yerr = data_std_e5_all, c = 'g', fmt = 'd', label = '$\widetilde \kappa = 10^{-5}$')
# plt.errorbar(num_center_list, data_all[1,:,1], yerr = data_all[1,:,2], c = 'g', fmt = 'd', label = '$\widetilde \kappa = 10^{-5}$')
# plt.errorbar(num_center_list, data_all[2,:,1], yerr = data_all[2,:,2], c = 'k', fmt = 's', label = '$\widetilde \kappa = 10^{-4}$')

# plt.errorbar(num_center_list, data_all[0,:,1], yerr = yerr[0,:], c = 'b', fmt = 'o', label = '$\widetilde \kappa = 10^{-6}$')
# plt.errorbar(num_center_list, data_all[1,:,1], yerr = yerr[1,:], c = 'g', fmt = 'd', label = '$\widetilde \kappa = 10^{-5}$')
# plt.errorbar(num_center_list, data_all[2,:,1], yerr = yerr[2,:], c = 'k', fmt = 's', label = '$\widetilde \kappa = 10^{-4}$')

plt.errorbar(packing_frac2, data_all[0,:,1], yerr = yerr[0,:], c = 'b', fmt = 'o', label = '$\widetilde \kappa = 10^{-6}$')
plt.errorbar(packing_frac2, data_all[1,:,1], yerr = yerr[1,:], c = 'g', fmt = 'd', label = '$\widetilde \kappa = 10^{-5}$')
plt.errorbar(packing_frac2, data_all[2,:,1], yerr = yerr[2,:], c = 'k', fmt = 's', label = '$\widetilde \kappa = 10^{-4}$')

# x_linspace = np.linspace(num_center_list[0],num_center_list[-1],100) 
x_linspace = np.linspace(packing_frac2[0],packing_frac2[-1],100) 
# x_linspace_e5 = np.linspace(num_center_e5[0],num_center_e5[-1],100) 
# line_fit1 = data_all[0,0,1]*x_linspace/num_center_list[0]
line_fit1 = data_all[0,0,1]*x_linspace/packing_frac2[0]
# for e-5
# line_fit2 = data_all[1,0,1]*x_linspace/num_center_list[0]
line_fit2 = data_all[1,0,1]*x_linspace/packing_frac2[0]
# line_fit2 = data_e5_all[0]*x_linspace_e5/num_center_e5[0]
# line_fit3 = data_all[2,0,1]*x_linspace/num_center_list[0]
line_fit3 = data_all[2,0,1]*x_linspace/packing_frac2[0]

plt.plot(x_linspace, line_fit1, c = 'b') 
plt.plot(x_linspace, line_fit2, c = 'g') 
# plt.plot(x_linspace_e5, line_fit2, c = 'g') 

plt.plot(x_linspace, line_fit3, c = 'k') 


plt.xlabel("Packing Fraction", fontsize = 15)
plt.ylabel("$<D_{far}>$", fontsize = 15)
plt.xticks([5, 10, 15], fontsize = 15)
plt.yticks(fontsize = 15)
plt.yscale('log')
plt.xscale('log')

plt.subplots_adjust(right=0.96,left=0.18,top=0.94,bottom=0.12)

plt.legend(loc = 'best')
plt.savefig(base+plot_folder+"dip_mom_vs_N_packing_frac_"+pbond_string+"_"+tol_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+".png")

plt.clf()
plt.close()


# plotting the LOCAL dipole moments vs number of dipoles for all kappa
plt.figure(figsize=(8, 5), dpi=300)

# plt.errorbar(num_center_list, data_all[2,:,5], yerr = data_all[2,:,6], c = 'k', fmt = 's', alpha = 0.5, label = '$\widetilde \kappa = 10^{-4}$')
# plt.errorbar(num_center_list, data_all[1,:,5], yerr = data_all[1,:,6], c = 'g', fmt = 'd', alpha = 0.5, label = '$\widetilde \kappa = 10^{-5}$')
# plt.errorbar(num_center_list, data_all[0,:,5], yerr = data_all[0,:,6], c = 'b', fmt = 'o', alpha = 0.5, label = '$\widetilde \kappa = 10^{-6}$')

plt.errorbar(num_center_list, data_all[2,:,5], yerr = yerr[2,:], ms = 15, c = 'k', fmt = 's', alpha = 0.5, label = '$\widetilde \kappa = 10^{-4}$')
plt.errorbar(num_center_list, data_all[1,:,5], yerr = yerr[1,:], ms = 15, c = 'g', fmt = 'd', alpha = 0.7, label = '$\widetilde \kappa = 10^{-5}$')
plt.errorbar(num_center_list, data_all[0,:,5], yerr = yerr[0,:], ms = 15, c = 'r', fmt = 'o', alpha = 0.5, label = '$\widetilde \kappa = 10^{-6}$')

x_linspace = np.linspace(num_center_list[0],num_center_list[-1],100) 
line_fit1 = data_all[0,0,5]*x_linspace/num_center_list[0]
# for e-5
line_fit2 = data_all[1,0,5]*x_linspace/num_center_list[0]
# line_fit2 = data_e5_all[0]*x_linspace_e5/num_center_e5[0]
line_fit3 = data_all[2,0,5]*x_linspace/num_center_list[0]

plt.plot(x_linspace, line_fit1, c = 'b', lw = 5, label = '$ {N}^{1} $') 
# plt.plot(x_linspace, line_fit2, c = 'g') 
# plt.plot(x_linspace_e5, line_fit2, c = 'g') 
# plt.plot(x_linspace, line_fit3, c = 'k') 

plt.xlabel("Number of Dipoles, $N_d$", fontsize = 30)
plt.ylabel("$<D_{loc}>$", fontsize = 30)
plt.xticks(num_center_list, fontsize = 30)
plt.yticks(fontsize = 30)
plt.yscale('log')
plt.xscale('log')

plt.subplots_adjust(right=0.98,left=0.17,top=0.98,bottom=0.20)

plt.legend(loc = 'best', fontsize = 24)
plt.savefig(base+plot_folder+"dip_loc_mom_vs_N_"+pbond_string+"_"+tol_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+".png")

plt.clf()
plt.close()

# plotting the dipole moments vs number of dipoles for all kappa - on a linear scale
plt.figure(figsize=(8, 5), dpi=300)

# plt.errorbar(num_center_list, data_all[0,:,1], yerr = data_all[0,:,2], c = 'b', fmt = 'o', label = '$\widetilde \kappa = 10^{-6}$')
# # for e-5::
# # plt.errorbar(num_center_e5, data_e5_all, yerr = data_std_e5_all, c = 'g', fmt = 'd', label = '$\widetilde \kappa = 10^{-5}$')
# plt.errorbar(num_center_list, data_all[1,:,1], yerr = data_all[1,:,2], c = 'g', fmt = 'd', label = '$\widetilde \kappa = 10^{-5}$')
# plt.errorbar(num_center_list, data_all[2,:,1], yerr = data_all[2,:,2], c = 'k', fmt = 's', label = '$\widetilde \kappa = 10^{-4}$')

plt.errorbar(num_center_list, data_all[0,:,1], yerr = yerr[0,:], c = 'b', fmt = 'o', label = '$\widetilde \kappa = 10^{-6}$')
plt.errorbar(num_center_list, data_all[1,:,1], yerr = yerr[1,:], c = 'g', fmt = 'd', label = '$\widetilde \kappa = 10^{-5}$')
plt.errorbar(num_center_list, data_all[2,:,1], yerr = yerr[2,:], c = 'k', fmt = 's', label = '$\widetilde \kappa = 10^{-4}$')

# plt.plot(x_linspace, line_fit1, c = 'b') 
# plt.plot(x_linspace, line_fit2, c = 'g') 
# # plt.plot(x_linspace_e5, line_fit2, c = 'g') 
# plt.plot(x_linspace, line_fit3, c = 'k') 

plt.xlabel("Number of Dipoles, $N_d$", fontsize = 15)
plt.ylabel("$<D_{far}>$", fontsize = 15)
plt.xticks([5, 10, 15], fontsize = 15)
plt.yticks(fontsize = 15)
# plt.yscale('log')
# plt.xscale('log')

plt.subplots_adjust(right=0.96,left=0.18,top=0.94,bottom=0.12)

plt.legend(loc = 'best')
plt.savefig(base+plot_folder+"dip_mom_vs_N_linear_"+pbond_string+"_"+tol_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+".png")

plt.clf()
plt.close()

# plotting the dipole moments vs number of dipoles for all kappa - on a semi-log scale
plt.figure(figsize=(8, 5), dpi=300)
# plt.errorbar(num_center_list, data_all[0,:,1], yerr = data_all[0,:,2], c = 'b', fmt = 'o', label = '$\widetilde \kappa = 10^{-6}$')
# plt.errorbar(num_center_list, data_all[1,:,1], yerr = data_all[1,:,2], c = 'g', fmt = 'd', label = '$\widetilde \kappa = 10^{-5}$')
# plt.errorbar(num_center_list, data_all[2,:,1], yerr = data_all[2,:,2], c = 'k', fmt = 's', label = '$\widetilde \kappa = 10^{-4}$')

plt.errorbar(num_center_list, data_all[0,:,1], yerr = yerr[0,:], c = 'b', fmt = 'o', label = '$\widetilde \kappa = 10^{-6}$')
plt.errorbar(num_center_list, data_all[1,:,1], yerr = yerr[1,:], c = 'g', fmt = 'd', label = '$\widetilde \kappa = 10^{-5}$')
plt.errorbar(num_center_list, data_all[2,:,1], yerr = yerr[2,:], c = 'k', fmt = 's', label = '$\widetilde \kappa = 10^{-4}$')

line_fit_exp1 = 0.5*data_all[0,0,1]*pow(1.9,x_linspace/num_center_list[0])
line_fit_exp2 = 0.5*data_all[1,0,1]*pow(1.9,x_linspace/num_center_list[0])
line_fit_exp3 = 0.5*data_all[2,0,1]*pow(1.7,x_linspace/num_center_list[0])

# plt.plot(x_linspace, line_fit_exp1, c = 'b', label = '$1.9^{x}$') 
# plt.plot(x_linspace, line_fit_exp2, c = 'g', label = '$1.9^{x}$') 
# plt.plot(x_linspace, line_fit_exp3, c = 'k', label = '$1.7^{x}$') 
plt.xlabel("Number of Dipoles", fontsize = 15)
plt.ylabel("$<D_{far}>$", fontsize = 15)
plt.xticks(num_center_list, fontsize = 15)
plt.yticks(fontsize = 15)
plt.yscale('log')
# plt.xscale('log')
plt.subplots_adjust(right=0.96,left=0.18,top=0.94,bottom=0.12)
plt.legend(loc = 'best')
plt.savefig(base+plot_folder+"dip_mom_vs_N_semi_log_"+pbond_string+"_"+tol_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+".png")
plt.clf()
plt.close()

#============================================================================================

# plotting the dipole moments vs number of dipoles for all kappa - on a linear scale - plotting separately - first kappa 1e-6
plt.figure(figsize=(8, 5), dpi=300)
# plt.errorbar(num_center_list, data_all[0,:,1], yerr = data_all[0,:,2], c = 'b', fmt = 'o', label = '$\widetilde \kappa = 10^{-6}$')
plt.errorbar(num_center_list, data_all[0,:,1], yerr = yerr[0,:], c = 'b', fmt = 'o', label = '$\widetilde \kappa = 10^{-6}$')
plt.plot(x_linspace, line_fit1, c = 'b') 
plt.xlabel("Number of Dipoles", fontsize = 15)
plt.ylabel("$<D_{far}>$", fontsize = 15)
plt.xticks(num_center_list, fontsize = 15)
plt.yticks(fontsize = 15)
# plt.yscale('log')
# plt.xscale('log')
plt.subplots_adjust(right=0.96,left=0.18,top=0.94,bottom=0.12)
plt.legend(loc = 'best')
plt.savefig(base+plot_folder+"dip_mom_vs_N_linear_kappa_1e-6_"+pbond_string+"_"+tol_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+".png")
plt.clf()
plt.close()

# plotting the dipole moments vs number of dipoles for all kappa - on a linear scale - plotting separately - second kappa 1e-5
plt.figure(figsize=(8, 5), dpi=300)
# plt.errorbar(num_center_list, data_all[1,:,1], yerr = data_all[1,:,2], c = 'b', fmt = 'o', label = '$\widetilde \kappa = 10^{-5}$')
plt.errorbar(num_center_list, data_all[1,:,1], yerr = yerr[1,:], c = 'b', fmt = 'o', label = '$\widetilde \kappa = 10^{-5}$')
plt.plot(x_linspace, line_fit2, c = 'b') 
plt.xlabel("Number of Dipoles", fontsize = 15)
plt.ylabel("$<D_{far}>$", fontsize = 15)
plt.xticks(num_center_list, fontsize = 15)
plt.yticks(fontsize = 15)
# plt.yscale('log')
# plt.xscale('log')
plt.subplots_adjust(right=0.96,left=0.18,top=0.94,bottom=0.12)
plt.legend(loc = 'best')
plt.savefig(base+plot_folder+"dip_mom_vs_N_linear_kappa_1e-5_"+pbond_string+"_"+tol_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+".png")
plt.clf()
plt.close()

# plotting the dipole moments vs number of dipoles for all kappa - on a linear scale - plotting separately - third kappa 1e-4
plt.figure(figsize=(8, 5), dpi=300)
# plt.errorbar(num_center_list, data_all[2,:,1], yerr = data_all[2,:,2], c = 'b', fmt = 'o', label = '$\widetilde \kappa = 10^{-4}$')
plt.errorbar(num_center_list, data_all[2,:,1], yerr = yerr[2,:], c = 'b', fmt = 'o', label = '$\widetilde \kappa = 10^{-4}$')
plt.plot(x_linspace, line_fit3, c = 'b') 
plt.xlabel("Number of Dipoles", fontsize = 15)
plt.ylabel("$<D_{far}>$", fontsize = 15)
plt.xticks(num_center_list, fontsize = 15)
plt.yticks(fontsize = 15)
# plt.yscale('log')
# plt.xscale('log')
plt.subplots_adjust(right=0.96,left=0.18,top=0.94,bottom=0.12)
plt.legend(loc = 'best')
plt.savefig(base+plot_folder+"dip_mom_vs_N_linear_kappa_1e-4_"+pbond_string+"_"+tol_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+".png")
plt.clf()
plt.close()


# plotting the dipole moments vs number of dipoles for all kappa - on a log scale - plotting separately - first kappa 1e-6
plt.figure(figsize=(8, 5), dpi=300)
# plt.errorbar(num_center_list, data_all[0,:,1], yerr = data_all[0,:,2], c = 'b', fmt = 'o', label = '$\widetilde \kappa = 10^{-6}$')
plt.errorbar(num_center_list, data_all[0,:,1], yerr = yerr[0,:], c = 'b', fmt = 'o', label = '$\widetilde \kappa = 10^{-6}$')
plt.plot(x_linspace, line_fit1, c = 'b') 
plt.xlabel("Number of Dipoles", fontsize = 15)
plt.ylabel("$<D_{far}>$", fontsize = 15)
plt.xticks(num_center_list, fontsize = 15)
plt.yticks(fontsize = 15)
plt.yscale('log')
plt.xscale('log')
plt.subplots_adjust(right=0.96,left=0.18,top=0.94,bottom=0.12)
plt.legend(loc = 'best')
plt.savefig(base+plot_folder+"dip_mom_vs_N_log_kappa_1e-6_"+pbond_string+"_"+tol_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+".png")
plt.clf()
plt.close()

# plotting the dipole moments vs number of dipoles for all kappa - on a semi-log scale - plotting separately - first kappa 1e-6
plt.figure(figsize=(8, 5), dpi=300)
# plt.errorbar(num_center_list, data_all[0,:,1], yerr = data_all[0,:,2], c = 'b', fmt = 'o', label = '$\widetilde \kappa = 10^{-6}$')
plt.errorbar(num_center_list, data_all[0,:,1], yerr = yerr[0,:], c = 'b', fmt = 'o', label = '$\widetilde \kappa = 10^{-6}$')
line_fit_exp1 = 0.5*data_all[0,0,1]*pow(1.9,x_linspace/num_center_list[0])
plt.plot(x_linspace, line_fit_exp1, c = 'b') 
plt.xlabel("Number of Dipoles", fontsize = 15)
plt.ylabel("$<D_{far}>$", fontsize = 15)
plt.xticks(num_center_list, fontsize = 15)
plt.yticks(fontsize = 15)
plt.yscale('log')
# plt.xscale('log')
plt.subplots_adjust(right=0.96,left=0.18,top=0.94,bottom=0.12)
plt.legend(loc = 'best')
plt.savefig(base+plot_folder+"dip_mom_vs_N_semi_log_kappa_1e-6_"+pbond_string+"_"+tol_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+".png")
plt.clf()
plt.close()

# plotting the dipole moments vs number of dipoles for all kappa - on a log scale - plotting separately - second kappa 1e-5
plt.figure(figsize=(8, 5), dpi=300)
# plt.errorbar(num_center_list, data_all[1,:,1], yerr = data_all[1,:,2], c = 'b', fmt = 'o', label = '$\widetilde \kappa = 10^{-5}$')
plt.errorbar(num_center_list, data_all[1,:,1], yerr = yerr[1,:], c = 'b', fmt = 'o', label = '$\widetilde \kappa = 10^{-5}$')
plt.plot(x_linspace, line_fit2, c = 'b') 
plt.xlabel("Number of Dipoles", fontsize = 15)
plt.ylabel("$<D_{far}>$", fontsize = 15)
plt.xticks(num_center_list, fontsize = 15)
plt.yticks(fontsize = 15)
plt.yscale('log')
plt.xscale('log')
plt.subplots_adjust(right=0.96,left=0.18,top=0.94,bottom=0.12)
plt.legend(loc = 'best')
plt.savefig(base+plot_folder+"dip_mom_vs_N_log_kappa_1e-5_"+pbond_string+"_"+tol_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+".png")
plt.clf()
plt.close()

# plotting the dipole moments vs number of dipoles for all kappa - on a log scale - plotting separately - third kappa 1e-4
plt.figure(figsize=(8, 5), dpi=300)
# plt.errorbar(num_center_list, data_all[2,:,1], yerr = data_all[2,:,2], c = 'b', fmt = 'o', label = '$\widetilde \kappa = 10^{-4}$')
plt.errorbar(num_center_list, data_all[2,:,1], yerr = yerr[2,:], c = 'b', fmt = 'o', label = '$\widetilde \kappa = 10^{-4}$')
line_fit3_2 = 0.1*data_all[2,0,1]*pow(x_linspace/num_center_list[0],2.9)
plt.plot(x_linspace, line_fit3, c = 'b', label = 'x') 
plt.plot(x_linspace, line_fit3_2, c = 'r', label = '$x^{2.9}$') 
plt.xlabel("Number of Dipoles", fontsize = 15)
plt.ylabel("$<D_{far}>$", fontsize = 15)
plt.xticks(num_center_list, fontsize = 15)
plt.yticks(fontsize = 15)
plt.yscale('log')
plt.xscale('log')
plt.subplots_adjust(right=0.96,left=0.18,top=0.94,bottom=0.12)
plt.legend(loc = 'best')
plt.savefig(base+plot_folder+"dip_mom_vs_N_log_kappa_1e-4_"+pbond_string+"_"+tol_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+".png")
plt.clf()
plt.close()

#============================================================================================


# plotting the dipole moments vs kappa for all dipole moments
plt.figure(figsize=(8, 5), dpi=300)

# plt.errorbar(kappa_list, data_all[:,0,1], yerr = data_all[:,0,2], c = 'b', fmt = 'o', label = 'N = 1')
# plt.errorbar(kappa_list, data_all[:,1,1], yerr = data_all[:,1,2], c = 'g', fmt = 'd', label = 'N = 2')
# plt.errorbar(kappa_list, data_all[:,2,1], yerr = data_all[:,2,2], c = 'k', fmt = 's', label = 'N = 5')
# plt.errorbar(kappa_list, data_all[:,3,1], yerr = data_all[:,3,2], c = 'r', fmt = 'o', label = 'N = 10')
# plt.errorbar(kappa_list, data_all[:,4,1], yerr = data_all[:,4,2], c = 'y', fmt = 'd', label = 'N = 15')
# plt.errorbar(kappa_list, data_all[:,5,1], yerr = data_all[:,5,2], c = 'm', fmt = 's', label = 'N = 20')
# plt.errorbar(kappa_list, data_all[:,6,1], yerr = data_all[:,6,2], c = 'c', fmt = 'o', alpha = 0.3, label = 'N = 25')
# plt.errorbar(kappa_list, data_all[:,7,1], yerr = data_all[:,7,2], c = 'coral', fmt = 'd', alpha = 0.3, label = 'N = 30')
# plt.errorbar(kappa_list, data_all[:,8,1], yerr = data_all[:,8,2], c = 'orange', fmt = 'o', alpha = 0.3, label = 'N = 35')
# plt.errorbar(kappa_list, data_all[:,9,1], yerr = data_all[:,9,2], c = 'darkblue', fmt = 'd', alpha = 0.3, label = 'N = 40')

plt.errorbar(kappa_list, data_all[:,0,1], yerr = yerr[:,0], c = 'b', fmt = 'o', label = 'N = 1')
# plt.errorbar(kappa_list, data_all[:,1,1], yerr = yerr[:,1], c = 'g', fmt = 'd', label = 'N = 2')
# plt.errorbar(kappa_list, data_all[:,2,1], yerr = yerr[:,2], c = 'k', fmt = 's', label = 'N = 5')
# plt.errorbar(kappa_list, data_all[:,3,1], yerr = yerr[:,3], c = 'r', fmt = 'o', label = 'N = 10')
# plt.errorbar(kappa_list, data_all[:,4,1], yerr = yerr[:,4], c = 'y', fmt = 'd', label = 'N = 15')
# plt.errorbar(kappa_list, data_all[:,5,1], yerr = yerr[:,5], c = 'm', fmt = 's', label = 'N = 20')
# plt.errorbar(kappa_list, data_all[:,6,1], yerr = yerr[:,6], c = 'c', fmt = 'o', alpha = 0.3, label = 'N = 25')
# plt.errorbar(kappa_list, data_all[:,7,1], yerr = yerr[:,7], c = 'coral', fmt = 'd', alpha = 0.3, label = 'N = 30')
# plt.errorbar(kappa_list, data_all[:,8,1], yerr = yerr[:,8], c = 'orange', fmt = 'o', alpha = 0.3, label = 'N = 35')
# plt.errorbar(kappa_list, data_all[:,9,1], yerr = yerr[:,9], c = 'darkblue', fmt = 'd', alpha = 0.3, label = 'N = 40')
plt.errorbar(kappa_list, data_all[:,1,1], yerr = yerr[:,1], c = 'k', fmt = 's', label = 'N = 5')
plt.errorbar(kappa_list, data_all[:,2,1], yerr = yerr[:,2], c = 'r', fmt = 'o', label = 'N = 10')
plt.errorbar(kappa_list, data_all[:,3,1], yerr = yerr[:,3], c = 'y', fmt = 'd', label = 'N = 15')
plt.errorbar(kappa_list, data_all[:,4,1], yerr = yerr[:,4], c = 'm', fmt = 's', label = 'N = 20')

x_linspace = np.linspace(kappa_list[0],kappa_list[-1],100) 
line_fit1 = data_all[0,0,1]*x_linspace/kappa_list[0]
# line_fit2 = data_all[0,1,1]*x_linspace/kappa_list[0]
# line_fit3 = data_all[0,2,1]*x_linspace/kappa_list[0]
# line_fit4 = data_all[0,3,1]*x_linspace/kappa_list[0]
# line_fit5 = data_all[0,4,1]*x_linspace/kappa_list[0]
# line_fit6 = data_all[0,5,1]*x_linspace/kappa_list[0]
# line_fit7 = data_all[0,6,1]*x_linspace/kappa_list[0]
# line_fit8 = data_all[0,7,1]*x_linspace/kappa_list[0]
# line_fit9 = data_all[0,8,1]*x_linspace/kappa_list[0]
# line_fit10 = data_all[0,9,1]*x_linspace/kappa_list[0]
line_fit3 = data_all[0,1,1]*x_linspace/kappa_list[0]
line_fit4 = data_all[0,2,1]*x_linspace/kappa_list[0]
line_fit5 = data_all[0,3,1]*x_linspace/kappa_list[0]
line_fit6 = data_all[0,4,1]*x_linspace/kappa_list[0]

plt.plot(x_linspace, line_fit1, c = 'b') 
# plt.plot(x_linspace, line_fit2, c = 'g') 
plt.plot(x_linspace, line_fit3, c = 'k') 
plt.plot(x_linspace, line_fit4, c = 'r') 
plt.plot(x_linspace, line_fit5, c = 'y') 
plt.plot(x_linspace, line_fit6, c = 'm') 
# plt.plot(x_linspace, line_fit7, c = 'c') 
# plt.plot(x_linspace, line_fit8, c = 'coral') 
# plt.plot(x_linspace, line_fit9, c = 'orange') 
# plt.plot(x_linspace, line_fit10, c = 'darkblue') 

plt.xlabel("Bending Modulus", fontsize = 15)
plt.ylabel("$<D_{far}>$", fontsize = 15)
plt.xticks([1e-6, 1e-5, 1e-4], fontsize = 15)
plt.yticks(fontsize = 15)
plt.yscale('log')
plt.xscale('log')

plt.subplots_adjust(right=0.96,left=0.18,top=0.94,bottom=0.12)

plt.legend(loc = 'best', fontsize = 8)
plt.savefig(base+plot_folder+"dip_mom_vs_kappa_"+pbond_string+"_"+tol_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+".png")

plt.clf()
plt.close()


# plotting the dipole moments vs kappa for all dipole moments after taking logs
plt.figure(figsize=(8, 5), dpi=300)

# plt.errorbar(kappa_list, data_all[:,0,1], yerr = data_all[:,0,2], c = 'b', fmt = 'o', label = 'N = 5')
# plt.errorbar(kappa_list, data_all[:,1,1], yerr = data_all[:,1,2], c = 'g', fmt = 'd', label = 'N = 10')
# plt.errorbar(kappa_list, data_all[:,2,1], yerr = data_all[:,2,2], c = 'k', fmt = 's', label = 'N = 15')

plt.errorbar(kappa_list, data_all[:,0,1], yerr = yerr[:,0], c = 'b', fmt = 'o', label = 'N = 5')
plt.errorbar(kappa_list, data_all[:,1,1], yerr = yerr[:,1], c = 'g', fmt = 'd', label = 'N = 10')
plt.errorbar(kappa_list, data_all[:,2,1], yerr = yerr[:,2], c = 'k', fmt = 's', label = 'N = 15')

x_linspace = np.linspace(kappa_list[0],kappa_list[-1],100) 
pow_val = 0.95
line_fit1 = data_all[0,0,1]*pow(x_linspace/kappa_list[0],pow_val)   # scaling for fit done here
line_fit2 = data_all[0,1,1]*pow(x_linspace/kappa_list[0],pow_val)
line_fit3 = data_all[0,2,1]*pow(x_linspace/kappa_list[0],pow_val)

plt.plot(x_linspace, line_fit1, c = 'b', label = ('$ {\\frac{x}{x_{0}}}^{0.95} $')) 
plt.plot(x_linspace, line_fit2, c = 'g', label = ('$ {\\frac{x}{x_{0}}}^{0.95} $')) 
plt.plot(x_linspace, line_fit3, c = 'k', label = ('$ {\\frac{x}{x_{0}}}^{0.95} $')) 


plt.xlabel("Bending Modulus", fontsize = 15)
plt.ylabel("$<D_{far}>$", fontsize = 15)
plt.xticks([1e-6, 1e-5, 1e-4], fontsize = 15)
plt.yticks(fontsize = 15)
plt.yscale('log')
plt.xscale('log')

plt.subplots_adjust(right=0.96,left=0.18,top=0.94,bottom=0.12)

plt.legend(loc = 'best')
plt.savefig(base+plot_folder+"dip_mom_vs_kappa_log_fit_scaled_"+pbond_string+"_"+tol_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+".png")

plt.clf()
plt.close()


# plotting the dipole moments vs kappa for all dipole moments after taking logs but not scaling the fits by x0
plt.figure(figsize=(8, 5), dpi=300)

# plt.errorbar(kappa_list, data_all[:,0,1], yerr = data_all[:,0,2], c = 'b', fmt = 'o', label = 'N = 5')
# plt.errorbar(kappa_list, data_all[:,1,1], yerr = data_all[:,1,2], c = 'g', fmt = 'd', label = 'N = 10')
# plt.errorbar(kappa_list, data_all[:,2,1], yerr = data_all[:,2,2], c = 'k', fmt = 's', label = 'N = 15')

plt.errorbar(kappa_list, data_all[:,0,1], yerr = yerr[:,0], c = 'b', fmt = 'o', label = 'N = 5')
plt.errorbar(kappa_list, data_all[:,1,1], yerr = yerr[:,1], c = 'g', fmt = 'd', label = 'N = 10')
plt.errorbar(kappa_list, data_all[:,2,1], yerr = yerr[:,2], c = 'k', fmt = 's', label = 'N = 15')

x_linspace = np.linspace(kappa_list[0],kappa_list[-1],100) 
pow_val = 0.95
line_fit1 = 52*pow(x_linspace,pow_val)
line_fit2 = 90*pow(x_linspace,pow_val)
line_fit3 = 155*pow(x_linspace,pow_val)

plt.plot(x_linspace, line_fit1, c = 'b', label = ('$ {x}^{0.95} $')) 
plt.plot(x_linspace, line_fit2, c = 'g', label = ('$ {x}^{0.95} $')) 
plt.plot(x_linspace, line_fit3, c = 'k', label = ('$ {x}^{0.95} $')) 


plt.xlabel("Bending Modulus", fontsize = 15)
plt.ylabel("$<D_{far}>$", fontsize = 15)
plt.xticks([1e-6, 1e-5, 1e-4], fontsize = 15)
plt.yticks(fontsize = 15)
plt.yscale('log')
plt.xscale('log')

plt.subplots_adjust(right=0.96,left=0.18,top=0.94,bottom=0.12)

plt.legend(loc = 'best')
plt.savefig(base+plot_folder+"dip_mom_vs_kappa_log_fit_"+pbond_string+"_"+tol_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+".png")

plt.clf()
plt.close()

#------------------------------------------------------------------------------
# test of log-log plots and linearity
#------------------------------------------------------------------------------
a = np.linspace(0,100,100)
b = a
cval = a*2
d = a**2
e = np.exp(a)
f = a**3
g = a**4

# plotting the dipole moments vs kappa for all dipole moments after taking logs
plt.figure(figsize=(8, 5), dpi=300)

plt.plot(a, b, c = 'b', label = ('y = x')) 
plt.plot(a, cval, c = 'r', label = ('y = 2x')) 
plt.scatter(a, d, c = 'k', label = ('y = x^2')) 
# plt.scatter(a, e, c = 'y', label = ('y = e^x')) 
plt.scatter(a, f, c = 'y', label = ('y = x^3')) 
plt.scatter(a, g, c = 'm', label = ('y = x^4')) 


plt.xlabel("Dummy data x", fontsize = 15)
plt.ylabel("Dummy data y", fontsize = 15)
# plt.xticks([1e-6, 1e-5, 1e-4], fontsize = 15)
plt.yticks(fontsize = 15)
plt.yscale('log')
plt.xscale('log')

plt.subplots_adjust(right=0.96,left=0.18,top=0.94,bottom=0.12)

plt.legend(loc = 'best')
plt.savefig(base+plot_folder+"test_"+pbond_string+"_"+tol_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+".png")

plt.clf()
plt.close()




#***********************************************************************************************************************************************
# dipole moment slopes for p = 1
#***********************************************************************************************************************************************
# num_c = [1,2,3,4,5]#,10,20]
num_c = [1,2,3,4,5,10,15,20,25,30,35,40]
# num_c = [1,2,5,10]    # for 10 different dipole placements
num_c = np.array(num_c)

L = 64
# L = 128

if L == 128:
    num_c = np.array([1,5,20])

srand_flag = 1
# srand = 112
# srand = 113
# srand = 114

fname_base  = '/home/abhinav/david/lattice_nopbd_bndry_all_clamp_restlength_circle_inner_fixed_radial_arp/txt/area/'
fname_plot = '/home/abhinav/david/lattice_nopbd_bndry_all_clamp_restlength_circle_inner_fixed_radial_arp/png/strain_hist/'
# for radial only bonds
radial_only_flag = 0
if radial_only_flag == 1:
    fname_base  = '/home/abhinav/david/lattice_nopbd_bndry_all_clamp_restlength_circle_inner_fixed_only_radial_arp/txt/area/'
    fname_plot = '/home/abhinav/david/lattice_nopbd_bndry_all_clamp_restlength_circle_inner_fixed_only_radial_arp/png/strain_hist/'

if L != 64:
    fname_plot = fname_plot + str(L)+'_'
    

pbond_string = "%.2f" % 1 # to write the filename
kappa = 1e-6
kappa_str = "%.2e" % kappa # to write the filename

d_far_p1 = []
d_loc_p1 = []
d_ratio_p1 = []
for i in range(0,len(num_c)):
    fname_p1 = fname_base + 'dipole_moment_ratio_'+str(num_c[i])+'_'+str(num_c[i]*6)+'_'+pbond_string+ "_"+tol_str+"_"+ kappa_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+'_10.txt'
    if srand_flag == 1 and num_c[i] > 1 and L == 64:
        fname_p1 = fname_base + 'dipole_moment_ratio_'+ 'srand_'+ str(srand) +'_'+str(num_c[i])+'_'+str(num_c[i]*6)+'_'+pbond_string+ "_"+tol_str+"_"+ kappa_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+'_10.txt'
    elif srand_flag == 1 and L != 64:
        fname_p1 = fname_base + 'dipole_moment_ratio_'+ 'srand_'+ str(srand) +'_'+ str(L) +'_'+str(num_c[i])+'_'+str(num_c[i]*6)+'_'+pbond_string+ "_"+tol_str+"_"+ kappa_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+'_10.txt'
    elif srand_flag == 1 and radial_only_flag == 1:
        fname_p1 = fname_base + 'dipole_moment_ratio_'+ 'srand_'+ str(srand) +'_'+str(num_c[i])+'_'+str(num_c[i]*6)+'_'+pbond_string+ "_"+tol_str+"_"+ kappa_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+'_10.txt'

    print('Fname_p1: ',fname_p1)
    p1_data_temp = np.loadtxt(fname_p1)
    d_far_p1.append(p1_data_temp[2])
    d_loc_p1.append(p1_data_temp[1])
    d_ratio_p1.append(p1_data_temp[3])
    
# d_far_p1 = [6.0463340e-01, 1.2133648e+00, 1.8212706e+00, 2.4224860e+00, 3.0219202e+00]
# d_loc_p1 = [5.7661642e-01, 1.1533254e+00, 1.7300822e+00, 2.3069075e+00, 2.8838063e+00]
d_far_p1 = np.array(d_far_p1)
d_loc_p1 = np.array(d_loc_p1)
d_ratio_p1 = np.array(d_ratio_p1)

# d_loc_p1 = []
# d_loc_p1 = np.array(d_loc_p1)

#  dipole moment
plt.figure(figsize=(8, 5), dpi=300)
x_linspace = np.linspace(num_c[0],num_c[-1],100) 
line_fit = d_far_p1[0]*x_linspace/num_c[0]
plt.plot(x_linspace, line_fit, c = 'b') 
plt.scatter(num_c, d_far_p1, marker = '.', s = 500, alpha=0.5)#, label = pbond_string)# Network 1')    
plt.xlabel("Number of Dipoles", fontsize = 15)
plt.ylabel("$D_{far}$", fontsize = 15)
plt.yticks(fontsize = 15)
plt.subplots_adjust(right=0.96,left=0.18,top=0.94,bottom=0.12)
plt.legend(loc = 'lower right')
if srand_flag == 1 and L == 64:
    plt.savefig(fname_plot+"dip_moment_p1_scatter_"+ 'srand_'+ str(srand) +"_" +str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+kappa_str+"_"+tol_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+".png")
else:
    plt.savefig(fname_plot+"dip_moment_p1_scatter_"+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+kappa_str+"_"+tol_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+".png")

plt.clf()
plt.close()


#  local dipole moment
plt.figure(figsize=(8, 5), dpi=300)
x_linspace = np.linspace(num_c[0],num_c[-1],100) 
line_fit = d_loc_p1[0]*x_linspace/num_c[0]
plt.plot(x_linspace, line_fit, c = 'b') 
plt.scatter(num_c, d_loc_p1, marker = '.', s = 500, alpha=0.5)#, label = pbond_string)# Network 1')    
plt.xlabel("Number of Dipoles", fontsize = 15)
plt.ylabel("$D_{loc}$", fontsize = 15)
plt.yticks(fontsize = 15)
plt.subplots_adjust(right=0.96,left=0.18,top=0.94,bottom=0.12)
plt.legend(loc = 'lower right')
if srand_flag == 1 and L == 64:
    plt.savefig(fname_plot+"dip_loc_moment_p1_scatter_"+ 'srand_'+ str(srand) +"_" +str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+kappa_str+"_"+tol_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+".png")
else:
    plt.savefig(fname_plot+"dip_loc_moment_p1_scatter_"+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+kappa_str+"_"+tol_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+".png")

plt.clf()
plt.close()

#  dipole moment divided by local dipole moment
# loc_dip_moment = 6*1.0*(1.0-0.9)*num_c 
# d_ratio_p1 = d_far_p1/loc_dip_moment
# d_ratio_p1 = d_far_p1/d_loc_p1

plt.figure(figsize=(8, 5), dpi=300)
plt.scatter(num_c, d_ratio_p1, marker = '.', s = 500, alpha=0.5)#, label = pbond_string)# Network 1')    
#plt.ylim([1.6,1.7])
plt.xlabel("Number of Dipoles", fontsize = 15)
plt.ylabel("$D_{far}$/$D_{loc}$", fontsize = 15)
plt.yticks(fontsize = 15)
plt.subplots_adjust(right=0.96,left=0.18,top=0.94,bottom=0.12)
plt.legend(loc = 'best')
if srand_flag == 1 and L == 64:
    plt.savefig(fname_plot+"dip_moment_ratio_p1_scatter_"+ 'srand_'+ str(srand) +"_" +str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+kappa_str+"_"+tol_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+".png")
else:
    plt.savefig(fname_plot+"dip_moment_ratio_p1_scatter_"+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+kappa_str+"_"+tol_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+".png")

plt.clf()
plt.close()


line_val = np.array([d_loc_p1[0]*x for x in num_c])
plt.figure(figsize=(8, 5), dpi=300)
plt.scatter(num_c, d_loc_p1-line_val, marker = '.', s = 500, alpha=0.5, color = 'r')#, label = pbond_string)# Network 1')    
#plt.plot(num_c, line_val)
#plt.ylim([1.6,1.7])
plt.xlabel("Number of Dipoles", fontsize = 15)
plt.ylabel("$D_{loc} - D_{loc,th}$", fontsize = 15)
plt.yticks(fontsize = 15)
plt.subplots_adjust(right=0.96,left=0.18,top=0.94,bottom=0.12)
plt.legend(loc = 'lower right')
if srand_flag == 1 and L == 64 :
    plt.savefig(fname_plot+"dip_moment_loc_deviation_p1_scatter_"+ 'srand_'+ str(srand) +"_" +str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+kappa_str+"_"+tol_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+".png")
else:
    plt.savefig(fname_plot+"dip_moment_loc_deviation_p1_scatter_"+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+kappa_str+"_"+tol_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+".png")

plt.clf()
plt.close()

# deviation percentage
plt.figure(figsize=(8, 5), dpi=300)
plt.scatter(num_c, (d_loc_p1-line_val)/line_val, marker = '.', s = 500, alpha=0.5, color = 'r')#, label = pbond_string)# Network 1')    
#plt.plot(num_c, line_val)
#plt.ylim([1.6,1.7])
plt.xlabel("Number of Dipoles", fontsize = 15)
plt.ylabel("$(D_{loc} - D_{loc,th})$/$D_{loc,th}$", fontsize = 15)
plt.yticks(fontsize = 15)
plt.subplots_adjust(right=0.96,left=0.18,top=0.94,bottom=0.12)
plt.legend(loc = 'lower right')
if srand_flag == 1 and L == 64 :
    plt.savefig(fname_plot+"dip_moment_loc_deviation_ratio_p1_scatter_"+ 'srand_'+ str(srand) +"_" +str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+kappa_str+"_"+tol_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+".png")
else:
    plt.savefig(fname_plot+"dip_moment_loc_deviation_ratio_p1_scatter_"+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+kappa_str+"_"+tol_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+".png")

plt.clf()
plt.close()


# '''
#==============================================================================
# putting dfar for p < 1 and p = 1 cases on the same plot
#==============================================================================

norm_d_far_p1 = d_far_p1/d_far_p1[0]
norm_d_far_e6 = data_all[0,:,1]/data_all[0,0,1]
norm_d_far_e5 = data_all[1,:,1]/data_all[1,0,1]
norm_d_far_e4 = data_all[2,:,1]/data_all[2,0,1]

plt.figure(figsize=(8, 5), dpi=300)

# fitting p = 1 network
x_linspace = np.linspace(num_c[0],num_c[-1],100) 
line_fit = norm_d_far_p1[0]*x_linspace/num_c[0]
plt.plot(x_linspace, line_fit, c = 'b') 

# fitting p < 1 network
x_linspace_k = np.linspace(num_center_list[5],num_center_list[-1],100) 
line_fit_k = (norm_d_far_e4[5]/num_center_list[5])*np.power(x_linspace_k,3)*0.003
# plt.plot(x_linspace_k, line_fit_k, c = 'r') 

plt.scatter(num_c, norm_d_far_p1, marker = '.', s = 500, alpha=0.5, label = '$p = 1$')# Network 1')  

# plt.errorbar(num_center_list, data_all[0,:,1], yerr = yerr[0,:], c = 'b', fmt = 'o', label = '$\widetilde \kappa = 10^{-6}$')
# plt.errorbar(num_center_list, data_all[1,:,1], yerr = yerr[1,:], c = 'g', fmt = 'd', label = '$\widetilde \kappa = 10^{-5}$')
# plt.errorbar(num_center_list, data_all[2,:,1], yerr = yerr[2,:], c = 'k', fmt = 's', label = '$\widetilde \kappa = 10^{-4}$')

plt.errorbar(num_center_list, norm_d_far_e6, c = 'b', fmt = 'o', label = '$\widetilde \kappa = 10^{-6}$')
plt.errorbar(num_center_list, norm_d_far_e5, c = 'g', fmt = 'd', label = '$\widetilde \kappa = 10^{-5}$')
plt.errorbar(num_center_list, norm_d_far_e4, c = 'k', fmt = 's', label = '$\widetilde \kappa = 10^{-4}$')

line_fit1 = norm_d_far_e6[0]*x_linspace/num_c[0]
# for e-5
line_fit2 = norm_d_far_e5[0]*x_linspace/num_c[0]
line_fit3 = norm_d_far_e4[0]*x_linspace/num_c[0]

plt.plot(x_linspace, line_fit1, c = 'b') 
plt.plot(x_linspace, line_fit2, c = 'g') 
plt.plot(x_linspace, line_fit3, c = 'k') 

plt.text(6, 3, '$\sim N_d^1$', fontsize = 20)
# plt.text(19, 100, '$\sim N_d^3$', fontsize = 20)

plt.xscale('log')
plt.yscale('log')
plt.xlabel("Number of Dipoles, $N_d$", fontsize = 20)
plt.ylabel("$<D_{far}>/<D_{far}>_{N=1}$", fontsize = 20)
plt.yticks(fontsize = 20)
plt.xticks(fontsize = 20)
ax = plt.gca()
# Increase tick width and length for major and minor ticks
ax.tick_params(axis='both', which='major', width=3, length=10)  # Major ticks
ax.tick_params(axis='both', which='minor', width=2, length=5)   # Minor ticks

plt.subplots_adjust(right=0.98,left=0.16,top=0.96,bottom=0.16)
plt.legend(loc = 'best', fontsize = 15)
if srand_flag == 1 and L == 64:
    plt.savefig(fname_plot+"dip_moment_all_p_scatter_norm_"+ 'srand_'+ str(srand) +"_" +str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+kappa_str+"_"+tol_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+".png")
else:
    plt.savefig(fname_plot+"dip_moment_all_p_scatter_norm_"+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+kappa_str+"_"+tol_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+".png")

plt.clf()
plt.close()

#------------------------------------------------------------------------------
# packing fraction
num_center_list = np.array(num_center_list)
r_in = 12
packing_frac1 = num_c*( 6*(np.sqrt(3)/4) ) / (np.pi * r_in * r_in)
packing_frac2 = num_center_list*( 6*(np.sqrt(3)/4) ) / (np.pi * r_in * r_in)

norm_d_far_p1 = d_far_p1/d_far_p1[0]
norm_d_far_e6 = data_all[0,:,1]/data_all[0,0,1]
norm_d_far_e5 = data_all[1,:,1]/data_all[1,0,1]
norm_d_far_e4 = data_all[2,:,1]/data_all[2,0,1]


plt.figure(figsize=(8, 5), dpi=300)

# x_linspace = np.linspace(num_c[0],num_c[-1],100) 
x_linspace = np.linspace(packing_frac2[0],packing_frac2[-1],100) 
line_fit = norm_d_far_p1[0]*x_linspace/packing_frac2[0]
plt.plot(x_linspace, line_fit, c = 'b') 
# plt.scatter(num_c, norm_d_far_p1, marker = '.', s = 500, alpha=0.5, label = '$p = 1$')# Network 1')  
plt.scatter(packing_frac1, norm_d_far_p1, marker = '.', s = 500, alpha=0.5, label = '$p = 1$')# Network 1')  

# plt.errorbar(num_center_list, data_all[0,:,1], yerr = yerr[0,:], c = 'b', fmt = 'o', label = '$\widetilde \kappa = 10^{-6}$')
# plt.errorbar(num_center_list, data_all[1,:,1], yerr = yerr[1,:], c = 'g', fmt = 'd', label = '$\widetilde \kappa = 10^{-5}$')
# plt.errorbar(num_center_list, data_all[2,:,1], yerr = yerr[2,:], c = 'k', fmt = 's', label = '$\widetilde \kappa = 10^{-4}$')

# plt.errorbar(num_center_list, norm_d_far_e6, c = 'b', fmt = 'o', label = '$\widetilde \kappa = 10^{-6}$')
# plt.errorbar(num_center_list, norm_d_far_e5, c = 'g', fmt = 'd', label = '$\widetilde \kappa = 10^{-5}$')
# plt.errorbar(num_center_list, norm_d_far_e4, c = 'k', fmt = 's', label = '$\widetilde \kappa = 10^{-4}$')


plt.errorbar(packing_frac2, norm_d_far_e6, c = 'b', fmt = 'o', label = '$\widetilde \kappa = 10^{-6}$')
plt.errorbar(packing_frac2, norm_d_far_e5, c = 'g', fmt = 'd', label = '$\widetilde \kappa = 10^{-5}$')
plt.errorbar(packing_frac2, norm_d_far_e4, c = 'k', fmt = 's', label = '$\widetilde \kappa = 10^{-4}$')

line_fit1 = norm_d_far_e6[0]*x_linspace/num_c[0]
# for e-5
line_fit2 = norm_d_far_e5[0]*x_linspace/num_c[0]
line_fit3 = norm_d_far_e4[0]*x_linspace/num_c[0]

# plt.plot(x_linspace, line_fit1, c = 'b') 
# plt.plot(x_linspace, line_fit2, c = 'g') 
# plt.plot(x_linspace, line_fit3, c = 'k') 

plt.xscale('log')
plt.yscale('log')
# plt.xlim([5e-3,5e-1])
plt.xlabel("Packing Fraction", fontsize = 15)
plt.ylabel("$<D_{far}>/<D_{far}>_{N=1}$", fontsize = 15)
plt.xticks(fontsize = 15)
plt.yticks(fontsize = 15)
plt.subplots_adjust(right=0.96,left=0.18,top=0.94,bottom=0.12)
plt.legend(loc = 'best')
if srand_flag == 1 and L == 64:
    plt.savefig(fname_plot+"dip_moment_all_p_scatter_norm_packing_frac_"+ 'srand_'+ str(srand) +"_" +str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+kappa_str+"_"+tol_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+".png")
else:
    plt.savefig(fname_plot+"dip_moment_all_p_scatter_norm_packing_frac_"+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+kappa_str+"_"+tol_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+".png")

plt.clf()
plt.close()

# '''
#==============================================================================
# finding alpha_m using dfar for p = 1 and for p = 0.55 networks
# we can then manually copy paste this data to network_circular_caller.py to finish the analysis
#==============================================================================
# for kappa = 1e-6
# creating the mask to remove N = 3 and N= 4 cases for which we do not have p = 0.55 results
d_far_p1_for_alpha_m = []
for i in range(0,len(num_c)):
    if num_c[i] != 3 and num_c[i] != 4:
        d_far_p1_for_alpha_m.append(d_far_p1[i])
        
alpha_m_kappae_6 = data_all[0,:,1]/d_far_p1_for_alpha_m   # for kappa = 1e-6 only
yerr_alpha_m_kappae_6 = yerr[0,:]/d_far_p1_for_alpha_m    # for kappa = 1e-6 only
# plotting alpha_m vs N
plt.figure(figsize=(8, 5), dpi=300)
# plt.scatter(num_center_list, alpha_m_kappae_6, marker = '.', s = 500, alpha=0.5, color = 'r')#, label = pbond_string)# Network 1')    
plt.errorbar(num_center_list, alpha_m_kappae_6, yerr = yerr_alpha_m_kappae_6, marker = '.', ms = 50, alpha=0.5, color = 'r')#, label = pbond_string)# Network 1')    
#plt.plot(num_c, line_val)
# plt.ylim([0,0.002])
plt.xlabel("Number of Dipoles", fontsize = 15)
plt.ylabel("$ \\alpha _{m} $", fontsize = 15)
plt.yticks(fontsize = 15)
plt.xticks(fontsize = 15)
plt.subplots_adjust(right=0.96,left=0.18,top=0.94,bottom=0.12)
plt.legend(loc = 'lower right')
if srand_flag == 1 and L == 64 :
    plt.savefig(fname_plot+"alpha_m_vs_N_kappa_e-6_"+ 'srand_'+ str(srand) +"_" +str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+kappa_str+"_"+tol_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+".png")
else:
    plt.savefig(fname_plot+"alpha_m_vs_N_kappa_e-6_"+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+kappa_str+"_"+tol_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+".png")

plt.clf()
plt.close()

