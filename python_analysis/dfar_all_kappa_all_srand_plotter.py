#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  9 07:01:39 2025

it reads the output from dfar_all_kappa_srand.py for all kappa
and then makes a plot of dfar vs N and packing fraction for the paper
the result is figure 2 in the contraction paper.

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
    
cluster_new_flag = 1   # for results from the cluster
hex_flag = 0           # for dipoles with hexagonal bonds forced to be present
hex_rand_flag = 0      # for dipoles with hexagonal bonds randomly removed as per p value
cluster_radial_flag = 1   # for dipoles with just radial bonds

# num_center = 1
num_center_list = [1, 2, 5, 10, 15, 20, 25, 30, 35, 40]
# srand_list = [112,113, 114, 116, 117, 118, 119, 120, 121, 122]                         # 115 is BAD!!!!!
kappa_list = [1e-6,1e-5,1e-4]

# num_dip = num_center*6
num = 5+1                                                                        # the total number of force steps
num = 10+1

pbond = 1
# pbond = 0.5
pbond = 0.55
# pbond = 0.6


mu = 1
mu_c = 1
tol = 2e-6
tol = 1.0e-6
tol = 1.0e-7

# kappa = 1e-6
# kappa = 1e-5
# kappa = 1e-4


rlen = 0.9
rlen_txt = "%.4f" % rlen
rlen_txt_ring = "%.9f" % 0.1


pbond_string = "%.2f" % pbond # to write the filename
tol_str = "%.2e" % tol # to write the filename
    
mu_str = "%.4f" % mu # to write the filename
mu_c_str = "%.4f" % mu_c # to write the filename

base = "/home/abhinav/david/"

all_dfar = []
for i in range(0,len(kappa_list)):      
    kappa = kappa_list[i]
    kappa_str = "%.2e" % kappa # to write the filename
    if kappa == 1e-6:
        kappa_fname = 'kappa2_e-6/'
    elif kappa == 1e-5:
        kappa_fname = 'kappa2_e-5/'
    elif kappa == 1e-4:
        kappa_fname = 'kappa2_e-4/'

    folder = 'cluster_download/lattice_nopbd_bndry_all_clamp_restlength_circle_inner_fixed_radial_arp_bash_radial/'+ kappa_fname +ranseed+"/means/"    
    plot_folder = 'cluster_download/lattice_nopbd_bndry_all_clamp_restlength_circle_inner_fixed_radial_arp_bash_radial/'
            
    dfar_fname = base+folder+"mean_dip_mom_bend_dom_vs_N_values_srand_113_122_"+pbond_string+"_"+tol_str+"_"+kappa_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+".txt"        
    dfar_temp = np.loadtxt(dfar_fname)
    all_dfar.append(dfar_temp)
    
all_dfar = np.array(all_dfar)
mean_dfar_kappa = all_dfar[:,:,1]
sem_dfar_kappa = all_dfar[:,:,-1]


#========================================================================================================================
# plotting
#========================================================================================================================
# plotting the dipole moments vs number of dipoles for all kappa
plt.figure(figsize=(8, 5), dpi=300)

plt.errorbar(num_center_list, mean_dfar_kappa[0], yerr = sem_dfar_kappa[0], c = 'b', fmt = 'o', ms = 10, label = '$\widetilde \kappa = 10^{-6}$')
plt.errorbar(num_center_list, mean_dfar_kappa[1], yerr = sem_dfar_kappa[1], c = 'g', fmt = 'd', ms = 10, label = '$\widetilde \kappa = 10^{-5}$')
plt.errorbar(num_center_list, mean_dfar_kappa[2], yerr = sem_dfar_kappa[2], c = 'k', fmt = 's', ms = 10, label = '$\widetilde \kappa = 10^{-4}$')

x_linspace = np.linspace(num_center_list[0],num_center_list[-1],100) 
line_fit1 = mean_dfar_kappa[0,0]*x_linspace/num_center_list[0]
plt.plot(x_linspace, line_fit1, c = 'b') 
line_fit2 = mean_dfar_kappa[1,0]*x_linspace/num_center_list[0]
line_fit3 = mean_dfar_kappa[2,0]*x_linspace/num_center_list[0]
plt.plot(x_linspace, line_fit2, c = 'g') 
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
plt.savefig(base+plot_folder+"bend_dom_112_through_122_dip_mom_vs_N_"+pbond_string+"_"+tol_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+".png")

plt.clf()
plt.close()

#==========================================================================================================================
# Reading p = 1 L = 64 network to plot with L = 128 network
#==========================================================================================================================
fname_base  = '/home/abhinav/david/lattice_nopbd_bndry_all_clamp_restlength_circle_inner_fixed_radial_arp/txt/area/'

srand_list_64 = ['113', '114','116','117','118','119','120','121','122']                         # 115 is BAD!!!!!
num_center_list_64 = np.array([1, 2, 5, 10, 15, 20, 25, 30, 35, 40])

kappa = 1e-6
kappa_str = "%.2e" % kappa # to write the filename

d_far_p1_64 = []
d_loc_p1_64 = []
for i in range(0,len(num_center_list_64)):
    d_far_p1_64.append([])
    d_loc_p1_64.append([])
    for j in range(0,len(srand_list_64)):    
        srand = srand_list_64[j]
        fname_p1 = fname_base + 'dipole_moment_ratio_'+ 'srand_'+ str(srand) +'_'+str(num_center_list_64[i])+'_'+str(num_center_list_64[i]*6)+'_'+'1.00'+ "_"+tol_str+"_"+ kappa_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+'_10.txt'    
        # print('Fname_p1: ',fname_p1)
        p1_data_temp = np.loadtxt(fname_p1)
        d_far_p1_64[i].append(p1_data_temp[2])
        d_loc_p1_64[i].append(p1_data_temp[1])

mean_dfar_p1_64 = np.mean(d_far_p1_64, axis=1)
std_dfar_p1_64 = np.std(d_far_p1_64, axis=1)
sem_dfar_p1_64 = std_dfar_p1_64/len(srand_list_64)

norm_mean_dfar_p1_64 = mean_dfar_p1_64/mean_dfar_p1_64[0]
norm_std_dfar_p1_64 = std_dfar_p1_64/mean_dfar_p1_64[0]

mean_dloc_p1_64 = np.mean(d_loc_p1_64, axis=1)
std_dloc_p1_64 = np.std(d_loc_p1_64, axis=1)

#==============================================================================
# putting dfar for p < 1 and p = 1 cases on the same plot
#==============================================================================

norm_d_far_p1 = mean_dfar_p1_64/mean_dfar_p1_64[0]
norm_d_far_p1_sem = sem_dfar_p1_64/mean_dfar_p1_64[0]
norm_d_far_e6 = mean_dfar_kappa[0]/mean_dfar_kappa[0,0]
norm_d_far_e6_sem = sem_dfar_kappa[0]/mean_dfar_kappa[0,0]
norm_d_far_e5 = mean_dfar_kappa[1,]/mean_dfar_kappa[1,0]
norm_d_far_e5_sem = sem_dfar_kappa[1]/mean_dfar_kappa[1,0]
norm_d_far_e4 = mean_dfar_kappa[2]/mean_dfar_kappa[2,0]
norm_d_far_e4_sem = sem_dfar_kappa[2]/mean_dfar_kappa[2,0]


plt.figure(figsize=(8, 5), dpi=300)

# # fitting p = 1 network
# x_linspace = np.linspace(num_center_list[0],num_center_list[-1],100) 
# line_fit = norm_d_far_p1[0]*x_linspace/num_center_list[0]
# plt.plot(x_linspace, line_fit, c = 'b') 

# fitting p < 1 network : N_d^3
x_linspace_k = np.linspace(num_center_list[5],num_center_list[-1],100) 
line_fit_k = (norm_d_far_e6[5]/num_center_list[5])*np.power(x_linspace_k,3)*0.003
plt.plot(x_linspace_k, line_fit_k, c = 'r') 
plt.text(20, 200, '$\sim N_d^3$', fontsize = 20)

# fitting p < 1 network : N_d^2
x_linspace_k = np.linspace(num_center_list[5],num_center_list[-1],100) 
line_fit_k = (norm_d_far_e6[5]/num_center_list[5])*np.power(x_linspace_k,2)*0.04
plt.plot(x_linspace_k, line_fit_k, c = 'b') 
plt.text(30, 60, '$\sim N_d^2$', fontsize = 20)

plt.errorbar(num_center_list, norm_d_far_p1, yerr=norm_d_far_p1_sem, marker = '.', ms = 25, alpha=0.5, label = '$p = 1$')# Network 1')  

plt.errorbar(num_center_list, norm_d_far_e6, yerr=norm_d_far_e6_sem, c = 'b', fmt = 'o', ms = 10, alpha=0.5, zorder = 10, label = '$\widetilde \kappa = 10^{-6}$')
plt.errorbar(num_center_list, norm_d_far_e5, yerr=norm_d_far_e5_sem, c = 'g', fmt = 'd', ms = 10, alpha=0.5, zorder = 8, label = '$\widetilde \kappa = 10^{-5}$')
plt.errorbar(num_center_list, norm_d_far_e4, yerr=norm_d_far_e4_sem, c = 'k', fmt = 's', ms = 10, alpha=0.5, zorder = 2, label = '$\widetilde \kappa = 10^{-4}$')

line_fit1 = norm_d_far_e6[0]*x_linspace/num_center_list[0]
# line_fit2 = norm_d_far_e5[0]*x_linspace/num_center_list[0]
# line_fit3 = norm_d_far_e4[0]*x_linspace/num_center_list[0]

plt.plot(x_linspace, line_fit1, c = 'k') 
# plt.plot(x_linspace, line_fit2, c = 'g') 
# plt.plot(x_linspace, line_fit3, c = 'k') 

plt.text(6, 3, '$\sim N_d^1$', fontsize = 20)   # for line_fit1

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
plt.savefig(plot_folder+"bend_dom_dip_moment_all_p_scatter_norm_"+ 'srand_112_through_122_' +pbond_string+"_"+kappa_str+"_"+tol_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+".png")

plt.clf()
plt.close()
