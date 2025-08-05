#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May  1 10:57:32 2025

focussed on three models; plots for paper

This code reads the values of dfar and energies written out by all_srand_all_nd_diff_network.py code
and reads them to make plots.

this is where the plot for the paper is created for dfar vs Nd and en ratio vs Nd!!

it also calculates the p_eff of all three models 

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

bending_flag = 0         # set this to 1 for looking at only bending dominated cases

# num_center_list = [5, 10, 15]
num_center_list = [1, 2, 5, 10, 15, 20, 25, 30, 35, 40]
# num_center_list = [1, 5, 10, 15, 20]    # for paper plot of dfar vs kappa
# num_center_list = [1, 2, 5, 10]

if L == 128:
    # num_center_list = [1, 8, 20, 40, 60, 80, 100, 120, 140, 160]
    num_center_list = [4, 8, 20, 40, 60, 80, 100, 120, 140, 160]

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
kappa_list = [1e-6]
if num_center_list[-1] == 20:
    kappa_list = [1e-6, 5e-6, 1e-5, 5e-5, 1e-4]
    
kappa = 1e-6

rlen = 0.9
rlen_txt = "%.4f" % rlen
rlen_txt_ring = "%.9f" % 0.1

srand = 112

pbond_string = "%.2f" % pbond # to write the filename
tol_str = "%.2e" % tol # to write the filename
    
mu_str = "%.4f" % mu # to write the filename
mu_c_str = "%.4f" % mu_c # to write the filename

base = "/home/abhinav/david/"

data_all = []
i = 0

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

kappa_str = "%.2e" % kappa # to write the filename
num_center = 40
num_dip = num_center*6

ranseed = '667720,601210'    # -- network 1

#==================================================================================================================================
# reading files*
#format of dfar  : heading = 'N         mean dfar         median dfar       std dfar         sem dfar'
#format of energy: heading = '        N          #sim   mean el en      median el en    std el en       mean bend en    median bend en   std bend en    mean en ratio   median en ratio  std en ratio'

folder1 = 'cluster_download/lattice_nopbd_bndry_all_clamp_restlength_circle_inner_fixed_radial_arp_bash_radial/'+ kappa_fname +ranseed+"/means/"
plot_folder1 = 'cluster_download/lattice_nopbd_bndry_all_clamp_restlength_circle_inner_fixed_radial_arp_bash_radial/'+ kappa_fname +ranseed+"/plots/"
fname1 = base + folder1 + "mean_dip_mom_vs_N_values_srand_113_122_"+pbond_string+"_"+tol_str+"_"+kappa_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+".txt"
print('reading dfar means : ', fname1)    
dip_mom_1 = np.loadtxt(fname1)
en_fname1 = base + folder1 + "mean_energy_vs_N_values_srand_113_122_"+pbond_string+"_"+tol_str+"_"+kappa_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+".txt"
print('reading energy means and medians : ', en_fname1)    
en_1 = np.loadtxt(en_fname1)

#==================================================================================================================================
folder2 = 'cluster_download/lattice_nopbd_bndry_all_clamp_restlength_circle_inner_fixed_radial_arp_bash_random_dipole_hex_bonds/'+ kappa_fname +ranseed+"/means/"
plot_folder2 = 'cluster_download/lattice_nopbd_bndry_all_clamp_restlength_circle_inner_fixed_radial_arp_bash_random_dipole_hex_bonds/'+ kappa_fname +ranseed+"/plots/"
fname2 = base + folder2 + "mean_dip_mom_vs_N_values_srand_113_122_"+pbond_string+"_"+tol_str+"_"+kappa_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+".txt"
print('reading dfar means : ', fname2)    
dip_mom_2 = np.loadtxt(fname2)
en_fname2 = base + folder2 + "mean_energy_vs_N_values_srand_113_122_"+pbond_string+"_"+tol_str+"_"+kappa_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+".txt"
print('reading energy means and medians : ', en_fname2)    
en_2 = np.loadtxt(en_fname2)

#==================================================================================================================================
folder3 = 'cluster_download/lattice_nopbd_bndry_all_clamp_restlength_circle_inner_fixed_radial_arp_bash_hex/'+ kappa_fname +ranseed+"/means/"
plot_folder3 = 'cluster_download/lattice_nopbd_bndry_all_clamp_restlength_circle_inner_fixed_radial_arp_bash_hex/'+ kappa_fname +ranseed+"/plots/"
fname3 = base + folder3 + "mean_dip_mom_vs_N_values_srand_113_122_"+pbond_string+"_"+tol_str+"_"+kappa_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+".txt"
print('reading dfar means : ', fname3)    
dip_mom_3 = np.loadtxt(fname3)
en_fname3 = base + folder3 + "mean_energy_vs_N_values_srand_113_122_"+pbond_string+"_"+tol_str+"_"+kappa_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+".txt"
print('reading energy means and medians : ', en_fname3)    
en_3 = np.loadtxt(en_fname3)


sem_en_1 = en_1[:,-1]/np.sqrt(en_1[:,1])
sem_en_2 = en_2[:,-1]/np.sqrt(en_2[:,1])
sem_en_3 = en_3[:,-1]/np.sqrt(en_3[:,1])


num_center_list = np.array(num_center_list)

# packing fraction for L = 64
unit_dipole_area = 6 * (np.sqrt(3)/4) * (1.5*1.5)     # 6 * area of unit traingle that makes a unit hexagon (note the side lenght  = 1.5 and not 1)
# unit_dipole_area = np.pi * ((2.5/2)*(2.5/2))     # 2.5 is the actual radius of no overlap
# area_of_inner_region_64 = 509.22   # this is gotten by 6* area of a triangle of side = 14. it is approximate total inner area
area_of_dipoles = num_center_list*unit_dipole_area
area_of_inner_region = np.pi * (12*12)   # this is almost circular

packing_frac = area_of_dipoles/area_of_inner_region

#==================================================================================================================================
# plotting the mean dipole moments vs number of dipoles for all kappa
plt.figure(figsize=(8, 5), dpi=300)
# plt.errorbar(num_center_list, dip_mom_1[:,1], yerr = dip_mom_1[:,-1], c = 'b', fmt = 's', ms = 12, alpha = 0.6, label = 'Model 1')
# plt.errorbar(num_center_list, dip_mom_2[:,1], yerr = dip_mom_2[:,-1], c = 'g', fmt = 'd', ms = 12, alpha = 0.6, label = 'Model 2')
# plt.errorbar(num_center_list, dip_mom_3[:,1], yerr = dip_mom_3[:,-1], c = 'r', fmt = 'o', ms = 12 , alpha = 0.6, label = 'Model 3')

plt.errorbar(packing_frac, dip_mom_1[:,1], yerr = dip_mom_1[:,-1], c = 'b', fmt = 's', ms = 12, alpha = 0.6, label = 'Model 1')
plt.errorbar(packing_frac, dip_mom_2[:,1], yerr = dip_mom_2[:,-1], c = 'g', fmt = 'd', ms = 12, alpha = 0.6, label = 'Model 2')
plt.errorbar(packing_frac, dip_mom_3[:,1], yerr = dip_mom_3[:,-1], c = 'r', fmt = 'o', ms = 12 , alpha = 0.6, label = 'Model 3')

# plt.xlabel("Number of Dipoles", fontsize = 20)
plt.xlabel("Dipole Packing Fraction, $\phi$", fontsize = 20)
plt.ylabel("$<D_{far}>$", fontsize = 20)
plt.xticks(fontsize = 20)
plt.yticks(fontsize = 20)
plt.yscale('log')
plt.xscale('log')

plt.subplots_adjust(right=0.98,left=0.18,top=0.98,bottom=0.16)
ax = plt.gca()
# Increase tick width and length for major and minor ticks
ax.tick_params(axis='both', which='major', width=3, length=10)  # Major ticks
ax.tick_params(axis='both', which='minor', width=2, length=5)   # Minor ticks

plt.legend(loc = 'best', fontsize = 15)
plt.savefig(base+"cluster_download/"+"dip_mom_means_vs_N_"+ str(srand) +"_" +pbond_string+"_"+tol_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+".png")

plt.clf()
plt.close()


# plotting the median dipole moments vs number of dipoles for all kappa
plt.figure(figsize=(8, 5), dpi=300)
plt.errorbar(num_center_list, dip_mom_1[:,2], yerr = dip_mom_1[:,-1], c = 'b', fmt = 's', ms = 12, alpha = 0.6, label = 'Model 1')
plt.errorbar(num_center_list, dip_mom_2[:,2], yerr = dip_mom_2[:,-1], c = 'g', fmt = 'd', ms = 12, alpha = 0.6, label = 'Model 2')
plt.errorbar(num_center_list, dip_mom_3[:,2], yerr = dip_mom_3[:,-1], c = 'r', fmt = 'o', ms = 12 , alpha = 0.6, label = 'Model 3')

plt.xlabel("Number of Dipoles", fontsize = 20)
plt.ylabel("Median  $D_{far}$", fontsize = 20)
plt.xticks(fontsize = 20)
plt.yticks(fontsize = 20)
plt.yscale('log')
plt.xscale('log')

plt.subplots_adjust(right=0.98,left=0.18,top=0.98,bottom=0.16)
ax = plt.gca()
# Increase tick width and length for major and minor ticks
ax.tick_params(axis='both', which='major', width=3, length=10)  # Major ticks
ax.tick_params(axis='both', which='minor', width=2, length=5)   # Minor ticks

plt.legend(loc = 'best', fontsize = 15)
plt.savefig(base+"cluster_download/"+"dip_mom_medians_vs_N_"+ str(srand) +"_" +pbond_string+"_"+tol_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+".png")

plt.clf()
plt.close()


# plotting the median energy ratios vs number of dipoles for all kappa
plt.figure(figsize=(8, 5), dpi=300)
plt.errorbar(num_center_list, en_1[:,9], yerr = sem_en_1, c = 'b', fmt = 's', ms = 12, alpha = 0.6, label = 'Model 1')
plt.errorbar(num_center_list, en_2[:,9], yerr = sem_en_2, c = 'g', fmt = 'd', ms = 12, alpha = 0.6, label = 'Model 2')
plt.errorbar(num_center_list, en_3[:,9], yerr = sem_en_3, c = 'r', fmt = 'o', ms = 12 , alpha = 0.6, label = 'Model 3')

plt.xlabel("Number of Dipoles", fontsize = 20)
plt.ylabel("Median $E_{st}/E_{bend}$", fontsize = 20)
plt.xticks(fontsize = 20)
plt.yticks(fontsize = 20)
plt.yscale('log')
plt.xscale('log')

plt.subplots_adjust(right=0.98,left=0.18,top=0.98,bottom=0.16)
ax = plt.gca()
# Increase tick width and length for major and minor ticks
ax.tick_params(axis='both', which='major', width=3, length=10)  # Major ticks
ax.tick_params(axis='both', which='minor', width=2, length=5)   # Minor ticks

plt.legend(loc = 'best', fontsize = 15)
plt.savefig(base+"cluster_download/"+"en_ratio_medians_vs_N_"+ str(srand) +"_" +pbond_string+"_"+tol_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+".png")

plt.clf()
plt.close()


# plotting the mean energy ratios vs number of dipoles for all kappa
plt.figure(figsize=(8, 5), dpi=300)
plt.errorbar(num_center_list, en_1[:,8], yerr = sem_en_1, c = 'b', fmt = 's', ms = 12, alpha = 0.6, label = 'Model 1')
plt.errorbar(num_center_list, en_2[:,8], yerr = sem_en_2, c = 'g', fmt = 'd', ms = 12, alpha = 0.6, label = 'Model 2')
plt.errorbar(num_center_list, en_3[:,8], yerr = sem_en_3, c = 'r', fmt = 'o', ms = 12 , alpha = 0.6, label = 'Model 3')

plt.xlabel("Number of Dipoles", fontsize = 20)
plt.ylabel("$<E_{st}/E_{bend}>$", fontsize = 20)
plt.xticks(fontsize = 20)
plt.yticks(fontsize = 20)
plt.yscale('log')
plt.xscale('log')

plt.subplots_adjust(right=0.98,left=0.18,top=0.98,bottom=0.16)
ax = plt.gca()
# Increase tick width and length for major and minor ticks
ax.tick_params(axis='both', which='major', width=3, length=10)  # Major ticks
ax.tick_params(axis='both', which='minor', width=2, length=5)   # Minor ticks

plt.legend(loc = 'best', fontsize = 15)
plt.savefig(base+"cluster_download/"+"en_ratio_means_vs_N_"+ str(srand) +"_" +pbond_string+"_"+tol_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+".png")

plt.clf()
plt.close()

# plotting the mean energy ratios vs number of dipoles for all kappa
plt.figure(figsize=(8, 5), dpi=300)
plt.errorbar(num_center_list, en_1[:,8], yerr = sem_en_1, c = 'b', fmt = 's', ms = 12, alpha = 0.6, label = '$\widetilde\kappa = 10^{-6}$')
plt.xlabel("Number of Dipoles", fontsize = 20)
plt.ylabel("$<E_{st}/E_{bend}>$", fontsize = 20)
plt.xticks(fontsize = 20)
plt.yticks(fontsize = 20)
plt.yscale('log')
plt.xscale('log')
plt.subplots_adjust(right=0.98,left=0.18,top=0.98,bottom=0.16)
ax = plt.gca()
# Increase tick width and length for major and minor ticks
ax.tick_params(axis='both', which='major', width=3, length=10)  # Major ticks
ax.tick_params(axis='both', which='minor', width=2, length=5)   # Minor ticks

plt.legend(loc = 'best', fontsize = 15)
plt.savefig(base+"cluster_download/"+"en_ratio_means_vs_N_model1_kapae-6_"+ str(srand) +"_" +pbond_string+"_"+tol_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+".png")

plt.clf()
plt.close()


# plotting the median energy ratios vs number of dipoles for all kappa
plt.figure(figsize=(8, 5), dpi=300)
plt.errorbar(num_center_list, en_1[:,9], yerr = sem_en_1, c = 'b', fmt = 's', ms = 12, alpha = 0.6, label = 'Model 1')
plt.errorbar(num_center_list, en_2[:,9], yerr = sem_en_2, c = 'g', fmt = 'd', ms = 12, alpha = 0.6, label = 'Model 2')
plt.errorbar(num_center_list, en_3[:,9], yerr = sem_en_3, c = 'r', fmt = 'o', ms = 12 , alpha = 0.6, label = 'Model 3')

plt.xlabel("Number of Dipoles", fontsize = 20)
plt.ylabel("Median $E_{st}/E_{bend}$", fontsize = 20)
plt.xticks(fontsize = 20)
plt.yticks(fontsize = 20)
plt.yscale('log')
plt.xscale('log')

plt.subplots_adjust(right=0.98,left=0.18,top=0.98,bottom=0.16)
ax = plt.gca()
# Increase tick width and length for major and minor ticks
ax.tick_params(axis='both', which='major', width=3, length=10)  # Major ticks
ax.tick_params(axis='both', which='minor', width=2, length=5)   # Minor ticks

plt.legend(loc = 'best', fontsize = 15)
plt.savefig(base+"cluster_download/"+"en_ratio_medians_vs_N_"+ str(srand) +"_" +pbond_string+"_"+tol_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+".png")

plt.clf()
plt.close()


#==========================================================================================================================
# Reading p = 1 L = 64 network to find p_eff
#==========================================================================================================================
fname_base  = '/home/abhinav/david/lattice_nopbd_bndry_all_clamp_restlength_circle_inner_fixed_radial_arp/txt/area/'
# fname_plot = fname_plot + str(L)+'_'   # works for L = 128 case only
srand_list = srand_list = [113, 114, 116, 117, 118, 119, 120, 121, 122]                         # 115 is BAD!!!!!

pbond_string = "%.2f" % 1 # to write the filename
kappa = 1e-6
kappa_str = "%.2e" % kappa # to write the filename

d_far_p1 = []
for i in range(0,len(num_center_list)):
    d_far_p1.append([])
    for j in range(0,len(srand_list)):
        srand = srand_list[j]
        fname_p1 = fname_base + 'dipole_moment_ratio_'+ 'srand_'+ str(srand) +'_'+str(num_center_list[i])+'_'+str(num_center_list[i]*6)+'_'+pbond_string+ "_"+tol_str+"_"+ kappa_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+'_10.txt'    
        p1_data_temp = np.loadtxt(fname_p1)
        d_far_p1[i].append(p1_data_temp[2])

mean_dfar_p1 = np.mean(d_far_p1, axis=1)
std_dfar_p1 = np.std(d_far_p1, axis=1)

#==========================================================================================================================
# need to use EMT to get p_eff from mu_m
#==========================================================================================================================
from emt_func import compute_p_eff

# for model 1
mu_m1 = dip_mom_1[:,1]/mean_dfar_p1
sigma_mu_m1 = mu_m1 * np.sqrt(np.square(dip_mom_1[:,-1]/dip_mom_1[:,1]) + np.square(std_dfar_p1/mean_dfar_p1))
p_eff_1, sigma_p_1 = compute_p_eff(mu_m1, sigma_mu_m1)

# for model 2
mu_m2 = dip_mom_2[:,1]/mean_dfar_p1
sigma_mu_m2 = mu_m2 * np.sqrt(np.square(dip_mom_2[:,-1]/dip_mom_2[:,1]) + np.square(std_dfar_p1/mean_dfar_p1))
p_eff_2, sigma_p_2 = compute_p_eff(mu_m2, sigma_mu_m2)

# for model 3
mu_m3 = dip_mom_3[:,1]/mean_dfar_p1
sigma_mu_m3 = mu_m3 * np.sqrt(np.square(dip_mom_3[:,-1]/dip_mom_3[:,1]) + np.square(std_dfar_p1/mean_dfar_p1))
p_eff_3, sigma_p_3 = compute_p_eff(mu_m3, sigma_mu_m3)

#==========================================================================================================================
# plotting
#==========================================================================================================================

fig = plt.figure(figsize=(8, 5), dpi=300)
plt.errorbar(packing_frac, mu_m1, yerr = sigma_mu_m1, c = 'b', ls = 'none', marker = 's', alpha = 0.7, ms = 13, label = 'Model 1')
plt.errorbar(packing_frac, mu_m2, yerr = sigma_mu_m2, c = 'g', ls = 'none', marker = 'd', alpha = 0.7, ms = 13, label = 'Model 2')
plt.errorbar(packing_frac, mu_m3, yerr = sigma_mu_m3, c = 'r', ls = 'none', marker = 'o', alpha = 0.7, ms = 13, label = 'Model 3')
# plt.axhline(y=2/3, color='k', linestyle='--', linewidth=2)
# plt.xlabel('$ N_d $', fontsize = 25)
plt.xlabel("Dipole Packing Fraction, $\phi$", fontsize = 20)
plt.ylabel(' $D_{far,p=0.55}/D_{far,p=1}$ ', fontsize = 20)
plt.xticks(fontsize = 20)
plt.yticks(fontsize = 20)
plt.subplots_adjust(top=0.96, left = 0.17, bottom = 0.16, right = 0.98)
plt.legend(loc = 'best', fontsize = 15)
plt.savefig(base+"cluster_download/"+"mu_m_srand_113_through_122_pooled_all_models_all_dom.png")  # just plotting the last step data
plt.close()
plt.clf()

fig = plt.figure(figsize=(8, 5), dpi=300)
# plt.errorbar(num_center_list, p_eff_1, yerr = sigma_p_1, c = 'b', ls = 'none', marker = 's', alpha = 0.7, ms = 13, label = 'Model 1')
# plt.errorbar(num_center_list, p_eff_2, yerr = sigma_p_2, c = 'g', ls = 'none', marker = 'd', alpha = 0.7, ms = 13, label = 'Model 2')
# plt.errorbar(num_center_list, p_eff_3, yerr = sigma_p_3, c = 'r', ls = 'none', marker = 'o', alpha = 0.7, ms = 13, label = 'Model 3')

plt.errorbar(packing_frac, p_eff_1, yerr = sigma_p_1, c = 'b', ls = 'none', marker = 's', alpha = 0.7, ms = 13, label = 'Model 1')
plt.errorbar(packing_frac, p_eff_2, yerr = sigma_p_2, c = 'g', ls = 'none', marker = 'd', alpha = 0.7, ms = 13, label = 'Model 2')
plt.errorbar(packing_frac, p_eff_3, yerr = sigma_p_3, c = 'r', ls = 'none', marker = 'o', alpha = 0.7, ms = 13, label = 'Model 3')
plt.axhline(y=2/3, color='k', linestyle='--', linewidth=2)
# plt.xlabel('$ N_d $', fontsize = 25)
plt.xlabel("Dipole Packing Fraction, $\phi$", fontsize = 20)
plt.ylabel(' $p_{eff}$ ', fontsize = 20)
plt.xticks(fontsize = 20)
plt.yticks(fontsize = 20)
plt.subplots_adjust(top=0.96, left = 0.17, bottom = 0.16, right = 0.98)
plt.legend(loc = 'best', fontsize = 15)
plt.savefig(base+"cluster_download/"+"peff_srand_113_through_122_pooled_all_models_all_dom.png")  # just plotting the last step data
plt.close()
plt.clf()



# finding the z_eff that corresponds to p_eff
z = 6   # for a triangular lattice at p=1
z_eff_1 = p_eff_1*z
z_eff_2 = p_eff_2*z
z_eff_3 = p_eff_3*z

sigma_z_1 = sigma_p_1*z
sigma_z_2 = sigma_p_2*z
sigma_z_3 = sigma_p_3*z

fig = plt.figure(figsize=(8, 5), dpi=300)
plt.errorbar(packing_frac, z_eff_1, yerr = sigma_z_1, c = 'b', ls = 'none', marker = 's', alpha = 0.7, ms = 13, label = 'Model 1')
plt.errorbar(packing_frac, z_eff_2, yerr = sigma_z_2, c = 'g', ls = 'none', marker = 'd', alpha = 0.7, ms = 13, label = 'Model 2')
plt.errorbar(packing_frac, z_eff_3, yerr = sigma_z_3, c = 'r', ls = 'none', marker = 'o', alpha = 0.7, ms = 13, label = 'Model 3')
plt.axhline(y=(2/3)*z, color='k', linestyle='--', linewidth=2)
# plt.xlabel('$ N_d $', fontsize = 25)
plt.xlabel("Dipole Packing Fraction, $\phi$", fontsize = 20)
plt.ylabel(' $z_{eff}$ ', fontsize = 25)
plt.xticks(fontsize = 25)
plt.yticks(fontsize = 25)
plt.subplots_adjust(top=0.96, left = 0.17, bottom = 0.15, right = 0.98)
plt.legend(loc = 'best', fontsize = 20)
plt.savefig(base+"cluster_download/"+"zeff_srand_113_through_122_pooled_all_models_all_dom.png")  # just plotting the last step data
plt.close()
plt.clf()


#==========================================================================================================================
# finding the slope of z_eff
#==========================================================================================================================
# Central difference for interior points :: slope is like d(z_eff)/d(N_d)
slope1 = np.empty_like(z_eff_1)
slope1_error = np.empty_like(z_eff_1)
delta_nd_central = num_center_list[2:] - num_center_list[:-2]
slope1[1:-1] = (z_eff_1[2:] - z_eff_1[:-2]) / delta_nd_central
# For the boundaries (forward/backward difference)
delta_nd_forward = num_center_list[1] - num_center_list[0]
slope1[0] = (z_eff_1[1] - z_eff_1[0]) / delta_nd_forward
# Backward difference (last point)
delta_nd_backward = num_center_list[-1] - num_center_list[-2]
slope1[-1] = (z_eff_1[-1] - z_eff_1[-2]) / delta_nd_backward
slope1_error[-1] = np.sqrt(sigma_z_1[-1]**2 + sigma_z_1[-2]**2) / np.abs(delta_nd_backward)
slope1_error[1:-1] = np.sqrt(sigma_z_1[2:]**2 + sigma_z_1[:-2]**2) / np.abs(delta_nd_central)
slope1_error[0] = np.sqrt(sigma_z_1[1]**2 + sigma_z_1[0]**2) / np.abs(delta_nd_forward)


# for model 2::
slope2 = np.empty_like(z_eff_2)
slope2_error = np.empty_like(z_eff_2)
# slope2[1:-1] = (z_eff_2[2:] - z_eff_2[:-2]) / (packing_frac[2:] - packing_frac[:-2])
slope2[1:-1] = (z_eff_2[2:] - z_eff_2[:-2]) / delta_nd_central
# For the boundaries (forward/backward difference)
slope2[0] = (z_eff_2[1] - z_eff_2[0]) / delta_nd_forward
slope2[-1] = (z_eff_2[-1] - z_eff_2[-2]) / delta_nd_backward
slope2_error[-1] = np.sqrt(sigma_z_2[-1]**2 + sigma_z_2[-2]**2) / np.abs(delta_nd_backward)
slope2_error[1:-1] = np.sqrt(sigma_z_2[2:]**2 + sigma_z_2[:-2]**2) / np.abs(delta_nd_central)
slope2_error[0] = np.sqrt(sigma_z_2[1]**2 + sigma_z_2[0]**2) / np.abs(delta_nd_forward)

# for model 3
slope3 = np.empty_like(z_eff_3)
slope3_error = np.empty_like(z_eff_3)
# slope3[1:-1] = (z_eff_3[2:] - z_eff_3[:-2]) / (packing_frac[2:] - packing_frac[:-2])
slope3[1:-1] = (z_eff_3[2:] - z_eff_3[:-2]) / delta_nd_central
# For the boundaries (forward/backward difference)
slope3[0] = (z_eff_3[1] - z_eff_3[0]) / delta_nd_forward
slope3[-1] = (z_eff_3[-1] - z_eff_3[-2]) / delta_nd_backward
slope3_error[-1] = np.sqrt(sigma_z_3[-1]**2 + sigma_z_3[-2]**2) / np.abs(delta_nd_backward)
slope3_error[1:-1] = np.sqrt(sigma_z_3[2:]**2 + sigma_z_3[:-2]**2) / np.abs(delta_nd_central)
slope3_error[0] = np.sqrt(sigma_z_3[1]**2 + sigma_z_3[0]**2) / np.abs(delta_nd_forward)


# plotting
fig = plt.figure(figsize=(8, 6), dpi=300)
plt.errorbar(num_center_list[2:], slope1[2:], yerr = slope1_error[2:], c = 'b', ls = 'none', marker = 's', alpha = 0.7, ms = 13, label = 'Model 1')
plt.errorbar(num_center_list[2:], slope2[2:], yerr = slope2_error[2:], c = 'g', ls = 'none', marker = 'd', alpha = 0.7, ms = 13, label = 'Model 2')
plt.errorbar(num_center_list[2:], slope3[2:], yerr = slope3_error[2:], c = 'r', ls = 'none', marker = 'o', alpha = 0.7, ms = 13, label = 'Model 3')
plt.xlabel("$N_d$", fontsize = 20)
plt.ylabel('Slope of $z_{eff}$', fontsize = 20)
plt.xticks(fontsize = 25)
plt.yticks(fontsize = 25)
plt.subplots_adjust(top=0.96, left = 0.17, bottom = 0.15, right = 0.98)
plt.savefig(base+"cluster_download/"+"zeff_slope_vs_Nd_srand_113_through_122_pooled_all_models_all_dom.png")  # just plotting the last step data
plt.close()
plt.clf()


#==========================================================================================================================
# converting the slope of z_eff to the actual constraint number
#==========================================================================================================================
N = 517   # nodes in the inner region
const = 2/(6*N)    # const * constraint = slope of peff vs Nd or zeff vs Nd

nc1 = slope1/const          # number of constraints for model 1
nc2 = slope2/const
nc3 = slope3/const

nc1_err = slope1_error/const
nc2_err = slope2_error/const
nc3_err = slope3_error/const

# plotting
fig = plt.figure(figsize=(8, 6), dpi=300)
plt.errorbar(num_center_list[2:], nc1[2:], yerr = nc1_err[2:], c = 'b', ls = 'none', marker = 's', alpha = 0.7, ms = 13, label = 'Model 1')
plt.errorbar(num_center_list[2:], nc2[2:], yerr = nc2_err[2:], c = 'g', ls = 'none', marker = 'd', alpha = 0.7, ms = 13, label = 'Model 2')
plt.errorbar(num_center_list[2:], nc3[2:], yerr = nc3_err[2:], c = 'r', ls = 'none', marker = 'o', alpha = 0.7, ms = 13, label = 'Model 3')
plt.xlabel("$N_d$", fontsize = 25)
plt.ylim([0,140])
plt.ylabel('Number of Constraints', fontsize = 25)
plt.xticks(fontsize = 25)
plt.yticks(fontsize = 25)
plt.subplots_adjust(top=0.96, left = 0.17, bottom = 0.15, right = 0.98)
plt.savefig(base+"cluster_download/"+"constraints_vs_Nd_srand_113_through_122_pooled_all_models_all_dom.png")  # just plotting the last step data
plt.close()
plt.clf()
