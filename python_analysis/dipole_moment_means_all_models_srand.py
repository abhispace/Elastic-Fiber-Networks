#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 14 21:11:58 2025

his code reads all of the dipole moments for all three models and plots them together.
First compare: srand = 112 data for all models
Then: use all srand values for the first model.

all_srand..._plotteer.y does the same thing but better!!
    
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

kappa_str = "%.2e" % kappa # to write the filename
num_center = 40
num_dip = num_center*6

ranseed = '667720,601210'    # -- network 1

#====================================================================================================================================================
# reading the data
#====================================================================================================================================================
# hex dipole data
folder = 'cluster_download/lattice_nopbd_bndry_all_clamp_restlength_circle_inner_fixed_radial_arp_bash_hex/'+ kappa_fname +ranseed+"/"+str(srand)+"/means/"
fname = base+folder+"all_dom_mean_dip_mom_vs_N_values_"+str(srand)+"_"+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+".txt"  
data_hex = np.loadtxt(fname)

# random transverse bonds dipole data
folder = 'cluster_download/lattice_nopbd_bndry_all_clamp_restlength_circle_inner_fixed_radial_arp_bash_random_dipole_hex_bonds/'+ kappa_fname +ranseed+"/"+str(srand)+"/means/"
fname = base+folder+"all_dom_mean_dip_mom_vs_N_values_"+str(srand)+"_"+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+".txt"
data_rand = np.loadtxt(fname)

# radial dipole data
folder = 'cluster_download/lattice_nopbd_bndry_all_clamp_restlength_circle_inner_fixed_radial_arp_bash_radial/'+ kappa_fname +ranseed+"/"+str(srand)+"/means/"
fname = base+folder+"all_dom_mean_dip_mom_vs_N_values_"+str(srand)+"_"+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+".txt"  
data_radial = np.loadtxt(fname)

plot_folder = 'cluster_download/'

data_hex = np.array(data_hex)
data_rand = np.array(data_rand)
data_radial = np.array(data_radial)

#====================================================================================================================================================
# reading the data
#====================================================================================================================================================


# # standard error of mean = standard deviation / sqrt(n)
# # yerr = data_all[:,:,2]/np.sqrt(data_all[:,:,0])
# yerr = data_all[:,:,2]

# # Following works only for kappae-6 case!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# # Calculate asymmetric errors in linear space
# y_lower = data_all[0,:,1] - yerr
# y_upper = data_all[0,:,1] + yerr
# # Ensure no negative values (log scale can't handle zero or negatives)
# y_lower = np.clip(y_lower, 1e-10, None)
# # Build asymmetric error array
# asymmetric_error = np.squeeze([data_all[0,:,1] - y_lower, y_upper - data_all[0,:,1]])

# plotting the dipole moments vs number of dipoles for all kappa
plt.figure(figsize=(8, 5), dpi=300)
plt.errorbar(num_center_list, data_hex[:,1], yerr = data_hex[:,2], c = 'r', fmt = 'o', label = 'Model 3')
plt.errorbar(num_center_list, data_rand[:,1], yerr = data_rand[:,2], c = 'g', fmt = 'd', label = 'Model 2')
plt.errorbar(num_center_list, data_radial[:,1], yerr = data_radial[:,2], c = 'b', fmt = 's', label = 'Model 1')

# x_linspace = np.linspace(num_center_list[0],num_center_list[-1],100) 
# line_fit1 = data_all[0,0,1]*x_linspace/num_center_list[0]
# line_fit2 = data_all[1,0,1]*x_linspace/num_center_list[0]
# line_fit3 = data_all[2,0,1]*x_linspace/num_center_list[0]

# plt.plot(x_linspace, line_fit1, c = 'b') 
# plt.plot(x_linspace, line_fit2, c = 'g') 
# plt.plot(x_linspace, line_fit3, c = 'k') 

plt.xlabel("Number of Dipoles", fontsize = 15)
plt.ylabel("$<D_{far}>$", fontsize = 15)
plt.xticks(fontsize = 15)
plt.yticks(fontsize = 15)
plt.yscale('log')
plt.xscale('log')

plt.subplots_adjust(right=0.96,left=0.18,top=0.94,bottom=0.12)

plt.legend(loc = 'best')
plt.savefig(base+plot_folder+"all_models_all_dom_dip_mom_vs_N_"+ str(srand) +"_" +pbond_string+"_"+tol_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+".png")

plt.clf()
plt.close()


log_num_dip = np.log10(num_center_list)
log_dfar_hex = np.log10(data_hex[:,1])
log_yerr_hex = data_hex[:,2] / (data_hex[:,1] * np.log10(10))  # error propagation in log10
log_yerr_sem_hex = data_hex[:,3] / (data_hex[:,1] * np.log10(10))  # error propagation in log10

log_dfar_rand = np.log10(data_rand[:,1])
log_yerr_rand = data_rand[:,2] / (data_rand[:,1] * np.log10(10))  # error propagation in log10
log_yerr_sem_rand = data_rand[:,3] / (data_rand[:,1] * np.log10(10))  # error propagation in log10

log_dfar_radial = np.log10(data_radial[:,1])
log_yerr_radial = data_radial[:,2] / (data_radial[:,1] * np.log10(10))  # error propagation in log10
log_yerr_sem_radial = data_radial[:,3] / (data_radial[:,1] * np.log10(10))  # error propagation in log10

# plotting the dipole moments vs number of dipoles for all kappa
plt.figure(figsize=(8, 5), dpi=300)
plt.errorbar(log_num_dip, log_dfar_hex, yerr = log_yerr_hex, c = 'r', fmt = 'o', label = 'Model 3')
plt.errorbar(log_num_dip, log_dfar_rand, yerr = log_yerr_rand, c = 'g', fmt = 'd', label = 'Model 2')
plt.errorbar(log_num_dip, log_dfar_radial, yerr = log_yerr_radial, c = 'b', fmt = 's', label = 'Model 1')

# x_linspace = np.linspace(num_center_list[0],num_center_list[-1],100) 
# line_fit1 = data_all[0,0,1]*x_linspace/num_center_list[0]

# plt.plot(x_linspace, line_fit1, c = 'b') 
# plt.plot(x_linspace, line_fit2, c = 'g') 
# plt.plot(x_linspace, line_fit3, c = 'k') 

plt.xlabel("$ log(N_{d}) $", fontsize = 15)
plt.ylabel("$ log(<D_{far}>) $", fontsize = 15)
plt.xticks(fontsize = 15)
plt.yticks(fontsize = 15)
plt.subplots_adjust(right=0.96,left=0.18,top=0.94,bottom=0.12)
plt.legend(loc = 'best')
plt.savefig(base+plot_folder+"all_models_all_dom_dip_mom_vs_N_logged_"+ str(srand) +"_" +pbond_string+"_"+tol_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+".png")
plt.clf()
plt.close()


# plotting the dipole moments vs number of dipoles for all kappa
plt.figure(figsize=(8, 5), dpi=300)
plt.errorbar(log_num_dip, log_dfar_hex, yerr = log_yerr_sem_hex, c = 'r', fmt = 'o', label = 'Model 3')
plt.errorbar(log_num_dip, log_dfar_rand, yerr = log_yerr_sem_rand, c = 'g', fmt = 'd', label = 'Model 2')
plt.errorbar(log_num_dip, log_dfar_radial, yerr = log_yerr_sem_radial, c = 'b', fmt = 's', label = 'Model 1')

# x_linspace = np.linspace(num_center_list[0],num_center_list[-1],100) 
# line_fit1 = data_all[0,0,1]*x_linspace/num_center_list[0]

# plt.plot(x_linspace, line_fit1, c = 'b') 
# plt.plot(x_linspace, line_fit2, c = 'g') 
# plt.plot(x_linspace, line_fit3, c = 'k') 

plt.xlabel("$ log_{10}(N_{d}) $", fontsize = 15)
plt.ylabel("$ log_{10}(<D_{far}>) $", fontsize = 15)
plt.xticks(fontsize = 15)
plt.yticks(fontsize = 15)
plt.subplots_adjust(right=0.96,left=0.18,top=0.94,bottom=0.12)
plt.legend(loc = 'best')
plt.savefig(base+plot_folder+"all_models_all_dom_dip_mom_vs_N_logged_SEM_"+ str(srand) +"_" +pbond_string+"_"+tol_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+".png")
plt.clf()
plt.close()


##======================================================================================================================================================================
# for all srand values!!
##======================================================================================================================================================================

srand_list = ['113', '114','116','117','118','119','120','121','122']                         # 115 is BAD!!!!!

data_hex_srand = []
data_radial_srand = []
data_rand_srand = []
data_hex_en_srand = []
data_radial_en_srand = []
data_rand_en_srand = []

for srand in srand_list:
    # reading the data
    # hex dipole data
    folder = 'cluster_download/lattice_nopbd_bndry_all_clamp_restlength_circle_inner_fixed_radial_arp_bash_hex/'+ kappa_fname +ranseed+"/"+str(srand)+"/means/"
    fname = base+folder+"all_dom_mean_dip_mom_vs_N_values_"+str(srand)+"_"+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+".txt"  
    data_hex_srand.append(np.loadtxt(fname))
    fname = base+folder+"all_dom_mean_energy_vs_N_values_"+str(srand)+"_"+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+".txt"  
    data_hex_en_srand.append(np.loadtxt(fname))
    
    # # random transverse bonds dipole data
    folder = 'cluster_download/lattice_nopbd_bndry_all_clamp_restlength_circle_inner_fixed_radial_arp_bash_random_dipole_hex_bonds/'+ kappa_fname +ranseed+"/"+str(srand)+"/means/"
    fname = base+folder+"all_dom_mean_dip_mom_vs_N_values_"+str(srand)+"_"+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+".txt"
    data_rand_srand.append(np.loadtxt(fname))
    fname = base+folder+"all_dom_mean_energy_vs_N_values_"+str(srand)+"_"+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+".txt"
    data_rand_en_srand.append(np.loadtxt(fname))
    
    # radial dipole data
    folder = 'cluster_download/lattice_nopbd_bndry_all_clamp_restlength_circle_inner_fixed_radial_arp_bash_radial/'+ kappa_fname +ranseed+"/"+str(srand)+"/means/"
    fname = base+folder+"all_dom_mean_dip_mom_vs_N_values_"+str(srand)+"_"+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+".txt"  
    data_radial_srand.append(np.loadtxt(fname))
    fname = base+folder+"all_dom_mean_energy_vs_N_values_"+str(srand)+"_"+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+".txt"  
    data_radial_en_srand.append(np.loadtxt(fname))
    
    plot_folder = 'cluster_download/'
    
data_hex_srand = np.array(data_hex_srand)
data_rand_srand = np.array(data_rand_srand)
data_radial_srand = np.array(data_radial_srand)

data_hex_en_srand = np.array(data_hex_en_srand)
data_rand_en_srand = np.array(data_rand_en_srand)
data_radial_en_srand = np.array(data_radial_en_srand)

# calculating mean values
mean_dfar_hex_srand = np.mean(data_hex_srand[:,:,1], axis = 0)
mean_dfar_radial_srand = np.mean(data_radial_srand[:,:,1], axis = 0)
mean_dfar_rand_srand = np.mean(data_rand_srand[:,:,1], axis = 0)
    
mean_en_ratio_hex_srand = np.mean(data_hex_en_srand[:,:,6], axis = 0)
mean_en_ratio_radial_srand = np.mean(data_radial_en_srand[:,:,6], axis = 0)
mean_en_ratio_rand_srand = np.mean(data_rand_en_srand[:,:,6], axis = 0)

# calculating sem values
sem_dfar_hex_srand = np.std(data_hex_srand[:,:,1], axis = 0)/np.sqrt(np.sum(data_hex_en_srand[:,:,1], axis = 0))
sem_dfar_radial_srand = np.std(data_radial_srand[:,:,1], axis = 0)/np.sqrt(np.sum(data_radial_en_srand[:,:,1], axis = 0))
sem_dfar_rand_srand = np.std(data_rand_srand[:,:,1], axis = 0)/np.sqrt(np.sum(data_rand_en_srand[:,:,1], axis = 0))

sem_en_ratio_hex_srand = np.std(data_hex_en_srand[:,:,6], axis = 0)/np.sqrt(np.sum(data_hex_en_srand[:,:,1], axis = 0))
sem_en_ratio_radial_srand = np.std(data_radial_en_srand[:,:,6], axis = 0)/np.sqrt(np.sum(data_radial_en_srand[:,:,1], axis = 0))
sem_en_ratio_rand_srand = np.std(data_rand_en_srand[:,:,6], axis = 0)/np.sqrt(np.sum(data_rand_en_srand[:,:,1], axis = 0))

#-----------------------------------------------------------------------------------------------------------------
# plotting the dipole moments vs number of dipoles for all kappa
plt.figure(figsize=(8, 5), dpi=300)
plt.errorbar(num_center_list, mean_dfar_hex_srand, yerr = sem_dfar_hex_srand, c = 'r', fmt = 'o', ms = 12, alpha = 0.6, label = 'Model 3')
plt.errorbar(num_center_list, mean_dfar_rand_srand, yerr = sem_dfar_rand_srand, c = 'g', fmt = 'd', ms = 12, alpha = 0.6, label = 'Model 2')
plt.errorbar(num_center_list, mean_dfar_radial_srand, yerr = sem_dfar_radial_srand, c = 'b', fmt = 's', ms = 12 , alpha = 0.6, label = 'Model 1')

plt.xlabel("Number of Dipoles", fontsize = 20)
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
plt.savefig(base+plot_folder+"all_models_all_dom_dip_mom_vs_N_srand113_through122_"+ str(srand) +"_" +pbond_string+"_"+tol_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+".png")

plt.clf()
plt.close()


# plotting the energy ratio vs number of dipoles for all kappa
plt.figure(figsize=(8, 5), dpi=300)
plt.errorbar(num_center_list, mean_en_ratio_hex_srand, yerr = sem_en_ratio_hex_srand, c = 'r', fmt = 'o', ms = 12, alpha = 0.6, label = 'Model 3')
plt.errorbar(num_center_list, mean_en_ratio_rand_srand, yerr = sem_en_ratio_rand_srand, c = 'g', fmt = 'd', ms = 12, alpha = 0.6, label = 'Model 2')
plt.errorbar(num_center_list, mean_en_ratio_radial_srand, yerr = sem_en_ratio_radial_srand, c = 'b', fmt = 's', ms = 12 , alpha = 0.6, label = 'Model 1')

plt.xlabel("Number of Dipoles", fontsize = 20)
plt.ylabel("$<E_{st}/E_{b}>$", fontsize = 20)
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
plt.savefig(base+plot_folder+"all_models_all_dom_energy_ratio_vs_N_srand113_through122_"+pbond_string+"_"+tol_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+".png")

plt.clf()
plt.close()
