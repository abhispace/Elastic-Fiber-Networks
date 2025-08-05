#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 10 11:09:00 2025

This code reads the dfar for L = 128 and calculates the p_eff.
It also reads dfar for L = 64 to make a master plot with both L data after normalization.

this code also verifies the previous p_eff plot made at the end of network_circular_caller.py.
however it does so a little differently than that code.
that code used the already calculated p_eff for different srands and took a mean to get the p_eff vs N_d plot.
this code actually reads in all the individual dfar values and then averages them.
then using this gloabl average dfar done over all srand cases, it calculates mu_m for each N_d case.
then after averaging the mu_m, it calls the P_eff calculating function and plots p_eff vs N_d.

NOTE: at the bottom it also makes plots for figure 3 in the paper!!

this code uses the means of the srand cases and then averages them!!

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

from matplotlib.ticker import FuncFormatter
# Define the formatter function
def scientific_format(x, pos):
    if x == 0:
        return "0"
    exponent = int(np.floor(np.log10(abs(x))))
    coefficient = x / 10**exponent
    return r"${:.1f} \times 10^{{{}}}$".format(coefficient, exponent)

L = 64
L = 128
num_pts = L*L

# bending_flag = 0    # set it to 1 to use only bending dominated cases: applies only for L = 64 for now!!

# num_center_list = [5, 10, 15]
num_center_list = [1, 2, 5, 10, 15, 20, 25, 30, 35, 40]
# num_center_list = [1, 5, 10, 15, 20]    # for paper plot of dfar vs kappa
# num_center_list = [1, 2, 5, 10]

if L == 64:
    # srand_list = ['112','113', '114','116','117','118','119','120','121','122']                         # 115 is BAD!!!!!
    srand_list = ['113', '114','116','117','118','119','120','121','122']                         # 115 is BAD!!!!!
if L == 128:
    num_center_list1 = [1, 4, 8, 20, 40, 60]
    num_center_list2 = [80, 100, 120, 140, 160]
    # srand_list1 = ['116','117','118','119','120']                         # 115 is BAD!!!!!
    # srand_list2 = ['116','117']
    # srand_list2 = ['116','117','118','119','120']
    # srand_list2 = ['114','115','116','117','118','119','120','121','122']
    # for the paper
    srand_list1 = ['113','114','115','116','117','118','119','120','121','122']                         # 115 is BAD!!!!!
    srand_list2 = ['114','115','116','117','121','122']
    diff_network_list = ['111111,101111','111111,111111']

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
# tol = 1.0e-8


rlen = 0.9
rlen_txt = "%.4f" % rlen
rlen_txt_ring = "%.9f" % 0.1


pbond_string = "%.2f" % pbond # to write the filename
tol_str = "%.2e" % tol # to write the filename
    
mu_str = "%.4f" % mu # to write the filename
mu_c_str = "%.4f" % mu_c # to write the filename

base = "/home/abhinav/david/"

cluster_new_flag = 1   # for dipoles with only radial bonds
hex_flag = 0           # for dipoles with hexagonal bonds forced to be present
hex_rand_flag = 0      # for dipoles with hexagonal bonds randomly removed as per p value
cluster_radial_flag = 1   # for dipoles with just radial bonds

i = 0

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

kappa_str = "%.2e" % kappa # to write the filename

#==========================================================================================================================
# for the L = 128; new method
#==========================================================================================================================
num_center_list = np.concatenate((num_center_list1, num_center_list2))
data_all = []
data_dip_mom = []
data_dip_mom1_bend = []

en_data_all1 = []
en_st_all1 = []
en_bend_all1 = []
en_ratio_all1 = []

i = 0
ranseed = '667720,601210'    # -- network 1

for i in range(0,len(num_center_list1)):

    data_all.append([])
    en_data_all1.append([])
    en_st_all1.append([])
    en_bend_all1.append([])
    en_ratio_all1.append([])
    data_dip_mom.append([])
    data_dip_mom1_bend.append([])
    num_center = num_center_list1[i]
    num_dip = num_center*6    

    for j in range(0,len(srand_list1)):        
        srand = srand_list1[j]
        data_all[i].append([])
        data_dip_mom[i].append([])
        data_dip_mom1_bend[i].append([])
        en_data_all1[i].append([])
        en_st_all1[i].append([])
        en_bend_all1[i].append([])
        en_ratio_all1[i].append([])

        for k in range(0,len(diff_network_list)):
            diff_network = diff_network_list[k]
            
            folder = 'cluster_download/lattice_nopbd_bndry_all_clamp_restlength_circle_inner_fixed_radial_arp_bash_radial/'+ kappa_fname +ranseed+"/"+str(srand)+"/"+ diff_network + "/"
            plot_folder = 'cluster_download/lattice_nopbd_bndry_all_clamp_restlength_circle_inner_fixed_radial_arp_bash_radial/'+ kappa_fname +ranseed+"/plots/"
        
            fname = base+folder+"txt/area/"+"dipole_moment_ratio_"+ str(L) + "_" +str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+".txt"
            # energy file
            en_fname = base+folder+"energy/strain/"+"Lattice_"+str(L)+"_"+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+"_force.txt"    

            if (num_center == 8 and srand == '119' and diff_network == '111111,111111' and tol == 1e-7):
                continue
            elif (num_center == 40 and (srand == '119') and diff_network == '111111,101111' and tol == 1e-8):
                continue
            else:
                data = np.loadtxt(fname)
                data_all[i][j].append(data)        
                data_dip_mom[i][j].append(data[2])

                en_data = np.loadtxt(en_fname)
                en_data_all1[i][j].append(en_data[-1])        
                en_st_all1[i][j].append(en_data[-1,1] + en_data[-1,2])        # stretch plus compression energies
                en_bend_all1[i][j].append(en_data[-1,3])        
                en_ratio = (en_data[-1,1] + en_data[-1,2])/en_data[-1,3]
                en_ratio_all1[i][j].append(en_ratio)        
                
                if en_ratio < 1:
                    data_dip_mom1_bend[i][j].append(data[2])
                # else:
                #     print('Not bending dominated!!', num_center_list1[i], srand_list1[j], diff_network_list[k])

dip_mom_1d_128_1 = []
median_dip_mom_128_1 = []
mean_dip_mom_128_1 = []
std_dip_mom_128_1 = []

dip_mom_bend_1d_128_1 = []
median_dip_mom_bend_128_1 = []
mean_dip_mom_bend_128_1 = []
std_dip_mom_bend_128_1 = []

en_ratio_1d_128_1 = []
median_en_ratio_128_1 = []
mean_en_ratio_128_1 = []
std_en_ratio_128_1 = []

for i in range(0,len(num_center_list1)):
    num_center = num_center_list1[i]
    num_dip = num_center*6
        
    dip_mom_1d_128_1.append([item for sublist in data_dip_mom[i] for item in sublist])    
    median_dip_mom_128_1.append(np.median(dip_mom_1d_128_1[i]))
    mean_dip_mom_128_1.append(np.mean(dip_mom_1d_128_1[i]))
    std_dip_mom_128_1.append(np.std(dip_mom_1d_128_1[i]))

    en_ratio_1d_128_1.append([item for sublist in en_ratio_all1[i] for item in sublist])    
    median_en_ratio_128_1.append(np.median(en_ratio_1d_128_1[i]))
    mean_en_ratio_128_1.append(np.mean(en_ratio_1d_128_1[i]))
    std_en_ratio_128_1.append(np.std(en_ratio_1d_128_1[i]))

    # just the bending dominated networks here
    dip_mom_bend_1d_128_1.append([item for sublist in data_dip_mom1_bend[i] for item in sublist])    
    median_dip_mom_bend_128_1.append(np.median(dip_mom_bend_1d_128_1[i]))
    mean_dip_mom_bend_128_1.append(np.mean(dip_mom_bend_1d_128_1[i]))
    std_dip_mom_bend_128_1.append(np.std(dip_mom_bend_1d_128_1[i]))

mean_dfar = []
std_dfar = []
for i in range(0,len(mean_dip_mom_128_1)):
    mean_dfar.append(mean_dip_mom_128_1[i])
    std_dfar.append(std_dip_mom_128_1[i])

mean_dfar_bend = []
std_dfar_bend = []
for i in range(0,len(mean_dip_mom_bend_128_1)):
    mean_dfar_bend.append(mean_dip_mom_bend_128_1[i])
    std_dfar_bend.append(std_dip_mom_bend_128_1[i])

mean_en_ratio = []
std_en_ratio = []
for i in range(0,len(mean_en_ratio_128_1)):
    mean_en_ratio.append(mean_en_ratio_128_1[i])
    std_en_ratio.append(std_en_ratio_128_1[i])

data_all2 = []
data_dip_mom2 = []
data_dip_mom2_bend = []

en_data_all2 = []
en_st_all2 = []
en_bend_all2 = []
en_ratio_all2 = []

for i in range(0,len(num_center_list2)):
    data_all2.append([])
    data_dip_mom2.append([])
    data_dip_mom2_bend.append([])
    en_data_all2.append([])
    en_st_all2.append([])
    en_bend_all2.append([])
    en_ratio_all2.append([])
    num_center = num_center_list2[i]
    num_dip = num_center*6    
    for j in range(0,len(srand_list2)):        
        srand = srand_list2[j]
        data_all2[i].append([])
        data_dip_mom2[i].append([])
        data_dip_mom2_bend[i].append([])
        en_data_all2[i].append([])
        en_st_all2[i].append([])
        en_bend_all2[i].append([])
        en_ratio_all2[i].append([])

        for k in range(0,len(diff_network_list)):
            diff_network = diff_network_list[k]
                        
            folder = 'cluster_download/lattice_nopbd_bndry_all_clamp_restlength_circle_inner_fixed_radial_arp_bash_radial/'+ kappa_fname +ranseed+"/"+str(srand)+"/"+ diff_network + "/"
            plot_folder = 'cluster_download/lattice_nopbd_bndry_all_clamp_restlength_circle_inner_fixed_radial_arp_bash_radial/'+ kappa_fname +ranseed+"/plots/"
        
            fname = base+folder+"txt/area/"+"dipole_moment_ratio_"+ str(L) + "_" +str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+".txt"
            # energy file
            en_fname = base+folder+"energy/strain/"+"Lattice_"+str(L)+"_"+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+"_force.txt"    
        
            if (num_center == 80 and srand == '116' and diff_network == '111111,111111'):
                continue
            elif (num_center == 160 and srand == '114' and diff_network == '111111,111111'):
                continue
            elif (num_center == 160 and srand == '115' and diff_network == '111111,101111'):
                continue
            elif (num_center == 160 and (srand == '114' or srand == '115' or srand == '121' or srand == '122')):
                continue
            else:
                data = np.loadtxt(fname)
                data_all2[i][j].append(data)        
                data_dip_mom2[i][j].append(data[2])
                
                en_data = np.loadtxt(en_fname)
                en_data_all2[i][j].append(en_data[-1])        
                en_st_all2[i][j].append(en_data[-1,1] + en_data[-1,2])        # stretch plus compression energies
                en_bend_all2[i][j].append(en_data[-1,3])        
                en_ratio = (en_data[-1,1] + en_data[-1,2])/en_data[-1,3]
                en_ratio_all2[i][j].append(en_ratio)        
                
                if en_ratio < 1:
                    data_dip_mom2_bend[i][j].append(data[2])
                # else:
                #     print('Not bending dominated!!', num_center_list2[i], srand_list2[j], diff_network_list[k])
                    

dip_mom_1d_128_2 = []
median_dip_mom_128_2 = []
mean_dip_mom_128_2 = []
std_dip_mom_128_2 = []

dip_mom_bend_1d_128_2 = []
median_dip_mom_bend_128_2 = []
mean_dip_mom_bend_128_2 = []
std_dip_mom_bend_128_2 = []

en_ratio_1d_128_2 = []
median_en_ratio_128_2 = []
mean_en_ratio_128_2 = []
std_en_ratio_128_2 = []

for i in range(0,len(num_center_list2)):
    num_center = num_center_list2[i]
    num_dip = num_center*6
        
    dip_mom_1d_128_2.append([item for sublist in data_dip_mom2[i] for item in sublist])    
    median_dip_mom_128_2.append(np.median(dip_mom_1d_128_2[i]))
    mean_dip_mom_128_2.append(np.mean(dip_mom_1d_128_2[i]))
    std_dip_mom_128_2.append(np.std(dip_mom_1d_128_2[i]))

    dip_mom_bend_1d_128_2.append([item for sublist in data_dip_mom2_bend[i] for item in sublist])    
    median_dip_mom_bend_128_2.append(np.median(dip_mom_bend_1d_128_2[i]))
    mean_dip_mom_bend_128_2.append(np.mean(dip_mom_bend_1d_128_2[i]))
    std_dip_mom_bend_128_2.append(np.std(dip_mom_bend_1d_128_2[i]))

    en_ratio_1d_128_2.append([item for sublist in en_ratio_all2[i] for item in sublist])    
    median_en_ratio_128_2.append(np.median(en_ratio_1d_128_2[i]))
    mean_en_ratio_128_2.append(np.mean(en_ratio_1d_128_2[i]))
    std_en_ratio_128_2.append(np.std(en_ratio_1d_128_2[i]))

# mean_dfar.append(mean_dip_mom_128_2)
# std_dfar.append(std_dip_mom_128_2)
for i in range(0,len(mean_dip_mom_128_2)):
    mean_dfar.append(mean_dip_mom_128_2[i])
    std_dfar.append(std_dip_mom_128_2[i])

mean_dfar = np.array(mean_dfar)
std_dfar = np.array(std_dfar)

for i in range(0,len(mean_dip_mom_bend_128_2)):
    mean_dfar_bend.append(mean_dip_mom_bend_128_2[i])
    std_dfar_bend.append(std_dip_mom_bend_128_2[i])

mean_dfar_bend = np.array(mean_dfar_bend)
std_dfar_bend = np.array(std_dfar_bend)

for i in range(0,len(mean_en_ratio_128_2)):
    mean_en_ratio.append(mean_en_ratio_128_2[i])
    std_en_ratio.append(std_en_ratio_128_2[i])

mean_en_ratio = np.array(mean_en_ratio)
std_en_ratio = np.array(std_en_ratio)

# plotting the dipole moments vs number of dipoles for all kappa
plt.figure(figsize=(8, 5), dpi=300)
plt.errorbar(num_center_list, mean_dfar, yerr = std_dfar, c = 'b', fmt = 'o', ms = 10, label = '$\widetilde \kappa = 10^{-6}$')
# plt.errorbar(num_center_list, dfar_norm, yerr = yerr, c = 'b', fmt = 'o', ms = 10, label = '$\widetilde \kappa = 10^{-6}$')

x_linspace = np.linspace(num_center_list[0],num_center_list[-1],100) 
line_fit1 = mean_dfar[0]*x_linspace/num_center_list[0]
plt.plot(x_linspace, line_fit1, c = 'b', lw = 3) 

x_linspace2 = np.linspace(num_center_list[-5],num_center_list[-1],100) 
line_fit2 = mean_dfar[-5]*np.power(x_linspace2,2)/num_center_list[-5]*0.008
plt.plot(x_linspace2, line_fit2, c = 'g', lw = 3) 

x_linspace3 = np.linspace(num_center_list[-5],num_center_list[-1],100) 
line_fit3 = mean_dfar[-5]*np.power(x_linspace3,3)/num_center_list[-5]*0.00015
plt.plot(x_linspace3, line_fit3, c = 'r', lw = 3) 

plt.text(10, 3e-5, '$\sim N_d^1$', fontsize = 20)
plt.text(60, 5e-3, '$\sim N_d^3$', fontsize = 20)
plt.text(110, 2e-3, '$\sim N_d^2$', fontsize = 20)

plt.xlabel("Number of Dipoles", fontsize = 20)
plt.ylabel("$<D_{far}>$", fontsize = 20)
plt.xticks(fontsize = 20)
plt.yticks(fontsize = 20)
plt.yscale('log')
plt.xscale('log')

plt.subplots_adjust(right=0.98,left=0.15,top=0.98,bottom=0.16)

ax = plt.gca()
# Increase tick width and length for major and minor ticks
ax.tick_params(axis='both', which='major', width=3, length=10)  # Major ticks
ax.tick_params(axis='both', which='minor', width=2, length=5)   # Minor ticks

plt.legend(loc = 'best', fontsize = 15)
plt.savefig(base+plot_folder+"dip_mom_vs_N_L_all_dom_"+ str(L) +"_srand_pooled_" +pbond_string+"_"+tol_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+".png")

plt.clf()
plt.close()



# plotting the bending dominated cases: dipole moments vs number of dipoles for all kappa
plt.figure(figsize=(8, 5), dpi=300)
plt.errorbar(num_center_list, mean_dfar_bend, yerr = std_dfar_bend, c = 'b', fmt = 'o', ms = 10, label = '$\widetilde \kappa = 10^{-6}$')
# plt.errorbar(num_center_list, dfar_norm, yerr = yerr, c = 'b', fmt = 'o', ms = 10, label = '$\widetilde \kappa = 10^{-6}$')
x_linspace = np.linspace(num_center_list[0],num_center_list[-1],100) 
line_fit1 = mean_dfar_bend[0]*x_linspace/num_center_list[0]
plt.plot(x_linspace, line_fit1, c = 'b', lw = 3) 

# x_linspace2 = np.linspace(num_center_list[-5],num_center_list[-1],100) 
# line_fit2 = mean_dfar_bend[-5]*np.power(x_linspace2,2)/num_center_list[-5]*0.008
# plt.plot(x_linspace2, line_fit2, c = 'g', lw = 3) 

# x_linspace3 = np.linspace(num_center_list[-5],num_center_list[-1],100) 
# line_fit3 = mean_dfar_bend[-5]*np.power(x_linspace3,3)/num_center_list[-5]*0.00015
# plt.plot(x_linspace3, line_fit3, c = 'r', lw = 3) 

plt.text(10, 3e-5, '$\sim N_d^1$', fontsize = 20)
# plt.text(60, 5e-3, '$\sim N_d^3$', fontsize = 20)
# plt.text(110, 2e-3, '$\sim N_d^2$', fontsize = 20)

plt.xlabel("Number of Dipoles", fontsize = 20)
plt.ylabel("$<D_{far}>$", fontsize = 20)
plt.xticks(fontsize = 20)
plt.yticks(fontsize = 20)
plt.yscale('log')
plt.xscale('log')
plt.subplots_adjust(right=0.98,left=0.15,top=0.98,bottom=0.16)
ax = plt.gca()
# Increase tick width and length for major and minor ticks
ax.tick_params(axis='both', which='major', width=3, length=10)  # Major ticks
ax.tick_params(axis='both', which='minor', width=2, length=5)   # Minor ticks
plt.legend(loc = 'best', fontsize = 15)
plt.savefig(base+plot_folder+"dip_mom_vs_N_L_bend_dom_"+ str(L) +"_srand_pooled_" +pbond_string+"_"+tol_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+".png")
plt.clf()
plt.close()




# plotting the scatter plot of all bending dominated cases: for each Nd separately
for i in range(0,len(num_center_list1)):
    plt.figure(figsize=(8, 5), dpi=300)
    plt.scatter(np.arange(0,len(dip_mom_bend_1d_128_1[i])), dip_mom_bend_1d_128_1[i], c = 'b', marker = 'o', s = 50)#, label = '$\widetilde \kappa = 10^{-6}$')    
    plt.xlabel("Simulations", fontsize = 20)
    plt.ylabel("$D_{far}$", fontsize = 20)
    plt.xticks(fontsize = 20)
    plt.yticks(fontsize = 20)
    # plt.yscale('log')
    plt.xticks([0,5,10,15,20])
    plt.subplots_adjust(right=0.98,left=0.18,top=0.94,bottom=0.16)
    ax = plt.gca()
    # Increase tick width and length for major and minor ticks
    ax.tick_params(axis='both', which='major', width=3, length=10)  # Major ticks
    ax.tick_params(axis='both', which='minor', width=2, length=5)   # Minor ticks
    # plt.legend(loc = 'best', fontsize = 15)
    plt.savefig(base+plot_folder+"bend_dom_dip_mom_scater_"+ str(L) +"_srand_pooled_" + str(num_center_list1[i]) +'_'+pbond_string+"_"+tol_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+".png")
    plt.clf()
    plt.close()

for i in range(0,len(num_center_list2)):
    plt.figure(figsize=(8, 5), dpi=300)
    plt.scatter(np.arange(0,len(dip_mom_bend_1d_128_2[i])), dip_mom_bend_1d_128_2[i], c = 'b', marker = 'o', s = 50)#, label = '$\widetilde \kappa = 10^{-6}$')    
    plt.xlabel("Simulations", fontsize = 20)
    plt.ylabel("$D_{far}$", fontsize = 20)
    plt.xticks(fontsize = 20)
    plt.yticks(fontsize = 20)
    # plt.yscale('log')
    plt.xlim([0,len(dip_mom_bend_1d_128_2[i])])
    plt.subplots_adjust(right=0.98,left=0.18,top=0.94,bottom=0.16)
    ax = plt.gca()
    # Increase tick width and length for major and minor ticks
    ax.tick_params(axis='both', which='major', width=3, length=10)  # Major ticks
    ax.tick_params(axis='both', which='minor', width=2, length=5)   # Minor ticks
    # plt.legend(loc = 'best', fontsize = 15)
    plt.savefig(base+plot_folder+"bend_dom_dip_mom_scater_"+ str(L) +"_srand_pooled_" + str(num_center_list2[i]) +'_'+pbond_string+"_"+tol_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+".png")
    plt.clf()
    plt.close()

#==========================================================================================================================
# for the L = 64 networks (p < 1)
#==========================================================================================================================
L = 64
srand_list_64 = ['112','113', '114','116','117','118','119','120','121','122']                         # 115 is BAD!!!!!
# srand_list_64 = ['113', '114','116','117','118','119','120','121','122']                         # 115 is BAD!!!!!
num_center_list_64 = np.array([1, 2, 5, 10, 15, 20, 25, 30, 35, 40])
tol64 = 1.0e-7
tol64_str = "%.2e" % tol64 # to write the filename

data_all_64 = []
i = 0
for srand in srand_list_64:
    data_all_64.append([])
    for num_center in num_center_list_64:
        num_dip = num_center*6
        
        ranseed = '667720,601210'    # -- network 1
            
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
    
        fname = base+folder+"mean_dip_mom_bndry_force_"+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol64_str+"_"+kappa_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+".txt"    
        
        # if bending_flag == 0:
        #     fname = base+folder+"all_dom_mean_dip_mom_bndry_force_"+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+".txt"    
        #     plot_folder = plot_folder + str(L)+"_all_dom_"

        data = np.loadtxt(fname)
        
        data_all_64[i].append(data)
    i += 1

data_all_64 = np.array(data_all_64)

dfar64 = np.mean(data_all_64[:,:,1], axis = 0)
dloc64 = np.mean(data_all_64[:,:,5], axis = 0)
yerr64 = data_all_64[:,:,2]/np.sqrt(data_all_64[:,:,0])
yerr64 = np.mean(yerr64, axis = 0)
yerr_dloc64 = data_all_64[:,:,-1]/np.sqrt(data_all_64[:,:,0])
yerr_dloc64 = np.mean(yerr_dloc64, axis = 0)



# plotting the dipole moments vs number of dipoles for all kappa
plt.figure(figsize=(8, 5), dpi=300)
plt.errorbar(num_center_list_64, dfar64, yerr = yerr64, c = 'b', fmt = 'o', ms = 10, label = '$\widetilde \kappa = 10^{-6}$')
# plt.errorbar(num_center_list, dfar_norm, yerr = yerr, c = 'b', fmt = 'o', ms = 10, label = '$\widetilde \kappa = 10^{-6}$')

x_linspace = np.linspace(num_center_list_64[0],num_center_list_64[-1],100) 
line_fit64 = dfar64[0]*x_linspace/num_center_list_64[0]
plt.plot(x_linspace, line_fit64, c = 'b', lw = 3) 

# x_linspace2 = np.linspace(num_center_list[5],num_center_list[-1],100) 
# line_fit2 = dfar[5]*np.power(x_linspace2,2)/num_center_list[5]*0.01
# plt.plot(x_linspace2, line_fit2, c = 'g', lw = 3) 

x_linspace3 = np.linspace(num_center_list_64[5],num_center_list_64[-1],100) 
line_fit3_64 = dfar64[5]*np.power(x_linspace3,3)/num_center_list_64[5]*0.003
plt.plot(x_linspace3, line_fit3_64, c = 'r', lw = 3) 

plt.text(7, 1e-4, '$\sim N_d^1$', fontsize = 20)
plt.text(20, 5e-3, '$\sim N_d^3$', fontsize = 20)
# plt.text(20, 5e-2, '$\sim N_d^2$', fontsize = 20)

plt.xlabel("Number of Dipoles", fontsize = 20)
plt.ylabel("$<D_{far}>$", fontsize = 20)
plt.xticks(fontsize = 20)
plt.yticks(fontsize = 20)
plt.yscale('log')
plt.xscale('log')

plt.subplots_adjust(right=0.98,left=0.18,top=0.98,bottom=0.14)

plt.legend(loc = 'best', fontsize = 15)
plt.savefig(base+plot_folder+"dip_mom_vs_N_L_"+ str(L) +"_srand_pooled_" +pbond_string+"_"+tol_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+".png")

plt.clf()
plt.close()

# Create figure and primary axis
fig = plt.subplots(figsize=(8, 5), dpi=300)
# Plot first dataset on ax1 (left y-axis)
plt.errorbar(num_center_list_64, dfar64, yerr = yerr64,
             c='b', fmt='o', ms=12, alpha=0.7, label='$D_{far}$')
plt.errorbar(num_center_list_64, dloc64, yerr = yerr_dloc64,
             c='r', fmt='*', ms=12, alpha=0.7, label='$D_{loc}$')
ax = plt.gca()
# Increase tick width and length for major and minor ticks
ax.tick_params(axis='both', which='major', width=3, length=10)  # Major ticks
ax.tick_params(axis='both', which='minor', width=2, length=5)   # Minor ticks
plt.xlabel("Number of Dipoles, $N_d$", fontsize = 20)
plt.ylabel("Dipole Moment", fontsize = 20)
# plt.ylim([0,25])
plt.xticks(fontsize = 20)
plt.yticks(fontsize = 20)
plt.subplots_adjust(right=0.98,left=0.15,top=0.98,bottom=0.16)
plt.legend(loc='best', fontsize=18)
# Save and close figure
plt.savefig(base + plot_folder + "dfar_dloc_vs_Nd_L64_p55_srand_pooled_diff_style_" +
            pbond_string + "_" + tol_str + "_" + rlen_txt + "_" +
            mu_str + "_" + mu_c_str + ".png")
plt.clf()
plt.close()


# '''
#==========================================================================================================================
# plotting L = 64 and L = 128 networks together
#==========================================================================================================================
# packing fraction for L = 64
unit_dipole_area = 6 * (np.sqrt(3)/4) * (1.5*1.5)     # 6 * area of unit traingle that makes a unit hexagon (note the side lenght  = 1.5 and not 1)
# unit_dipole_area = np.pi * ((2.5/2)*(2.5/2))     # 2.5 is the actual radius of no overlap
area_of_dipoles_64 = num_center_list_64*unit_dipole_area
# area_of_inner_region_64 = 509.22   # this is gotten by 6* area of a triangle of side = 14. it is approximate total inner area
area_of_inner_region_64 = np.pi * (12*12)   # this is almost circular
area_of_dipoles_128 = num_center_list*unit_dipole_area
area_of_inner_region_128 = np.pi * (24*24)   # this is almost circular

packing_frac64 = area_of_dipoles_64/area_of_inner_region_64
packing_frac128 = area_of_dipoles_128/area_of_inner_region_128

# normalizing dfars
packing_frac128_n4 = packing_frac128[1:]
mean_dfar_n4 = mean_dfar_bend[1:]
std_dfar_n4 = std_dfar_bend[1:]
dfar64_norm = dfar64/dfar64[0]
dfar128_norm = mean_dfar_n4/mean_dfar_n4[0]    # bad cariable name but correct otherwise; this does not have Nd = 1 case
yerr64_norm = yerr64/dfar64[0]
yerr128_norm = std_dfar_n4/mean_dfar_n4[0]

dfar128_norm_all = mean_dfar_bend/mean_dfar_bend[0]    # this has Nd = 1 case as well
yerr128_norm_all = std_dfar_bend/mean_dfar_bend[0]

# plotting the dipole moments vs number of dipoles for all kappa
plt.figure(figsize=(8, 5), dpi=300)
plt.errorbar(packing_frac64, dfar64, yerr = yerr64, c = 'b', fmt = 'o', ms = 10, alpha = 0.7, label = '$R_2 = 25$')
plt.errorbar(packing_frac128_n4, mean_dfar_n4, yerr = std_dfar_n4, c = 'm', fmt = 'd', ms = 10, alpha = 0.7, label = '$R_2 = 50$')

x_linspace = np.linspace(packing_frac64[0],packing_frac64[-1],100) 
line_fit64 = dfar64[0]*x_linspace/packing_frac64[0]
plt.plot(x_linspace, line_fit64, c = 'b', lw = 3) 

# x_linspace3 = np.linspace(packing_frac64[5],packing_frac64[-1],100) 
# line_fit3_64 = dfar64[5]*np.power(x_linspace3,3)/packing_frac64[5]*10
# plt.plot(x_linspace3, line_fit3_64, c = 'r', lw = 3) 

plt.text(0.1, 1e-4, '$\sim \phi^1$', fontsize = 20)
# plt.text(0.38, 1.5e-3, '$\sim \phi^3$', fontsize = 20)

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
plt.savefig(base+plot_folder+"dip_mom_vs_packing_frac_all_L_srand_pooled_" +pbond_string+"_"+tol_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+".png")

plt.clf()
plt.close()


# plotting the dipole moments vs number of dipoles for all kappa
plt.figure(figsize=(8, 5), dpi=300)
plt.errorbar(packing_frac64, dfar64, yerr = yerr64, c = 'b', fmt = 'o', ms = 10, label = '$\widetilde \kappa = 10^{-6}$')

x_linspace = np.linspace(packing_frac64[0],packing_frac64[-1],100) 
line_fit64 = dfar64[0]*x_linspace/packing_frac64[0]
plt.plot(x_linspace, line_fit64, c = 'b', lw = 3) 

x_linspace3 = np.linspace(packing_frac64[5],packing_frac64[-1],100) 
line_fit3_64 = dfar64[5]*np.power(x_linspace3,3)/packing_frac64[5]*30
plt.plot(x_linspace3, line_fit3_64, c = 'r', lw = 3) 

plt.text(0.1, 1e-4, '$\sim \phi^1$', fontsize = 20)
plt.text(0.22, 5e-3, '$\sim \phi^3$', fontsize = 20)

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
plt.savefig(base+plot_folder+"dip_mom_vs_packing_frac_L_"+ str(L) +"_srand_pooled_" +pbond_string+"_"+tol_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+".png")

plt.clf()
plt.close()


# plotting the normalized dipole moments vs number of dipoles for all kappa
plt.figure(figsize=(8, 5), dpi=300)
plt.errorbar(packing_frac64, dfar64_norm, yerr = yerr64_norm, c = 'b', fmt = 'o', ms = 10, alpha = 0.7, label = '$R_2 = 25$')
plt.errorbar(packing_frac128_n4, dfar128_norm, yerr = yerr128_norm, c = 'm', fmt = 'd', ms = 10, alpha = 0.7, label = '$R_2 = 50$')

x_linspace = np.linspace(packing_frac128_n4[0],packing_frac128_n4[-1],100) 
line_fit128 = dfar128_norm[0]*x_linspace/packing_frac128_n4[0]
plt.plot(x_linspace, line_fit128, c = 'b', lw = 3) 

# x_linspace2 = np.linspace(packing_frac128_n4[5],packing_frac128_n4[-1],100) 
# line_fit2 = dfar128_norm[5]*np.power(x_linspace2,2)/packing_frac128_n4[5]*2
# plt.plot(x_linspace2, line_fit2, c = 'g', lw = 3) 

# x_linspace3 = np.linspace(packing_frac128_n4[5],packing_frac128_n4[-1],100) 
# line_fit3_128 = dfar128_norm[5]*np.power(x_linspace3,3)/packing_frac128_n4[5]*10
# plt.plot(x_linspace3, line_fit3_128, c = 'r', lw = 3) 

plt.text(0.08, 4, '$\sim \phi^1$', fontsize = 20)
# plt.text(0.4, 60, '$\sim \phi^2$', fontsize = 20)
# plt.text(0.27, 200, '$\sim \phi^3$', fontsize = 20)

ax = plt.gca()
# Increase tick width and length for major and minor ticks
ax.tick_params(axis='both', which='major', width=3, length=10)  # Major ticks
ax.tick_params(axis='both', which='minor', width=2, length=5)   # Minor ticks

plt.xlabel("Dipole Packing Fraction, $\phi$", fontsize = 20)
plt.ylabel("$<D_{far}>/<D_{far,ref}>$", fontsize = 20)
plt.xticks(fontsize = 20)
plt.yticks(fontsize = 20)
plt.yscale('log')
plt.xscale('log')

plt.subplots_adjust(right=0.98,left=0.15,top=0.98,bottom=0.16)

plt.legend(loc = 'best', fontsize = 15)
plt.savefig(base+plot_folder+"norm_dip_mom_vs_packing_frac_L_"+ str(L) +"_srand_pooled_" +pbond_string+"_"+tol_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+".png")

plt.clf()
plt.close()
# '''


#==========================================================================================================================
# Reading p = 1 L = 128 network to find p_eff
#==========================================================================================================================
fname_base  = '/home/abhinav/david/lattice_nopbd_bndry_all_clamp_restlength_circle_inner_fixed_radial_arp/txt/area/'
fname_plot = '/home/abhinav/david/lattice_nopbd_bndry_all_clamp_restlength_circle_inner_fixed_radial_arp/png/strain_hist/'
fname_plot = fname_plot + str(L)+'_'   # works for L = 128 case only
srand_list = np.concatenate((srand_list1, srand_list2))
L = 128

pbond_string = "%.2f" % 1 # to write the filename
kappa = 1e-6
kappa_str = "%.2e" % kappa # to write the filename

d_far_p1_1 = []
for i in range(0,len(num_center_list1)):
    d_far_p1_1.append([])
    for j in range(0,len(srand_list1)):
        srand = srand_list1[j]
        fname_p1 = fname_base + 'dipole_moment_ratio_'+ 'srand_'+ str(srand) +'_'+ str(L) +'_'+str(num_center_list1[i])+'_'+str(num_center_list1[i]*6)+'_'+pbond_string+ "_"+tol_str+"_"+ kappa_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+'_10.txt'    
        p1_data_temp = np.loadtxt(fname_p1)
        d_far_p1_1[i].append(p1_data_temp[2])

        # if num_center_list[i] == 1:
        #     print(fname_p1)
        #     print(p1_data_temp[2])
    
d_far_p1_2 = []
for i in range(0,len(num_center_list2)):
    d_far_p1_2.append([])
    for j in range(0,len(srand_list2)):
        srand = srand_list2[j]
        fname_p1 = fname_base + 'dipole_moment_ratio_'+ 'srand_'+ str(srand) +'_'+ str(L) +'_'+str(num_center_list2[i])+'_'+str(num_center_list2[i]*6)+'_'+pbond_string+ "_"+tol_str+"_"+ kappa_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+'_10.txt'    
        # print('Fname_p1: ',fname_p1)
        p1_data_temp = np.loadtxt(fname_p1)
        d_far_p1_2[i].append(p1_data_temp[2])

mean_dfar_p1_1 = np.mean(d_far_p1_1, axis=1)
mean_dfar_p1_2 = np.mean(d_far_p1_2, axis=1)
std_dfar_p1_1 = np.std(d_far_p1_1, axis=1)
std_dfar_p1_2 = np.std(d_far_p1_2, axis=1)

mean_dfar_p1 = np.concatenate((mean_dfar_p1_1, mean_dfar_p1_2))
std_dfar_p1 = np.concatenate((std_dfar_p1_1, std_dfar_p1_2))

norm_mean_dfar_p1 = mean_dfar_p1/mean_dfar_p1[0]
norm_std_dfar_p1 = std_dfar_p1/mean_dfar_p1[0]

#==========================================================================================================================
#==========================================================================================================================
# Reading p = 1 L = 64 network to plot with L = 128 network
#==========================================================================================================================
#==========================================================================================================================

# srand_list_64 = ['113', '114','116','117','118','119','120','121','122']                         # 115 is BAD!!!!!
srand_list_64 = ['112', '113', '114','116','117','118','119','120','121','122']                         # 115 is BAD!!!!!
num_center_list_64 = np.array([1, 2, 5, 10, 15, 20, 25, 30, 35, 40])
    
d_far_p1_64 = []
d_loc_p1_64 = []
for i in range(0,len(num_center_list_64)):
    d_far_p1_64.append([])
    d_loc_p1_64.append([])
    for j in range(0,len(srand_list_64)):    
        srand = srand_list_64[j]
        fname_p1 = fname_base + 'dipole_moment_ratio_'+ 'srand_'+ str(srand) +'_'+str(num_center_list_64[i])+'_'+str(num_center_list_64[i]*6)+'_'+pbond_string+ "_"+tol64_str+"_"+ kappa_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+'_10.txt'    
        # print('Fname_p1: ',fname_p1)
        p1_data_temp = np.loadtxt(fname_p1)
        d_far_p1_64[i].append(p1_data_temp[2])
        d_loc_p1_64[i].append(p1_data_temp[1])

mean_dfar_p1_64 = np.mean(d_far_p1_64, axis=1)
std_dfar_p1_64 = np.std(d_far_p1_64, axis=1)

norm_mean_dfar_p1_64 = mean_dfar_p1_64/mean_dfar_p1_64[0]
norm_std_dfar_p1_64 = std_dfar_p1_64/mean_dfar_p1_64[0]

mean_dloc_p1_64 = np.mean(d_loc_p1_64, axis=1)
std_dloc_p1_64 = np.std(d_loc_p1_64, axis=1)

#==========================================================================================================================
# plotting  p = 1: L = 64 and L = 128 cases together 
#==========================================================================================================================
# scatter plot of all p = 1 cases: dfar vs dloc
colors_list = [
    "#1f77b4",  # muted blue
    "#ff7f0e",  # orange
    "#2ca02c",  # green
    "#d62728",  # red
    "#9467bd",  # purple
    "#8c564b",  # brown
    "#e377c2",  # pink
    "#7f7f7f",  # gray
    "#bcbd22",  # yellow-green
    "#17becf"   # cyan
]

markers_list = [
    "o",  # circle
    "s",  # square
    "^",  # triangle_up
    "v",  # triangle_down
    "D",  # diamond
    "X",  # x-filled
    "*",  # star
    "P",  # plus-filled
    "d",  # plus
    "h"   # x
]

plt.figure(figsize=(8, 5), dpi=300)
for i in range(0,len(num_center_list_64)):
    plt.scatter(d_loc_p1_64[i], d_far_p1_64[i], marker = markers_list[i], s = 100, alpha=0.4, c = colors_list[i], label = '$N_d = %d $' % num_center_list_64[i])# Network 1')    

plt.xlabel("$D_{loc}$", fontsize = 20)
plt.ylabel("$D_{far}$", fontsize = 20)
plt.xscale('log')
plt.yscale('log')
plt.xticks(fontsize = 20)
plt.yticks(fontsize = 20)
plt.subplots_adjust(right=0.98,left=0.14,top=0.97,bottom=0.14)
plt.legend(loc = 'best', fontsize = 12)
plt.savefig(plot_folder+"bend_dom_dfar_vs_dloc_scatter_all_srand_all_Nd_"+pbond_string+"_"+kappa_str+"_"+tol_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+".png")
plt.clf()
plt.close()


# plotting the dfar and dloc vs Nd for all kappa
# Create figure and primary axis
fig, ax1 = plt.subplots(figsize=(8, 5), dpi=300)
# Plot first dataset on ax1 (left y-axis)
ax1.errorbar(num_center_list_64, mean_dfar_p1_64, yerr=std_dfar_p1_64,
             c='b', fmt='o', ms=10, alpha=0.7, label='$D_{far}$')
ax1.set_ylabel("$<D_{far}>$", fontsize=20, color='b')
ax1.tick_params(axis='y', labelcolor='b', labelsize=20)
ax1.tick_params(axis='x', labelsize=20)
# Create second y-axis
ax2 = ax1.twinx()
# Plot second dataset on ax2 (right y-axis)
ax2.errorbar(num_center_list_64, mean_dloc_p1_64, yerr=std_dloc_p1_64,
             c='r', fmt='*', ms=10, alpha=0.7, label='$D_{loc}$')
ax2.set_ylabel("$<D_{loc}>$", fontsize=20, color='r')
ax2.tick_params(axis='y', labelcolor='r', labelsize=20)
# Set shared x-label
ax1.set_xlabel("$N_d$", fontsize=20)
# ---------- Synchronize Y-axis scaling ----------
# Get combined y-limits
# ymin = min(ax1.get_ylim()[0], ax2.get_ylim()[0])
# ymax = max(ax1.get_ylim()[1], ax2.get_ylim()[1])
# Apply same limits to both axes
ax1.set_ylim(0, 25)
ax2.set_ylim(0, 25)
# Optionally match ticks too
ax2.set_yticks(ax1.get_yticks())
# ---------- Final formatting ----------
fig.tight_layout()
fig.subplots_adjust(right=0.90, left=0.10, top=0.96, bottom=0.16)
# Combine legends from both plots
lines1, labels1 = ax1.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()
ax1.legend(lines1 + lines2, labels1 + labels2, loc='best', fontsize=15)
# Save and close figure
plt.savefig(base + plot_folder + "dfar_dloc_vs_Nd_L64_srand_pooled_" +
            pbond_string + "_" + tol_str + "_" + rlen_txt + "_" +
            mu_str + "_" + mu_c_str + ".png")
plt.clf()
plt.close()

# Create figure and primary axis
fig = plt.subplots(figsize=(8, 5), dpi=300)
# Plot first dataset on ax1 (left y-axis)
plt.errorbar(num_center_list_64, mean_dfar_p1_64, yerr=std_dfar_p1_64,
             c='b', fmt='o', ms=12, alpha=0.7, label='$D_{far}$')
plt.errorbar(num_center_list_64, mean_dloc_p1_64, yerr=std_dloc_p1_64,
             c='r', fmt='*', ms=12, alpha=0.7, label='$D_{loc}$')
ax = plt.gca()
# Increase tick width and length for major and minor ticks
ax.tick_params(axis='both', which='major', width=3, length=10)  # Major ticks
ax.tick_params(axis='both', which='minor', width=2, length=5)   # Minor ticks
plt.xlabel("Number of Dipoles, $N_d$", fontsize = 20)
plt.ylabel("Dipole Moment", fontsize = 20)
plt.xticks(fontsize = 20)
plt.yticks(fontsize = 20)
plt.subplots_adjust(right=0.98,left=0.15,top=0.98,bottom=0.16)
plt.legend(loc='best', fontsize=18)
# Save and close figure
plt.savefig(base + plot_folder + "dfar_dloc_vs_Nd_L64_srand_pooled_diff_style_" +
            pbond_string + "_" + tol_str + "_" + rlen_txt + "_" +
            mu_str + "_" + mu_c_str + ".png")
plt.clf()
plt.close()


# plotting the dloc of p=1 and dloc of p=0.55 vs number of dipoles for all kappa
plt.figure(figsize=(8, 5), dpi=300)
plt.errorbar(num_center_list_64, mean_dloc_p1_64, yerr = std_dloc_p1_64, c = 'b', fmt = 's', ms = 10, alpha = 0.7, label = '$p = 0.55$')
plt.errorbar(num_center_list_64, dloc64, yerr = yerr_dloc64, c = 'r', fmt = 'd', ms = 10, alpha = 0.7, label = '$p = 1$')

ax = plt.gca()
# Increase tick width and length for major and minor ticks
ax.tick_params(axis='both', which='major', width=3, length=10)  # Major ticks
ax.tick_params(axis='both', which='minor', width=2, length=5)   # Minor ticks

plt.xlabel("$N_d$", fontsize = 20)
plt.ylabel("$<D_{loc}>$", fontsize = 20)
plt.xticks(fontsize = 20)
plt.yticks(fontsize = 20)
plt.subplots_adjust(right=0.98,left=0.15,top=0.98,bottom=0.16)
plt.legend(loc = 'best', fontsize = 15)
plt.savefig(base+plot_folder+"dloc_vs_Nd_L64_L128_srand_pooled_" +pbond_string+"_"+tol_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+".png")
plt.clf()
plt.close()

# plotting the dipole moments vs number of dipoles for all kappa
plt.figure(figsize=(8, 5), dpi=300)
plt.errorbar(num_center_list_64, mean_dfar_p1_64, yerr = std_dfar_p1_64, c = 'b', fmt = 'o', ms = 10, alpha = 0.7, label = '$R_2 = 25$')
plt.errorbar(num_center_list, mean_dfar_p1, yerr = std_dfar_p1, c = 'm', fmt = 'd', ms = 10, alpha = 0.7, label = '$R_2 = 50$')

ax = plt.gca()
# Increase tick width and length for major and minor ticks
ax.tick_params(axis='both', which='major', width=3, length=10)  # Major ticks
ax.tick_params(axis='both', which='minor', width=2, length=5)   # Minor ticks

plt.xlabel("$N_d$", fontsize = 20)
plt.ylabel("$<D_{far}>$", fontsize = 20)
plt.xticks(fontsize = 20)
plt.yticks(fontsize = 20)
plt.yscale('log')
plt.xscale('log')
plt.subplots_adjust(right=0.98,left=0.15,top=0.98,bottom=0.16)
plt.legend(loc = 'best', fontsize = 15)
plt.savefig(base+plot_folder+"dip_mom_vs_Nd_L64_L128_srand_pooled_" +pbond_string+"_"+tol_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+".png")
plt.clf()
plt.close()

# plotting the dipole moments vs packing fraction for all kappa
plt.figure(figsize=(8, 5), dpi=300)
plt.errorbar(packing_frac64, mean_dfar_p1_64, yerr = std_dfar_p1_64, c = 'b', fmt = 'o', ms = 10, alpha = 0.7, label = '$R_2 = 25$')
plt.errorbar(packing_frac128, mean_dfar_p1, yerr = std_dfar_p1, c = 'm', fmt = 'd', ms = 10, alpha = 0.7, label = '$R_2 = 50$')

ax = plt.gca()
# Increase tick width and length for major and minor ticks
ax.tick_params(axis='both', which='major', width=3, length=10)  # Major ticks
ax.tick_params(axis='both', which='minor', width=2, length=5)   # Minor ticks

plt.xlabel("Packing fraction $(\phi)$", fontsize = 20)
plt.ylabel("$<D_{far}>$", fontsize = 20)
plt.xticks(fontsize = 20)
plt.yticks(fontsize = 20)
plt.yscale('log')
plt.xscale('log')
plt.subplots_adjust(right=0.98,left=0.15,top=0.98,bottom=0.16)
plt.legend(loc = 'best', fontsize = 15)
plt.savefig(base+plot_folder+"dip_mom_vs_packing_frac_L64_L128_srand_pooled_" +pbond_string+"_"+tol_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+".png")
plt.clf()
plt.close()

# plotting the normalized dipole moments vs packing fraction for both L
plt.figure(figsize=(8, 5), dpi=300)
plt.errorbar(packing_frac64, norm_mean_dfar_p1_64, yerr = norm_std_dfar_p1_64, c = 'b', fmt = 'o', ms = 10, alpha = 0.7, label = '$R_2 = 25$')
plt.errorbar(packing_frac128, norm_mean_dfar_p1, yerr = norm_std_dfar_p1, c = 'm', fmt = 'd', ms = 10, alpha = 0.7, label = '$R_2 = 50$')

# x_linspace = np.linspace(packing_frac128[0],packing_frac128[-1],100) 
# line_fit128 = dfar128_norm_all[0]*x_linspace/packing_frac128[0]
# plt.plot(x_linspace, line_fit128, c = 'b', lw = 3) 
# plt.text(0.08, 10, '$\sim \phi^1$', fontsize = 20)

# x_linspace = np.linspace(packing_frac128[5],packing_frac128[-1],100) 
# line_fit128 = dfar128_norm_all[5]*np.power(x_linspace,2)/packing_frac128[5]*10
# plt.plot(x_linspace, line_fit128, c = 'r', lw = 3) 
# plt.text(0.2, 1e3, '$\sim \phi^2$', fontsize = 20)

ax = plt.gca()
# Increase tick width and length for major and minor ticks
ax.tick_params(axis='both', which='major', width=3, length=10)  # Major ticks
ax.tick_params(axis='both', which='minor', width=2, length=5)   # Minor ticks

plt.xlabel("Packing Fraction ($\phi$)", fontsize = 20)
plt.ylabel("$<D_{far}>/<D_{far,0}>$", fontsize = 20)
plt.xticks(fontsize = 20)
plt.yticks(fontsize = 20)
plt.yscale('log')
plt.xscale('log')
plt.subplots_adjust(right=0.98,left=0.15,top=0.98,bottom=0.16)
plt.legend(loc = 'best', fontsize = 15)
plt.savefig(base+plot_folder+"dip_mom_norm_vs_packing_frac_L64_L128_srand_pooled_" +pbond_string+"_"+tol_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+".png")
plt.clf()
plt.close()

#==========================================================================================================================
# plotting  p = 0.55 and p = 1 cases together for L = 128
#==========================================================================================================================
# plotting the dipole moments vs number of dipoles for all kappa
plt.figure(figsize=(8, 5), dpi=300)
plt.errorbar(num_center_list, mean_dfar_p1, yerr = std_dfar_p1, c = 'b', fmt = 'o', ms = 10, alpha = 0.7, label = '$p = 1$')
plt.errorbar(num_center_list, mean_dfar_bend, yerr = std_dfar, c = 'm', fmt = 'd', ms = 10, alpha = 0.7, label = '$p = 0.55$')

ax = plt.gca()
# Increase tick width and length for major and minor ticks
ax.tick_params(axis='both', which='major', width=3, length=10)  # Major ticks
ax.tick_params(axis='both', which='minor', width=2, length=5)   # Minor ticks

plt.xlabel("$N_d$", fontsize = 20)
plt.ylabel("$<D_{far}>$", fontsize = 20)
plt.xticks(fontsize = 20)
plt.yticks(fontsize = 20)
plt.yscale('log')
plt.xscale('log')
plt.subplots_adjust(right=0.98,left=0.15,top=0.98,bottom=0.16)
plt.legend(loc = 'best', fontsize = 15)
plt.savefig(base+plot_folder+"dip_mom_vs_Nd_p1_p55_srand_pooled_" +pbond_string+"_"+tol_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+".png")
plt.clf()
plt.close()


# plotting the normalized dipole moments vs number of dipoles for all kappa
plt.figure(figsize=(8, 5), dpi=300)
plt.errorbar(num_center_list, norm_mean_dfar_p1, yerr = norm_std_dfar_p1, c = 'b', fmt = 'o', ms = 10, alpha = 0.7, label = '$p = 1$')
plt.errorbar(num_center_list, dfar128_norm_all, yerr = yerr128_norm_all, c = 'm', fmt = 'd', ms = 10, alpha = 0.7, label = '$p = 0.55$')

x_linspace = np.linspace(num_center_list[0],num_center_list[-1],100) 
line_fit128 = dfar128_norm_all[0]*x_linspace/num_center_list[0]
plt.plot(x_linspace, line_fit128, c = 'b', lw = 3) 
plt.text(40, 10, '$\sim \phi^1$', fontsize = 20)

x_linspace = np.linspace(num_center_list[5],num_center_list[-1],100) 
line_fit128 = dfar128_norm_all[5]*np.power(x_linspace,2)/num_center_list[5]*0.03
plt.plot(x_linspace, line_fit128, c = 'r', lw = 3) 
plt.text(80, 1e3, '$\sim \phi^2$', fontsize = 20)

ax = plt.gca()
# Increase tick width and length for major and minor ticks
ax.tick_params(axis='both', which='major', width=3, length=10)  # Major ticks
ax.tick_params(axis='both', which='minor', width=2, length=5)   # Minor ticks

plt.xlabel("$N_d$", fontsize = 20)
plt.ylabel("$<D_{far}>/<D_{far,0}>$", fontsize = 20)
plt.xticks(fontsize = 20)
plt.yticks(fontsize = 20)
plt.yscale('log')
plt.xscale('log')
plt.subplots_adjust(right=0.98,left=0.15,top=0.98,bottom=0.16)
plt.legend(loc = 'best', fontsize = 15)
plt.savefig(base+plot_folder+"dip_mom_norm_vs_Nd_p1_p55_srand_pooled_" +pbond_string+"_"+tol_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+".png")
plt.clf()
plt.close()


# plotting the normalized dipole moments vs acking frac for all L
plt.figure(figsize=(8, 5), dpi=300)
plt.errorbar(packing_frac128, norm_mean_dfar_p1, yerr = norm_std_dfar_p1, c = 'b', fmt = 'o', ms = 10, alpha = 0.7, label = '$p = 1$')
plt.errorbar(packing_frac128, dfar128_norm_all, yerr = yerr128_norm_all, c = 'm', fmt = 'd', ms = 10, alpha = 0.7, label = '$p = 0.55$')

x_linspace = np.linspace(packing_frac128[0],packing_frac128[-1],100) 
line_fit128 = dfar128_norm_all[0]*x_linspace/packing_frac128[0]
plt.plot(x_linspace, line_fit128, c = 'b', lw = 3) 
plt.text(0.08, 10, '$\sim \phi^1$', fontsize = 20)

x_linspace = np.linspace(packing_frac128[5],packing_frac128[-1],100) 
line_fit128 = dfar128_norm_all[5]*np.power(x_linspace,2)/packing_frac128[5]*10
plt.plot(x_linspace, line_fit128, c = 'r', lw = 3) 
plt.text(0.2, 1e3, '$\sim \phi^2$', fontsize = 20)

ax = plt.gca()
# Increase tick width and length for major and minor ticks
ax.tick_params(axis='both', which='major', width=3, length=10)  # Major ticks
ax.tick_params(axis='both', which='minor', width=2, length=5)   # Minor ticks

plt.xlabel("Packing Fraction ($\phi$)", fontsize = 20)
plt.ylabel("$<D_{far}>/<D_{far,0}>$", fontsize = 20)
plt.xticks(fontsize = 20)
plt.yticks(fontsize = 20)
plt.yscale('log')
plt.xscale('log')
plt.subplots_adjust(right=0.98,left=0.15,top=0.98,bottom=0.16)
plt.legend(loc = 'best', fontsize = 15)
plt.savefig(base+plot_folder+"dip_mom_norm_vs_phi_p1_p55_srand_pooled_" +pbond_string+"_"+tol_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+".png")
plt.clf()
plt.close()
#==========================================================================================================================
# Ratio of dfar values gives us alpha_m or mu_m
#==========================================================================================================================
# the following is mu_m = <dfar_p0.55>/<dfar_p1> :: it also uses all results without looking at energy ratio
mu_m = mean_dfar/mean_dfar_p1
sigma_mu_m = mu_m * np.sqrt(np.square(std_dfar/mean_dfar) + np.square(std_dfar_p1/mean_dfar_p1))

# the following use onnly bending dominted cases
mu_m_bend = mean_dfar_bend/mean_dfar_p1
sigma_mu_m_bend = mu_m_bend * np.sqrt(np.square(std_dfar_bend/mean_dfar_bend) + np.square(std_dfar_p1/mean_dfar_p1))


# plotting the normalized dipole moments vs acking frac for all L
plt.figure(figsize=(8, 5), dpi=300)
plt.errorbar(num_center_list, mu_m_bend, yerr = sigma_mu_m_bend, c = 'g', fmt = 'o', ms = 10, alpha = 0.7, label = '$p = 1$')
# ax = plt.gca()
# Increase tick width and length for major and minor ticks
# ax.tick_params(axis='both', which='major', width=3, length=10)  # Major ticks
# ax.tick_params(axis='both', which='minor', width=2, length=5)   # Minor ticks

plt.xlabel("$N_d$", fontsize = 20)
plt.ylabel("$<D_{far,p=0.55}>/<D_{far,p=1}>$", fontsize = 20)
plt.xticks(fontsize = 20)
plt.yticks(fontsize = 20)
# plt.yscale('log')
# plt.xscale('log')
plt.subplots_adjust(right=0.98,left=0.22,top=0.98,bottom=0.16)
plt.legend(loc = 'best', fontsize = 15)
plt.savefig(base+plot_folder+"mu_m_bend_vs_Nd_p1_p55_srand_pooled_" +pbond_string+"_"+tol_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+".png")
plt.clf()
plt.close()


# the following is mu_m = <dfar_p0.55/dfar_p1>. 
# this is subtly different from what is written above
dfar_ratio1 = []
mean_dfar_per_srand_list = []
for i in range(0,len(num_center_list1)):
    dfar_ratio1.append([])
    mean_dfar_per_srand_list.append([])
    for j in range(0,len(srand_list1)):
        mean_dfar_per_srand = np.mean(data_dip_mom[i][j])     
        mean_dfar_per_srand_list[i].append(mean_dfar_per_srand)
        dfar_ratio_val = mean_dfar_per_srand/d_far_p1_1[i][j]
        dfar_ratio1[i].append(dfar_ratio_val)

mean_dfar_ratio1 = np.mean(dfar_ratio1, axis = 1)
mean_dfar_1 = np.mean(mean_dfar_per_srand_list, axis = 1)

#==========================================================================================================================
# need to use EMT to get p_eff from mu_m
#==========================================================================================================================
from emt_func import compute_p_eff

p_eff, sigma_p = compute_p_eff(mu_m, sigma_mu_m)

p_eff_bend, sigma_p_bend = compute_p_eff(mu_m_bend, sigma_mu_m_bend)  # just bending dominated cases

#==========================================================================================================================
# applying a linear fit and plotting the data
#==========================================================================================================================
from math import log10, floor
def find_exp_base(number):
    exp = floor(log10(abs(number)))
    return round(number/10**exp, 2), exp 

# coefficients = np.polyfit(num_center_list, p_eff_bend, 1, w = 1./sigma_p_bend**2)
coefficients, cov = np.polyfit(num_center_list, p_eff_bend, 1, cov=True)#, w = 1./sigma_p_bend**2)
slope_err, intercept_err = np.sqrt(np.diag(cov))
m, b = coefficients
xaxis_n_d = np.linspace(num_center_list[0], num_center_list[-1], 100)
y_fit_peff = m * xaxis_n_d + b
# Plot the data and the fitted line
base_num, exp = find_exp_base(m)
label_txt = '$' + str(base_num) + '\\times 10^{'+ str(exp) + '}$' + ' $N_d + $' + r"$%0.2f$" % b
print('L = 128 p_eff slope = ', m, '    slope err = ', slope_err)

num_pts_circle_128 = 9271     # total nodes in the circular network
n_c_d_128 = 3*num_pts_circle_128*m           # m[0] is the slope of peff vs Nd; this calculates number of constraints per dipole
n_c_d_err_128 = 3*num_pts_circle_128*slope_err   # error in number of constraints per dipole
print('L = 128:    ', n_c_d_128, n_c_d_err_128)

# bending dominated networks: mu_m vs Nd for L = 128
plt.figure(figsize=(8, 5), dpi=300)
plt.errorbar(num_center_list, mu_m_bend, yerr = sigma_mu_m_bend, c = 'g', ls = 'none', marker = 'o', alpha = 1, ms = 10, label = '$R_2 = 50$')
plt.xlabel("Number of dipoles, $N_d$", fontsize = 20)
plt.ylabel("$D_{far,p=0.55}/D_{far,p=1}$", fontsize = 20)
plt.xticks(fontsize = 20)
plt.yticks(fontsize = 20)
formatter = FuncFormatter(scientific_format)
plt.gca().yaxis.set_major_formatter(formatter)
plt.subplots_adjust(right=0.98,left=0.24,top=0.97,bottom=0.15)
plt.legend(loc = 'upper left', fontsize = 15)
plt.savefig(base+plot_folder+"mu_m_vs_N_L128_" +pbond_string+"_"+tol_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+".png")
plt.clf()
plt.close()

# bending dominated networks: peff vs Nd for L = 128
plt.figure(figsize=(8, 5), dpi=300)
plt.errorbar(num_center_list, p_eff_bend, yerr = sigma_p_bend, c = 'g', ls = 'none', marker = 'o', alpha = 1, ms = 10, label = '$R_2 = 50$')
plt.plot(xaxis_n_d, y_fit_peff, c = 'b', ls = 'dashed', lw = 4, label = label_txt)
plt.xlabel("Number of dipoles, $N_d$", fontsize = 20)
plt.ylabel("$p_{eff}$", fontsize = 20)
plt.xticks(fontsize = 20)
plt.yticks(fontsize = 20)
plt.subplots_adjust(right=0.97,left=0.14,top=0.97,bottom=0.14)
plt.legend(loc = 'best', fontsize = 15)
plt.savefig(base+plot_folder+"peff_vs_N_L128_" +pbond_string+"_"+tol_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+".png")
plt.clf()
plt.close()

# bending dominated and stretching dominated networks: peff vs Nd for L = 128
coefficients = np.polyfit(num_center_list, p_eff, 1)#, w = 1./sigma_p_bend**2)
m, b = coefficients
xaxis_n_d = np.linspace(num_center_list[0], num_center_list[-1], 100)
y_fit_peff = m * xaxis_n_d + b
# Plot the data and the fitted line
base_num, exp = find_exp_base(m)
label_txt = '$' + str(base_num) + '\\times 10^{'+ str(exp) + '}$' + ' $N_d + $' + r"$%0.2f$" % b

plt.figure(figsize=(8, 5), dpi=300)
plt.errorbar(num_center_list, p_eff, yerr = sigma_p, c = 'g', ls = 'none', marker = 'o', alpha = 1, ms = 10, label = '$R_2 = 50$')
plt.plot(xaxis_n_d, y_fit_peff, c = 'b', ls = 'dashed', lw = 4, label = label_txt)
plt.xlabel("$N_d$", fontsize = 15)
plt.ylabel("$p_{eff}$", fontsize = 15)
plt.xticks(fontsize = 15)
plt.yticks(fontsize = 15)

plt.subplots_adjust(right=0.96,left=0.18,top=0.94,bottom=0.12)
plt.legend(loc = 'best', fontsize = 15)
plt.savefig(base+plot_folder+"peff_vs_N_L128_all_dom_" +pbond_string+"_"+tol_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+".png")
plt.clf()
plt.close()

#==========================================================================================================================
# redoing old master plot of p_eff vs N_d
# the old method in network_circular_caller.py averages all the p_eff values for each srand - which were in turn written by the same code
# this method actually reads in the alpha_m or mu_m values and averages them across different srand before calculating peff.
# this works only for L=64 for srand = 113 through 122
#==========================================================================================================================
emt_folder = "cluster_download/shear_modulus_EMT/"  # this is where the plot is saved

num_center_val = 40    # now we can go upto 40 although srand 114 and 122 do not work anymore
pbond_string_val = '1.00'
# srand_list = ['113', '114','116','117','118','119','120','121','122']                         # 115 is BAD!!!!!
srand_list = ['112','113', '114','116','117','118','119','120','121','122']                         # 115 is BAD!!!!!

mu_m_64_dfar = []
mu_m_sigma_64_dfar = []
for i in range(0,len(srand_list)):
    srand = srand_list[i]
    disp_folder = base + 'cluster_download/lattice_nopbd_bndry_all_clamp_restlength_circle_inner_fixed_radial_arp_bash_radial/'+ kappa_fname +ranseed+"/"+str(srand)+"/means/"
    fname_alpha = disp_folder +  "alpha_m_values_"+str(srand)+"_"+str(num_center_val)+"_"+str(num_center_val*6)+"_"+pbond_string_val+"_"+tol64_str+"_"+kappa_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+".txt"
    # print(fname_alpha)
    alpha_all_data = np.loadtxt(fname_alpha)
    x_new = alpha_all_data[:,1]
    sigma_x_new = alpha_all_data[:,2]
    
    mu_m_64_dfar.append(x_new)
    mu_m_sigma_64_dfar.append(sigma_x_new)
    
mean_mu_m_64_dfar = np.mean(mu_m_64_dfar, axis = 0)
mean_mu_m_sigma_64_dfar = np.mean(mu_m_sigma_64_dfar, axis = 0)

# finding p_eff based on the mu_m values
p_eff_64, sigma_p_64 = compute_p_eff(mean_mu_m_64_dfar, mean_mu_m_sigma_64_dfar)


# plotting
coefficients = np.polyfit(num_center_list_64, p_eff_64, 1, w = 1./sigma_p_64, cov=True)
m, b = coefficients
y_fit_peff = m[0] * num_center_list_64 + m[1]
# Plot the data and the fitted line
base_num, exp = find_exp_base(m[0])
label_txt = '$' + str(base_num) + '\\times 10^{'+ str(exp) + '}$' + ' $N_d + $' + r"$%0.2f$" % m[1]
slope_error = np.sqrt(b[0, 0])

    
fig = plt.figure(figsize=(8, 6), dpi=300)
plt.errorbar(num_center_list_64, p_eff_64, yerr = sigma_p_64, c = 'g', ls = 'none', marker = 'o', alpha = 1, ms = 10, label = '$p = 0.55$   $\widetilde\kappa = 10^{-6}$ ')
plt.plot(num_center_list_64, y_fit_peff, label=label_txt)
plt.xlabel('Number of dipoles, $ N_d $', fontsize = 25)
# plt.ylabel(' $ \\Sigma1 / \\alpha_{m}$ ', fontsize = 25)
plt.ylabel(' $p_{eff}$ ', fontsize = 25)
plt.xticks(fontsize = 20)
plt.yticks(fontsize = 20)
# plt.yscale('log')
# plt.xscale('log')
plt.subplots_adjust(top=0.96, left = 0.17, bottom = 0.15, right = 0.98)
plt.legend(loc = 'best', fontsize = 15)
plt.savefig(base+emt_folder+"peff_srand_113_through_122_pooled_mu_m_first_tol_"+tol64_str+".png")  # just plotting the last step data
plt.close()
plt.clf()


fig = plt.figure(figsize=(8, 6), dpi=300)
plt.errorbar(num_center_list_64, mean_mu_m_64_dfar, yerr = mean_mu_m_sigma_64_dfar, c = 'g', ls = 'none', marker = 'o', alpha = 1, ms = 10, label = '$p = 0.55$   $\widetilde\kappa = 10^{-6}$ ')
plt.xlabel('Number of dipoles, $N_d $', fontsize = 25)
# plt.ylabel(' $ \\Sigma1 / \\alpha_{m}$ ', fontsize = 25)
plt.ylabel(' $D_{far, p = 0.55}/D_{far, p = 1}$ ', fontsize = 25)
plt.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
plt.xticks(fontsize = 20)
plt.yticks(fontsize = 20)
# plt.yscale('log')
# plt.xscale('log')
formatter = FuncFormatter(scientific_format)
plt.gca().yaxis.set_major_formatter(formatter)
plt.subplots_adjust(top=0.96, left = 0.25, bottom = 0.15, right = 0.98)
plt.legend(loc = 'best', fontsize = 15)
plt.savefig(base+emt_folder+"mu_m_srand_113_through_122_pooled_mu_m_first_tol_"+tol64_str+".png")  # just plotting the last step data
plt.close()
plt.clf()


###############################################################################
# number of constraints per dipole
###############################################################################
# it works for L = 64!!
num_pts_circle = 2347     # total nodes in the circular network
n_c_d = 3*num_pts_circle*m[0]           # m[0] is the slope of peff vs Nd; this calculates number of constraints per dipole
n_c_d_err = 3*num_pts_circle*slope_err   # error in number of constraints per dipole

print('L = 64:    ', n_c_d, n_c_d_err)
