#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 13 16:19:15 2025
the code reads in the mean values of fitted sigma1/mu_m  (from the radial displacement plots)
it also reads in the mean values of measured sigma1 (from the stress in rings calculations)

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

# num_center_list = [5, 10, 15]
num_center_list = [1, 2, 5, 10, 15, 20, 25, 30, 35, 40]
# num_center_list = [1, 5, 10, 15, 20]    # for paper plot of dfar vs kappa
# num_center_list = [1, 2, 5, 10]

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

kappa = 1e-6
# kappa_list = [1e-6, 1e-5, 1e-4]
# kappa_list = [1e-6]
# if num_center_list[-1] == 20:
#     kappa_list = [1e-6, 5e-6, 1e-5, 5e-5, 1e-4]

rlen = 0.9
rlen_txt = "%.4f" % rlen
rlen_txt_ring = "%.9f" % 0.1

srand = 113

pbond_string = "%.2f" % pbond # to write the filename
tol_str = "%.2e" % tol # to write the filename
    
mu_str = "%.4f" % mu # to write the filename
mu_c_str = "%.4f" % mu_c # to write the filename

base = "/home/abhinav/david/"

cluster_new_flag = 1   # for dipoles with only radial bonds
hex_flag = 0           # for dipoles with hexagonal bonds forced to be present
hex_rand_flag = 0      # for dipoles with hexagonal bonds randomly removed as per p value
cluster_radial_flag = 1   # for dipoles with just radial bonds

data_sigma1_fitted = []
sigma1_stress = []
sigma1_fitted_stress = []

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
        plot_folder = 'cluster_download/lattice_nopbd_bndry_all_clamp_restlength_circle_inner_fixed_radial_arp_bash/'+ kappa_fname +ranseed+"/"+str(srand)+"/"
    elif cluster_new_flag == 1 and hex_flag == 1 and hex_rand_flag == 0 and cluster_radial_flag == 0:
        folder = 'cluster_download/lattice_nopbd_bndry_all_clamp_restlength_circle_inner_fixed_radial_arp_bash_hex/'+ kappa_fname +ranseed+"/"+str(srand)+"/means/"
        plot_folder = 'cluster_download/lattice_nopbd_bndry_all_clamp_restlength_circle_inner_fixed_radial_arp_bash_hex/'+ kappa_fname +ranseed+"/"+str(srand)+"/"
    elif cluster_new_flag == 1 and hex_flag == 0 and hex_rand_flag == 1 and cluster_radial_flag == 0:
        folder = 'cluster_download/lattice_nopbd_bndry_all_clamp_restlength_circle_inner_fixed_radial_arp_bash_random_dipole_hex_bonds/'+ kappa_fname +ranseed+"/"+str(srand)+"/means/"
        plot_folder = 'cluster_download/lattice_nopbd_bndry_all_clamp_restlength_circle_inner_fixed_radial_arp_bash_random_dipole_hex_bonds/'+ kappa_fname +ranseed+"/"+str(srand)+"/"
    elif cluster_new_flag == 1 and hex_flag == 0 and hex_rand_flag == 0 and cluster_radial_flag == 1:
        folder = 'cluster_download/lattice_nopbd_bndry_all_clamp_restlength_circle_inner_fixed_radial_arp_bash_radial/'+ kappa_fname +ranseed+"/"+str(srand)+"/means/"
        plot_folder = 'cluster_download/lattice_nopbd_bndry_all_clamp_restlength_circle_inner_fixed_radial_arp_bash_radial/'+ kappa_fname +ranseed+"/"+str(srand)+"/"

    # reading the fits to displacements which gives us sigma1/mu_m
    disp_fname = base+folder+"fit_radial_disp_"+str(srand)+"_"+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+".txt"
    data = np.loadtxt(disp_fname)    
    data_sigma1_fitted.append(data)    # its actually sigma1/mu_m

    # reading the measured sigma1 values (read fom rings)
    stress_fname = base+folder+"mean_ring_stress_"+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+".txt"    
    data = np.loadtxt(stress_fname)    
    sigma1_stress.append(data)    

    # reading the fitted sigma1 values (read from fits to stress_rr values)
    stress_fitted_fname = base+folder+"mean_sigma1_got_from_fit_to_radial_stress_"+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+".txt"    
    data = np.loadtxt(stress_fitted_fname)    
    sigma1_fitted_stress.append(data)    # its actually sigma1/mu_m

data_sigma1_fitted = np.array(data_sigma1_fitted)   # its actually sigma1/mu_m
data_sigma1_divided_by_mu_m = data_sigma1_fitted[:,-1]

sigma1_stress = np.array(sigma1_stress)   # its measured from the stress in rings concentric to the center
# sigma1_inner_stress = sigma1_stress[:,12,1]     # 12th row and 1st column has the stress (counting starts at 0)
sigma1_inner_stress = sigma1_stress[:,12,1]     # 12th row and 1st column has the stress (counting starts at 0)

sigma1_fitted_stress = np.array(sigma1_fitted_stress)
sigma1_fitted_stress = sigma1_fitted_stress[:,1]

# finding mu_m from measured stress on the rings (using one ring here)
mu_m = -1*(sigma1_inner_stress/data_sigma1_divided_by_mu_m)

# finding mu_m from sigma1 that we got by fitting radial stress to theory
mu_m3 = -1*(sigma1_fitted_stress/data_sigma1_divided_by_mu_m)

# plotting the mu_m vs N_d
plt.figure(figsize=(8, 5), dpi=300)
plt.scatter(num_center_list, mu_m, marker = 'o', s = 150, label = '$\widetilde \kappa = 10^{-6}$')
plt.xlabel("Number of Dipoles", fontsize = 15)
plt.ylabel("$\mu_m$", fontsize = 15)
plt.xticks(num_center_list, fontsize = 15)
plt.yticks(fontsize = 15)
# plt.yscale('log')
# plt.xscale('log')
plt.subplots_adjust(right=0.96,left=0.18,top=0.94,bottom=0.12)
plt.legend(loc = 'best')
plt.savefig(base+plot_folder+"mu_m_vs_dipole_num_kappa_1e-6_"+ str(srand) +"_" +pbond_string+"_"+tol_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+".png")
plt.clf()
plt.close()


# plotting the mu_m vs N_d
plt.figure(figsize=(8, 5), dpi=300)
plt.scatter(num_center_list, sigma1_inner_stress, marker = 'o', s = 150, label = '$\widetilde \kappa = 10^{-6}$')
plt.xlabel("Number of Dipoles", fontsize = 15)
plt.ylabel("$\Sigma_{1}$", fontsize = 15)
plt.xticks(num_center_list, fontsize = 15)
plt.yticks(fontsize = 15)
# plt.yscale('log')
# plt.xscale('log')
plt.subplots_adjust(right=0.96,left=0.18,top=0.94,bottom=0.12)
plt.legend(loc = 'best')
plt.savefig(base+plot_folder+"sigma_m_vs_dipole_num_kappa_1e-6_"+ str(srand) +"_" +pbond_string+"_"+tol_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+".png")
plt.clf()
plt.close()

# plotting the simga1 (from fits to radial stress vs radius) vs N_d
plt.figure(figsize=(8, 5), dpi=300)
plt.scatter(num_center_list, sigma1_fitted_stress, marker = 'o', s = 150, label = '$\widetilde \kappa = 10^{-6}$')
plt.xlabel("Number of Dipoles", fontsize = 15)
plt.ylabel("$\Sigma_{1}$", fontsize = 15)
plt.xticks(num_center_list, fontsize = 15)
plt.yticks(fontsize = 15)
# plt.yscale('log')
# plt.xscale('log')
plt.subplots_adjust(right=0.96,left=0.18,top=0.94,bottom=0.12)
plt.legend(loc = 'best')
plt.savefig(base+plot_folder+"sigma1_from_stress_fit_vs_dipole_num_kappa_1e-6_"+ str(srand) +"_" +pbond_string+"_"+tol_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+".png")
plt.clf()
plt.close()

# plotting the mu_m got from simga1 (from fits to radial stress vs radius) 
plt.figure(figsize=(8, 5), dpi=300)
plt.scatter(num_center_list, mu_m3, marker = 'o', s = 150, label = '$\widetilde \kappa = 10^{-6}$')
plt.xlabel("Number of Dipoles", fontsize = 15)
plt.ylabel("$\mu_m$", fontsize = 15)
plt.xticks(num_center_list, fontsize = 15)
plt.yticks(fontsize = 15)
# plt.yscale('log')
# plt.xscale('log')
plt.subplots_adjust(right=0.96,left=0.18,top=0.94,bottom=0.12)
plt.legend(loc = 'best')
plt.savefig(base+plot_folder+"mu_m_from_sigma1_from_stress_fit_vs_dipole_num_kappa_1e-6_"+ str(srand) +"_" +pbond_string+"_"+tol_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+".png")
plt.clf()
plt.close()

#########################################################################################################################
# saving the data!!
#########################################################################################################################
# writing out the file with sigma1 measured from stress in rings
outfname = base + folder + "simga1_measured_rings_values_"+pbond_string+"_"+tol_str+"_"+kappa_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+".txt"
print('here the sigma1/mu_m values found by fitting displacements are written to file : ',outfname)    
heading = '# Dipoles          Measured Sigma1'
fmt = '%12d', '%20.7e'
# np.savetxt(outfname, np.column_stack((num_center_list[i], mu_m, np.sqrt(np.diag(pcov))[0])), header = heading, fmt = fmt)
np.savetxt(outfname, np.column_stack((num_center_list, mu_m)), header = heading, fmt = fmt)


# writing out the file with sigma1 fitted from radial stress vs radial distance plots
outfname = base + folder + "simga1_fitted_radial_distance_values_"+pbond_string+"_"+tol_str+"_"+kappa_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+".txt"
print('here the sigma1/mu_m values found by fitting displacements are written to file : ',outfname)    
heading = '# Dipoles          Measured Sigma1'
fmt = '%12d', '%20.7e'
# np.savetxt(outfname, np.column_stack((num_center_list[i], mu_m, np.sqrt(np.diag(pcov))[0])), header = heading, fmt = fmt)
np.savetxt(outfname, np.column_stack((num_center_list, mu_m3)), header = heading, fmt = fmt)

'''
#########################################################################################################################
# finding the mean sigma1 divided by mu_m for all srand
#########################################################################################################################
# new_srand_list = ['113', '114', '116','117','118','119','120','121', '122']                         # 115 is BAD!!!!!
new_srand_list = ['113', '116', '117', '118', '119', '120', '121', '122']                         # 115 is BAD!!!!!
# new_srand_list = ['112','113','114', '116','117','118','119','120','121', '122']                         # 115 is BAD!!!!!
data_all_stress = []
data_all_stress_fit = []   # fits to stress gives sigma1
data_all_stress_rings = []   # fits to stress gives sigma1
data_all_disp = []
data_all_mu_m_from_dfar = []

for i in np.arange(0,len(num_center_list)):
    data_all_stress.append([])
    data_all_stress_fit.append([])
    data_all_stress_rings.append([])
    data_all_disp.append([])
    for srand_val in new_srand_list:    
        # reading fitted displacements to get sigma1/mu_m
        folder = 'cluster_download/lattice_nopbd_bndry_all_clamp_restlength_circle_inner_fixed_radial_arp_bash_radial/'+ kappa_fname +ranseed+"/"+str(srand_val)+"/means/"
        stress_fname = base+folder+"fit_radial_disp_"+str(srand_val)+"_"+str(num_center_list[i])+"_"+str(num_center_list[i]*6)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+".txt"    
        data_stress = np.loadtxt(stress_fname)
        data_all_stress[i].append(data_stress)
        
        # reading fitted displacements to get sigma1/mu_m from each mean for each srand value (the above file just looks at each mean fit value and is happy with that. here we will ool all the data together and then find the fit)
        disp_fname = base+folder+"bndry_node_radial_disp_mean_"+str(srand_val)+"_"+str(num_center_list[i])+"_"+str(num_center_list[i]*6)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+".txt"    
        data_disp = np.loadtxt(disp_fname)
        data_all_disp[i].append(data_disp)

        # reading the fitted sigma1 values (read from fits to stress_rr values)
        stress_fitted_fname = base+folder+"mean_sigma1_got_from_fit_to_radial_stress_"+str(num_center_list[i])+"_"+str(num_center_list[i]*6)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+".txt"    
        data = np.loadtxt(stress_fitted_fname)    
        data_all_stress_fit[i].append(data)    # its actually sigma1/mu_m

        # reading the measured sigma1 values (read from radial rings)
        stress_rings_fname = base+folder+"mean_ring_stress_"+str(num_center_list[i])+"_"+str(num_center_list[i]*6)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+".txt"    
        data = np.loadtxt(stress_rings_fname)    
        data_all_stress_rings[i].append(data[12,1])    # its actually sigma1/mu_m
        
for srand_val in new_srand_list:    
    folder = 'cluster_download/lattice_nopbd_bndry_all_clamp_restlength_circle_inner_fixed_radial_arp_bash_radial/'+ kappa_fname +ranseed+"/"+str(srand_val)+"/means/"
    mu_m_from_dfar_fname = base + folder + "alpha_m_values_"+str(srand_val)+"_"+str(num_center_list[-1])+"_"+str(num_center_list[-1]*6)+"_"+'1.00'+"_"+tol_str+"_"+kappa_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+".txt"
    data = np.loadtxt(mu_m_from_dfar_fname)    
    data_all_mu_m_from_dfar.append(data[:,1])    # its actually sigma1/mu_m


plot_folder = 'cluster_download/lattice_nopbd_bndry_all_clamp_restlength_circle_inner_fixed_radial_arp_bash_radial/'+ kappa_fname +ranseed+"/"

# for data from displacement vs r fits: gives us sigma1/mu_m
data_all_stress = np.array(data_all_stress)
data_all_stress = data_all_stress[:,:,-1]
mean_stresses_fit = np.mean(data_all_stress, axis = 1)

# for data from stress vs r fits: gives us sigma1
data_all_stress_fit = np.array(data_all_stress_fit)
data_all_stress_fit = data_all_stress_fit[:,:,-1]
mean_stresses_fit_from_stress = np.mean(data_all_stress_fit, axis = 1)

# for data from radial displacements averaged over multiple outer network configurations
data_all_disp = np.array(data_all_disp)
data_disp = data_all_disp[:,:,:,1]
disp_bins = data_all_disp[0,0,:,0]
mean_radial_disp = np.mean(data_disp, axis = 1)
std_radial_disp = np.std(data_disp, axis = 1)

# for data from measured stress  in rings: gives us sigma1 at the right radial distance
data_all_stress_rings = np.array(data_all_stress_rings)
# data_all_stress_rings = data_all_stress_rings[:,:,-1]
mean_stresses_from_rings = np.mean(data_all_stress_rings, axis = 1)
mu_m2 = -1*(mean_stresses_from_rings/mean_stresses_fit)  # from stress on rings
mu_m_stress_fits = -1*(mean_stresses_fit_from_stress/mean_stresses_fit)   # from stress vs radial fits


# for data from dfar rations of p =/055 and p=1: gives us mu_m
data_all_mu_m_from_dfar = np.array(data_all_mu_m_from_dfar)
mean_mu_m_from_dfar = np.mean(data_all_mu_m_from_dfar, axis = 0)


fig = plt.figure(figsize=(8, 6), dpi=300)
plt.scatter(num_center_list, mu_m2, s = 150, label='random dipole positions: measured stress')
plt.scatter(num_center_list, mu_m_stress_fits, s = 150, label='random dipole positions: fitted stress')
plt.xlabel('$ N_d $', fontsize = 25)
plt.ylabel(' $\\mu_m$ ', fontsize = 25)
plt.xticks(fontsize = 25)
plt.yticks(fontsize = 25)
# plt.yscale('log')
# plt.xscale('log')
plt.subplots_adjust(top=0.96, left = 0.17, bottom = 0.15, right = 0.98)
plt.legend(loc = 'best', fontsize = 15)
plt.savefig(base+plot_folder+"mu_m_srand_113_through_122_pooled.png")  # just plotting the last step data
plt.close()
plt.clf()


fig = plt.figure(figsize=(8, 6), dpi=300)
plt.scatter(num_center_list, mean_stresses_from_rings, s = 150, label='random dipole positions: measured stress')
# plt.scatter(num_center_list, mu_m_stress_fits, s = 150, label='random dipole positions: fitted stress')
plt.xlabel('$ N_d $', fontsize = 25)
plt.ylabel(' $\\Sigma_1$ ', fontsize = 25)
plt.xticks(fontsize = 25)
plt.yticks(fontsize = 25)
plt.yscale('log')
plt.xscale('log')
plt.subplots_adjust(top=0.96, left = 0.17, bottom = 0.15, right = 0.98)
plt.legend(loc = 'best', fontsize = 15)
plt.savefig(base+plot_folder+"sigma1_srand_113_through_122_pooled.png")  # just plotting the last step data
plt.close()
plt.clf()

fig = plt.figure(figsize=(8, 6), dpi=300)
plt.scatter(num_center_list[:6], mean_stresses_from_rings[:6], s = 150, label='random dipole positions: measured stress')
# plt.scatter(num_center_list, mu_m_stress_fits, s = 150, label='random dipole positions: fitted stress')
plt.xlabel('$ N_d $', fontsize = 25)
plt.ylabel(' $\\Sigma_1$ ', fontsize = 25)
plt.xticks(fontsize = 25)
plt.yticks(fontsize = 25)
# plt.yscale('log')
# plt.xscale('log')
plt.subplots_adjust(top=0.96, left = 0.17, bottom = 0.15, right = 0.98)
plt.legend(loc = 'best', fontsize = 15)
plt.savefig(base+plot_folder+"sigma1_srand_113_through_122_pooled_first6Nd.png")  # just plotting the last step data
plt.close()
plt.clf()

fig = plt.figure(figsize=(8, 6), dpi=300)
plt.scatter(num_center_list[-5:], mean_stresses_from_rings[-5:], s = 150, label='random dipole positions: measured stress')
# plt.scatter(num_center_list, mu_m_stress_fits, s = 150, label='random dipole positions: fitted stress')
plt.xlabel('$ N_d $', fontsize = 25)
plt.ylabel(' $\\Sigma_1$ ', fontsize = 25)
plt.xticks(fontsize = 25)
plt.yticks(fontsize = 25)
plt.yscale('log')
plt.xscale('log')
plt.subplots_adjust(top=0.96, left = 0.17, bottom = 0.15, right = 0.98)
plt.legend(loc = 'best', fontsize = 15)
plt.savefig(base+plot_folder+"sigma1_srand_113_through_122_pooled_last5Nd.png")  # just plotting the last step data
plt.close()
plt.clf()


#=============================================================================================================
# fitting displacements
#=============================================================================================================
from scipy.optimize import curve_fit
r1 = 13
r2 = 25
fit_theory_outer_xdata = np.linspace(r1, r2)

def func(x,a):
    return (a) * (r1**2) / (2*((np.sqrt(3)/2) * r1**2 + ((np.sqrt(3)/4) * r2**2)))  * ( (r2**2 - x**2)/x )

from math import log10, floor
def find_exp_base(number):
    exp = floor(log10(abs(number)))
    return round(number/10**exp, 2), exp 

EMT_folder = 'cluster_download/shear_modulus_EMT/'

# plotting the mean radial displacement vs radial distance
all_sigma1_by_mu_m = []
all_sigma1_by_mu_m_error = []

for i in range(0,len(num_center_list)):
    # finding the fit
    popt, pcov = curve_fit(func, disp_bins[r1-1:], mean_radial_disp[i,r1-1:], maxfev=10000)  # -1 beacuse counting starts from 0 
    # print('Fit parameters: ', popt)
    fit_theory_outer_ydata = func(fit_theory_outer_xdata, *popt)
    base_num, exp = find_exp_base(popt[0])
    base_num_err, exp_err = find_exp_base(np.sqrt(np.diag(pcov))[0])
    label_txt = '$\\Sigma_{1}/\\mu_{m} =' + str(base_num) + ' \\times 10^{'+ str(exp) + '}$' + ' $ \\pm '+ str(base_num_err) + '\\times 10^{'+ str(exp_err) + '}$'

    # plotting each of them
    fig = plt.figure(figsize=(8, 6), dpi=300)
    # plt.scatter(disp_bins[12:], mean_radial_disp[i,12:], s = 150, label='$N_d = %d $' % num_center_list[i])
    plt.errorbar(disp_bins[12:], mean_radial_disp[i,12:], yerr = std_radial_disp[i,12:], ms = 10, ls = 'None', marker = 'o', label='$N_d = %d $' % num_center_list[i])
    plt.plot(fit_theory_outer_xdata, fit_theory_outer_ydata, 'b-', label = label_txt)    
    # plt.scatter(num_center_list, mu_m_stress_fits, s = 150, label='random dipole positions: fitted stress')
    plt.xlabel('Radial Distance', fontsize = 25)
    plt.ylabel(' $ \\langle u_r \\rangle $ ', fontsize = 25)
    plt.xticks(fontsize = 25)
    plt.yticks(fontsize = 25)
    # plt.yscale('log')
    # plt.xscale('log')
    plt.subplots_adjust(top=0.96, left = 0.23, bottom = 0.15, right = 0.98)
    plt.legend(loc = 'best', fontsize = 15)
    plt.savefig(base+plot_folder+"radial_disp_srand_113_through_122_pooled_N_"+str(num_center_list[i])+".png")  # just plotting the last step data
    plt.close()
    plt.clf()

    # writing out the file with sigma1/mu_m fits from fits to radial displacement
    outfname = base + EMT_folder + "simga1_by_mu_m_values_srand113_122_pooled_"+str(num_center_list[i])+"_"+str(num_center_list[i])+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+".txt"
    print('here the sigma1/mu_m values found by fitting displacements are written to file : ',outfname)    
    heading = '# Dipoles         Sigma1/mu_m         error'
    fmt = '%12d', '%20.7e', '%20.7e'
    np.savetxt(outfname, np.column_stack((num_center_list[i], popt[0], np.sqrt(np.diag(pcov))[0])), header = heading, fmt = fmt)
    
    # saving the data in an array
    all_sigma1_by_mu_m.append(popt[0])
    all_sigma1_by_mu_m_error.append(np.sqrt(np.diag(pcov))[0])

all_sigma1_by_mu_m = np.array(all_sigma1_by_mu_m)
all_sigma1_by_mu_m_error = np.array(all_sigma1_by_mu_m_error)

#============================================================================================================================================================
#============================================================================================================================================================
# reading the data from EMT data thief file
#============================================================================================================================================================
#============================================================================================================================================================
emt_folder = "cluster_download/shear_modulus_EMT/"

fname_kappa_e_6 = base+emt_folder+"kappae-6_shear_modulus.txt"    
data_e_6 = np.loadtxt(fname_kappa_e_6, delimiter=',')
shear_e_6_emt = data_e_6[:,1]

fname_kappa_e_5 = base+emt_folder+"kappae-5_shear_modulus_manual.txt"    
data_e_5 = np.loadtxt(fname_kappa_e_5, delimiter=',')
shear_e_5_emt = data_e_5[:,1]

fname_kappa_e_4 = base+emt_folder+"kappae-4_shear_modulus_manual.txt"    
data_e_4 = np.loadtxt(fname_kappa_e_4, delimiter=',')
shear_e_4_emt = data_e_4[:,1]

# finding mu_m from shear modulus read from EMT file
alpha_m_list_e_6_emt = 4*shear_e_6_emt/np.sqrt(3)
alpha_m_list_e_5_emt = 4*shear_e_5_emt/np.sqrt(3)
alpha_m_list_e_4_emt = 4*shear_e_4_emt/np.sqrt(3)

# finding the bulk modulus from mu_m values - dont know why
bulk_m_list_e_6_emt = 0.5*alpha_m_list_e_6_emt*np.sqrt(3)
bulk_m_list_e_5_emt = 0.5*alpha_m_list_e_5_emt*np.sqrt(3)
bulk_m_list_e_4_emt = 0.5*alpha_m_list_e_4_emt*np.sqrt(3)

x = alpha_m_list_e_6_emt    # x is the alpha_m values
y = data_e_6[:,0]   # using only kappa e-6 for now; y is p value of the network

# Create the cubic spline interpolator
from scipy.interpolate import CubicSpline
f = CubicSpline(x, y, bc_type='natural')  # 'natural' for natural boundary conditions

# Interpolate y values using the spline
mu_m2 = -1*(mean_stresses_from_rings/all_sigma1_by_mu_m)  # from stress on rings
mu_m_stress_fits = -1*(mean_stresses_fit_from_stress/all_sigma1_by_mu_m)   # from stress vs radial fits

peff_stress_on_rings = f(mu_m2)
dy_dx_stress_on_rings = f(mu_m2, 1)                    # the derivative of the function
# sigma_peff_stress_on_rings = np.abs(dy_dx_stress_on_rings) * sigma_mu_m2  # the error in y

peff_stress_fits = f(mu_m_stress_fits)
dy_dx_stress_fits = f(mu_m_stress_fits, 1)                    # the derivative of the function
# sigma_peff_stress_on_rings = np.abs(dy_dx_stress_on_rings) * sigma_mu_m2  # the error in y

# plotting peff vs dipole number for the stress in rings case
fig = plt.figure(figsize=(8, 6), dpi=300)
plt.scatter(num_center_list, peff_stress_on_rings, s = 150, label='Stress on Rings')
plt.scatter(num_center_list, peff_stress_fits, s = 150, label='Stress fits')
plt.scatter(num_center_list, f(2*mu_m2), s = 150, label='test')
plt.xlabel('$ N_d $', fontsize = 25)
plt.ylabel(' $ P_{eff} $ ', fontsize = 25)
plt.xticks(fontsize = 25)
plt.yticks(fontsize = 25)
# plt.yscale('log')
# plt.xscale('log')
plt.subplots_adjust(top=0.96, left = 0.17, bottom = 0.15, right = 0.98)
plt.legend(loc = 'best', fontsize = 15)
plt.savefig(base+EMT_folder+"peff_vs_mu_m_stressFits_and_stress_rings.png")  # just plotting the last step data
plt.close()
plt.clf()

# plotting EMT p vs mu_m
fig = plt.figure(figsize=(8, 6), dpi=300)
plt.scatter(x, y, s = 150, label='EMT')
# plt.scatter(num_center_list, mu_m_stress_fits, s = 150, label='random dipole positions: fitted stress')
plt.xlabel('$ \\mu_m $', fontsize = 25)
plt.ylabel(' $ P_{eff} $ ', fontsize = 25)
plt.xticks(fontsize = 25)
plt.yticks(fontsize = 25)
plt.yscale('log')
plt.xscale('log')
plt.subplots_adjust(top=0.96, left = 0.17, bottom = 0.15, right = 0.98)
plt.legend(loc = 'best', fontsize = 15)
plt.savefig(base+EMT_folder+"emt_p_vs_mu_m.png")  # just plotting the last step data
plt.close()
plt.clf()

'''