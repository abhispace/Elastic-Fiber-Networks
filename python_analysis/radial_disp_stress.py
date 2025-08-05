#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May  2 22:25:43 2025

reads the mean radial displacement, fits them and extracts sigma1/mu_m
the mean radial displacement is written by all_srand_all_nd_diff_network.py

for p=1, it reads the results from mean_radial_disp_p1.py

it also reads all the individual stresses in concentric rings for each simulation
and fits the decay as a function of radial distance to get sigma1

using sigma1 and sigma1/mu_m (from displacement results), it calculates the p_eff

@author: abhinav
"""

import numpy as np
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
from matplotlib import cm
import matplotlib.pyplot as plt
import math
from scipy.optimize import curve_fit
from math import log10, floor

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

ranseed = '667720,601210'    # -- network 1

#==================================================================================================================================
# reading files*
#format of dfar  : heading = 'N         mean dfar         median dfar       std dfar         sem dfar'
#format of energy: heading = '        N          #sim   mean el en      median el en    std el en       mean bend en    median bend en   std bend en    mean en ratio   median en ratio  std en ratio'

folder1 = 'cluster_download/lattice_nopbd_bndry_all_clamp_restlength_circle_inner_fixed_radial_arp_bash_radial/'+ kappa_fname +ranseed+"/means/"
plot_folder1 = 'cluster_download/lattice_nopbd_bndry_all_clamp_restlength_circle_inner_fixed_radial_arp_bash_radial/'+ kappa_fname +ranseed+"/plots/"

radial_disp_all = []
radial_stress_all = []

for i in range(0,len(num_center_list)):
    num_center = num_center_list[i]
    radial_disp_fname = base + folder1 + "radial_disp_srand_113_122_N_"+str(num_center)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+".txt"
    # heading = 'radius         mean disp         median disp       std disp         sem disp'
    radial_disp_temp = np.loadtxt(radial_disp_fname)    
    radial_disp_all.append(radial_disp_temp)
    
    # reading radial stresses
    stress_fname = base + folder1 + "radial_stress_srand_113_122_N_"+str(num_center)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+".txt"
    # heading = 'radius         mean stress         median stress       std stress         sem stress'
    radial_stress_temp = np.loadtxt(stress_fname)
    radial_stress_all.append(radial_stress_temp)

radial_disp_all = np.array(radial_disp_all)
radial_bins = radial_disp_all[0,:,0]

radial_stress_all = np.array(radial_stress_all)

# fitting the data
# def func(x):
#     return ( (x * r1**2) / (2*(bulk_m * r1**2+ shear_m * r2**2)) ) * ( (r2**2 - fit_x_outer**2)/fit_x_outer )
def func(x,a):
    return (a) * (r1**2) / (2*((np.sqrt(3)/2) * r1**2 + ((np.sqrt(3)/4) * r2**2)))  * ( (r2**2 - x**2)/x )

def find_exp_base(number):
    exp = floor(log10(abs(number)))
    return round(number/10**exp, 2), exp 

r1 = 13
r2 = 25

sigma1_mu_m = []                 # stores the value of sigma1_mu_m extracted from displacement profiles
sigma1_mu_m_err = []                 # stores the value of sigma1_mu_m extracted from displacement profiles

for i in range(0,len(num_center_list)):
    num_center = num_center_list[i]
    
    mean_radial_disp = radial_disp_all[i,:,1]
    median_radial_disp = radial_disp_all[i,:,2]
    std_radial_disp = radial_disp_all[i,:,3]
    sem_radial_disp = radial_disp_all[i,:,4]

    xdata = radial_bins[r1-1:]
    ydata = mean_radial_disp[r1-1:]
    # yerr = sem_radial_disp[r1-1:]
    yerr = std_radial_disp[r1-1:]

    popt, pcov = curve_fit(func, xdata, ydata, maxfev=10000)
    fit_error = np.sqrt(np.diag(pcov))[0]    
    # print('Fit parameters: ', popt, 'Fit error: ', fit_error)
    sigma1_mu_m.append(popt[0])
    sigma1_mu_m_err.append(fit_error)
    
    fit_theory_outer_xdata = np.linspace(xdata[0], xdata[-1])
    fit_theory_outer_ydata = func(fit_theory_outer_xdata, *popt)
    base_num, exp = find_exp_base(popt[0])
    base_num_err, exp_err = find_exp_base(fit_error)
    label_txt = '$\\Sigma_{1}/\\mu_{m} =' + str(base_num) + ' \\times 10^{'+ str(exp) + '}$' + ' $ \\pm '+ str(base_num_err) + '\\times 10^{'+ str(exp_err) + '}$'
    # print(base_num, )

    fig = plt.figure(figsize=(8, 4), dpi=300)
    # plt.scatter(bins, mean_disp)
    plt.errorbar(xdata, ydata, yerr = yerr, marker = 'o', ls = 'none')
    plt.plot(fit_theory_outer_xdata, fit_theory_outer_ydata, 'b-', label = label_txt)    
    plt.xlabel('Distance from center', fontsize = 15)
    plt.ylabel('Radial Displacement', fontsize = 15)
    plt.xticks(fontsize = 15)
    plt.yticks(fontsize = 15)
    plt.ylim([np.min(ydata),0])   # for p = 1; N = 1

    # plt.xlim([xdata[0]+0.5,xdata[-1]+0.5])
    plt.subplots_adjust(top=0.92, left = 0.17, bottom = 0.15, right = 0.98)
    plt.legend(loc = 'best')

    plt.savefig(plot_folder1+"displacement_scatter_mean_outer_"+str(num_center)+"_"+pbond_string+"_"+mu_str+"_"+mu_c_str+"_"+kappa_str+"_"+tol_str+"_"+rlen_txt+"_"+str(num-1)+".png")  # just plotting the last step data
    plt.close()
    plt.clf()


sigma1_mu_m_err = np.array(sigma1_mu_m_err)
sigma1_mu_m = np.array(sigma1_mu_m)

#============================================================================================
# writing the sigma1_mu_m vs N to an output file
#============================================================================================
outfname = base + folder1 + "sigma1_mu_m_vs_N_values_srand_113_122_"+pbond_string+"_"+tol_str+"_"+kappa_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+".txt"
print('writing the sigma1_mu_m vs N to an output file : ',outfname)    
heading = 'N_d           sigma1/mu_m       error'
fmt = '%12d', '%15.7e', '%15.7e'
np.savetxt(outfname, np.column_stack((num_center_list, sigma1_mu_m, sigma1_mu_m_err)), header = heading, fmt = fmt)

#============================================================================================
# plotting the sigma1_mu_m vs N
#============================================================================================
fig = plt.figure(figsize=(8, 6), dpi=300)
plt.errorbar(num_center_list, -1*sigma1_mu_m, yerr= sigma1_mu_m_err, ls = "None", marker = 'o', ms = 10, c = 'b', label = '$p = 0.55$')
plt.xlabel('$N_d$', fontsize = 20)
plt.ylabel('$\Sigma_1/\mu_m$', fontsize = 20)
plt.xticks(fontsize = 20)
plt.yticks(fontsize = 20)
plt.subplots_adjust(top=0.96, left = 0.17, bottom = 0.15, right = 0.98)
plt.legend(loc='best', fontsize = 15)
plt.savefig(base + plot_folder1 + "sigma1_mum_m_vs_N_scatter_outer_region_mean_"+'srand_113_122_'+pbond_string+"_"+mu_str+"_"+mu_c_str+"_"+kappa_str+"_"+tol_str+"_"+rlen_txt+"_"+str(num-1)+".png")  # just plotting the last step data

plt.close()
plt.clf()


# ###############################################################################    
# reading different cases of radial displacement to plot on the same plot for a figure panel for the paper
# ###############################################################################    
temp_num_center1 = 5
temp_num_center2 = 10
temp_num_center3 = 15

# r1 = 13
# xdata = bins[r1:]
fit_x_outer = np.linspace(xdata[0], r2, 100)     # outer region x range
bins_outer = np.arange(xdata[0],r2+1)

radial_data1 = radial_disp_all[2,r1-1:,1]
radial_data1_std = radial_disp_all[2,r1-1:,3]

popt1, pcov1 = curve_fit(func, xdata, radial_data1, maxfev=10000)
# popt1, pcov1 = curve_fit(func, xdata, radial_data1, sigma=radial_data1_std, maxfev=10000)
u_r_outer_fit1 = func(fit_x_outer, popt1[0])

radial_data2 = radial_disp_all[3,r1-1:,1]
radial_data2_std = radial_disp_all[3,r1-1:,3]

popt2, pcov2 = curve_fit(func, xdata, radial_data2, maxfev=10000)
# popt2, pcov2 = curve_fit(func, xdata, radial_data2, sigma=radial_data2_std, maxfev=10000)
u_r_outer_fit2 = func(fit_x_outer, popt2[0])

radial_data3 = radial_disp_all[4,r1-1:,1]
radial_data3_std = radial_disp_all[4,r1-1:,3]

popt3, pcov3 = curve_fit(func, xdata, radial_data3, maxfev=10000)
u_r_outer_fit3 = func(fit_x_outer, popt3[0])

base_num, exp = find_exp_base(popt1[0])
label1_txt = '$\\Sigma_{1}/\\mu_{m} =' + str(base_num) + ' \\times 10^{'+ str(exp) + '}$'
base_num, exp = find_exp_base(popt2[0])
label2_txt = '$\\Sigma_{1}/\\mu_{m} =' + str(base_num) + ' \\times 10^{'+ str(exp) + '}$'
base_num, exp = find_exp_base(popt3[0])
label3_txt = '$\\Sigma_{1}/\\mu_{m} =' + str(base_num) + ' \\times 10^{'+ str(exp) + '}$'


# plotting them together
fig = plt.figure(figsize=(8, 6), dpi=300)
plt.errorbar(bins_outer, radial_data1, yerr= radial_data1_std, ls = "None", marker = 'o', ms = 10, c = 'b', label = '$N_d = 5$')
plt.errorbar(bins_outer, radial_data2, yerr= radial_data2_std, ls = "None", marker = 'd', ms = 10, c = 'orange', label = '$N_d = 10$')
plt.errorbar(bins_outer, radial_data3, yerr= radial_data3_std, ls = "None", marker = 's', ms = 10, c = 'g', label = '$N_d = 15$')
plt.plot(fit_x_outer, u_r_outer_fit1, c = 'b', label = label1_txt)
plt.plot(fit_x_outer, u_r_outer_fit2, c = 'orange', label = label2_txt)
plt.plot(fit_x_outer, u_r_outer_fit3, c = 'g', label = label3_txt)
plt.xlabel('Radial Distance', fontsize = 20)
plt.ylabel('Radial Displacement', fontsize = 20)
plt.xticks(fontsize = 20)
plt.yticks(fontsize = 20)

plt.subplots_adjust(top=0.96, left = 0.19, bottom = 0.15, right = 0.98)
plt.legend(loc='best', fontsize = 15)
plt.savefig(base + plot_folder1 + "mean_displacement_scatter_outer_region_mean_"+'srand_113_122_' +str(temp_num_center1)+"_"+str(temp_num_center2)+"_"+str(temp_num_center3)+"_"+pbond_string+"_"+mu_str+"_"+mu_c_str+"_"+kappa_str+"_"+tol_str+"_"+rlen_txt+"_"+str(num-1)+".png")  # just plotting the last step data

plt.close()
plt.clf()



###############################################################################
# reading p=1 results
###############################################################################
p1_folder = 'lattice_nopbd_bndry_all_clamp_restlength_circle_inner_fixed_radial_arp/'
p1_data = []
for i in range(0,len(num_center_list)):
    num_center = num_center_list[i]
    fname = base + p1_folder + "txt/displacement/radial_disp_srand_113_122_N_"+str(num_center)+"_1.00_"+tol_str+"_"+kappa_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+".txt"
    # heading = 'radius         mean disp         median disp       std disp         sem disp'
    p1_data.append(np.loadtxt(fname))

p1_data = np.array(p1_data)

# ###############################################################################    
# finding the radial displacement for p=1 networks
# ###############################################################################    
sigma1_mu_m_p1 = []                 # stores the value of sigma1_mu_m extracted from displacement profiles
sigma1_mu_m_err_p1 = []                 # stores the value of sigma1_mu_m extracted from displacement profiles

for i in range(0,len(num_center_list)):
    num_center = num_center_list[i]
    
    mean_radial_disp_p1 = p1_data[i,:,1]
    median_radial_disp_p1 = p1_data[i,:,2]
    std_radial_disp_p1 = p1_data[i,:,3]
    sem_radial_disp_p1 = p1_data[i,:,4]

    ydata_p1 = mean_radial_disp_p1[r1-1:]
    yerr_p1 = std_radial_disp_p1[r1-1:]

    popt_p1, pcov_p1 = curve_fit(func, xdata, ydata_p1, maxfev=10000)
    fit_error_p1 = np.sqrt(np.diag(pcov_p1))[0]    
    sigma1_mu_m_p1.append(popt_p1[0])
    sigma1_mu_m_err_p1.append(fit_error_p1)
    
    fit_theory_outer_ydata_p1 = func(fit_theory_outer_xdata, *popt_p1)
    base_num_p1, exp_p1 = find_exp_base(popt_p1[0])
    base_num_err_p1, exp_err_p1 = find_exp_base(fit_error_p1)
    label_txt_p1 = '$\\Sigma_{1}/\\mu_{m} =' + str(base_num_p1) + ' \\times 10^{'+ str(exp_p1) + '}$' + ' $ \\pm '+ str(base_num_err_p1) + '\\times 10^{'+ str(exp_err_p1) + '}$'

    fig = plt.figure(figsize=(8, 4), dpi=300)
    plt.errorbar(xdata, ydata_p1, yerr = yerr_p1, marker = 'o', ls = 'none')
    plt.plot(fit_theory_outer_xdata, fit_theory_outer_ydata_p1, 'b-', label = label_txt_p1)    
    plt.xlabel('Distance from center', fontsize = 15)
    plt.ylabel('Radial Displacement', fontsize = 15)
    plt.xticks(fontsize = 15)
    plt.yticks(fontsize = 15)
    plt.ylim([np.min(ydata_p1),0])   # for p = 1; N = 1

    plt.subplots_adjust(top=0.92, left = 0.17, bottom = 0.15, right = 0.98)
    plt.legend(loc = 'best')

    plt.savefig(plot_folder1+"displacement_scatter_mean_outer_p1_"+str(num_center)+"_"+pbond_string+"_"+mu_str+"_"+mu_c_str+"_"+kappa_str+"_"+tol_str+"_"+rlen_txt+"_"+str(num-1)+".png")  # just plotting the last step data
    plt.close()
    plt.clf()


sigma1_mu_m_err_p1 = np.array(sigma1_mu_m_err_p1)
sigma1_mu_m_p1 = np.array(sigma1_mu_m_p1)

#============================================================================================
# plotting the sigma1_mu_m vs N: p = 0.55 (kappa == 1e-6) and p = 1 together
#============================================================================================
fig = plt.figure(figsize=(8, 6), dpi=300)
plt.errorbar(num_center_list, -1*sigma1_mu_m, yerr= sigma1_mu_m_err, ls = "None", marker = 'o', ms = 10, c = 'b', label = '$p = 0.55$')
plt.errorbar(num_center_list, -1*sigma1_mu_m_p1, yerr= sigma1_mu_m_err_p1, ls = "None", marker = 's', ms = 10, c = 'r', label = '$p = 1$')
plt.xlabel('$N_d$', fontsize = 20)
plt.ylabel('$\Sigma_1/\mu_m$', fontsize = 20)
plt.xticks(fontsize = 20)
plt.yticks(fontsize = 20)
plt.subplots_adjust(top=0.96, left = 0.17, bottom = 0.15, right = 0.98)
plt.legend(loc='best', fontsize = 15)
plt.savefig(base + plot_folder1 + "sigma1_mum_m_vs_N_scatter_outer_region_mean_p1_55_"+'srand_113_122_'+pbond_string+"_"+mu_str+"_"+mu_c_str+"_"+kappa_str+"_"+tol_str+"_"+rlen_txt+"_"+str(num-1)+".png")  # just plotting the last step data

plt.close()
plt.clf()


# normalized plot
sigma1_mu_m_norm = sigma1_mu_m/sigma1_mu_m[0]
sigma1_mu_m_err_norm = sigma1_mu_m_err/abs(sigma1_mu_m[0])
sigma1_mu_m_p1_norm = sigma1_mu_m_p1/sigma1_mu_m_p1[0]
sigma1_mu_m_err_p1_norm = sigma1_mu_m_err_p1/abs(sigma1_mu_m_p1[0])
xfit = np.linspace(num_center_list[0], num_center_list[-1], 100)
line_fit = sigma1_mu_m_p1_norm[0]*xfit/xfit[0]

fig = plt.figure(figsize=(8, 6), dpi=300)
plt.errorbar(num_center_list, sigma1_mu_m_p1_norm, yerr= sigma1_mu_m_err_p1_norm, ls = "None", marker = 'o', ms = 12, c = 'r', alpha = 0.5, label = '$p = 1$')
plt.errorbar(num_center_list, sigma1_mu_m_norm, yerr= sigma1_mu_m_err_norm, ls = "None", marker = 'd', ms = 12, c = 'b', alpha = 0.5, label = '$p = 0.55$')
plt.plot(xfit, line_fit, 'r-')    
plt.text(27, 23, '$\sim N_d^1$', fontsize = 20)
plt.xlabel('$N_d$', fontsize = 20)
plt.ylabel('Normalized $\Sigma_1/\mu_m$', fontsize = 20)
plt.xticks(fontsize = 20)
plt.yticks(fontsize = 20)
plt.subplots_adjust(top=0.96, left = 0.12, bottom = 0.15, right = 0.98)
plt.legend(loc='best', fontsize = 15)
plt.savefig(base + plot_folder1 + "sigma1_mum_m_vs_N_scatter_outer_region_mean_p1_55_norm_"+'srand_113_122_'+pbond_string+"_"+mu_str+"_"+mu_c_str+"_"+kappa_str+"_"+tol_str+"_"+rlen_txt+"_"+str(num-1)+".png")  # just plotting the last step data

plt.close()
plt.clf()








#========================================================================================================================================================================================
# Fitting radial stresses now
#========================================================================================================================================================================================
def func(x,sigma1_fit):
    sigma_rr_theory_fit = -1*(sigma1_fit)*((-outer_radius**2/(2*x**2))-1) * (inner_radius**2)/(inner_radius**2 + (1/2)*outer_radius**2)   # maybe sign error; multiplying by -1
    return sigma_rr_theory_fit

sigma1_list = []                 # stores the value of sigma1_list extracted from displacement profiles
sigma1_list_err = []                 # stores the value of sigma1_list extracted from displacement profiles
cutoff_index = r1-1

sigma1_stress = []  # stores sigma1 using direct stress calculation

for i in range(0,len(num_center_list)):
    num_center = num_center_list[i]
    
    mean_radial_stress = radial_stress_all[i,:,1]
    median_radial_stress = radial_stress_all[i,:,2]
    std_radial_stress = radial_stress_all[i,:,3]
    sem_radial_stress = radial_stress_all[i,:,4]

    xdata = radial_bins[r1-1:]
    ydata = mean_radial_stress[r1-1:]
    # yerr = sem_radial_stress[r1-1:]
    yerr = sem_radial_stress[r1-1:]

    # sigma_rr from theory (using measured sigma1)
    sigma1 = mean_radial_stress[cutoff_index]
    r = np.linspace(xdata[0], xdata[-1],100)
    # outer_radius = 24
    # inner_radius = 13
    outer_radius = r2
    inner_radius = r1
    const = (sigma1*inner_radius**2)/(2*((inner_radius**2*(np.sqrt(3)/2))  + (outer_radius**2*(np.sqrt(3)/4))))
    # sigma_rr_theory = -1*( (-2*(np.sqrt(3)/4)*const*(outer_radius**2)/r**2) - (2*(np.sqrt(3)/2)*const) )   # maybe sign error; multiplying by -1
    
    # getting sigma1 as a fit parameter from the plot data (because mu_m cancels out)    
    popt, pcov = curve_fit(func, xdata, ydata, maxfev=10000)
    # popt, pcov = curve_fit(func, ring_stress[0,cutoff_index:-2,0], mean_ring_stress[cutoff_index:-2], maxfev=10000)
    sigma_rr_outer_fit = func(r, popt[0])
    fit_error = np.sqrt(np.diag(pcov))[0]    
    # print('Fit parameters: ', popt, 'Fit error: ', fit_error)
    sigma1_list.append(popt[0])
    sigma1_list_err.append(fit_error)
    
    fit_theory_outer_xdata = np.linspace(xdata[0], xdata[-1])
    fit_theory_outer_ydata = func(fit_theory_outer_xdata, *popt)
    base_num, exp = find_exp_base(popt[0])
    base_num_err, exp_err = find_exp_base(fit_error)
    label_txt = '$\\Sigma_{1} =' + str(base_num) + ' \\times 10^{'+ str(exp) + '}$' + ' $ \\pm '+ str(base_num_err) + '\\times 10^{'+ str(exp_err) + '}$'
    # print(base_num, )

    fig = plt.figure(figsize=(8, 4), dpi=300)
    # plt.scatter(bins, mean_stress)
    plt.errorbar(xdata, ydata, yerr = yerr, marker = 'o', ls = 'none')
    plt.plot(fit_theory_outer_xdata, fit_theory_outer_ydata, 'b-', label = label_txt)    
    plt.xlabel('Distance from center', fontsize = 15)
    plt.ylabel('Radial stress', fontsize = 15)
    plt.xticks(fontsize = 15)
    plt.yticks(fontsize = 15)
    # plt.ylim([np.min(ydata),0])   # for p = 1; N = 1

    # plt.xlim([xdata[0]+0.5,xdata[-1]+0.5])
    plt.subplots_adjust(top=0.92, left = 0.17, bottom = 0.15, right = 0.98)
    plt.legend(loc = 'best')

    plt.savefig(plot_folder1+"stress_scatter_mean_outer_"+str(num_center)+"_"+pbond_string+"_"+mu_str+"_"+mu_c_str+"_"+kappa_str+"_"+tol_str+"_"+rlen_txt+"_"+str(num-1)+".png")  # just plotting the last step data
    plt.close()
    plt.clf()
    
    # sgima1 directly from radial stress calculation (no need to fit anything)
    sigma1_stress.append(mean_radial_stress[12])

sigma1_list_err = np.array(sigma1_list_err)
sigma1_list = np.array(sigma1_list)

sigma1_stress = np.array(sigma1_stress)

#============================================================================================
# writing the sigma1_list vs N to an output file
#============================================================================================
outfname = base + folder1 + "sigma1_radial_stress_vs_N_values_srand_113_122_"+pbond_string+"_"+tol_str+"_"+kappa_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+".txt"
print('writing the sigma1_list vs N to an output file : ',outfname)    
heading = 'N_d           sigma1/mu_m       error'
fmt = '%12d', '%15.7e', '%15.7e'
np.savetxt(outfname, np.column_stack((num_center_list, sigma1_list, sigma1_list_err)), header = heading, fmt = fmt)

#============================================================================================
# plotting the sigma1_list vs N
#============================================================================================
fig = plt.figure(figsize=(8, 6), dpi=300)
plt.errorbar(num_center_list, sigma1_list, yerr= sigma1_list_err, ls = "None", marker = 'o', ms = 10, c = 'b', label = '$p = 0.55$')
plt.xlabel('$N_d$', fontsize = 20)
plt.ylabel('$\Sigma_1/\mu_m$', fontsize = 20)
plt.xticks(fontsize = 20)
plt.yticks(fontsize = 20)
plt.subplots_adjust(top=0.96, left = 0.17, bottom = 0.15, right = 0.98)
plt.legend(loc='best', fontsize = 15)
plt.savefig(base + plot_folder1 + "sigma1_radial_stress_vs_N_scatter_outer_region_mean_"+'srand_113_122_'+pbond_string+"_"+mu_str+"_"+mu_c_str+"_"+kappa_str+"_"+tol_str+"_"+rlen_txt+"_"+str(num-1)+".png")  # just plotting the last step data

plt.close()
plt.clf()

#==========================================================================================================================
# need to use EMT to get p_eff from mu_m
#==========================================================================================================================
from emt_func import compute_p_eff

mu_m = abs(sigma1_list/sigma1_mu_m)
sigma_mu_m = mu_m * np.sqrt(np.square(sigma1_list_err/sigma1_list) + np.square(sigma1_mu_m_err/sigma1_mu_m))
p_eff, sigma_p = compute_p_eff(mu_m, sigma_mu_m)

# for direct stress measurements
mu_m_stress = abs(sigma1_stress/sigma1_mu_m)
# sigma_mu_m = mu_m * np.sqrt(np.square(sigma1_list_err/sigma1_list) + np.square(sigma1_mu_m_err/sigma1_mu_m))
p_eff_stress, sigma_p_stress = compute_p_eff(mu_m_stress)

#==========================================================================================================================
# applying a linear fit and plotting the data
#==========================================================================================================================
# xaxis_n_d = np.linspace(num_center_list[0], num_center_list[-1], 100)
# coefficients = np.polyfit(num_center_list, p_eff, 1)#, w = 1./peff_err)

xaxis_n_d = np.linspace(num_center_list[4], num_center_list[-1], 100)
coefficients = np.polyfit(num_center_list[4:], p_eff[4:], 1)#, w = 1./peff_err)

m, b = coefficients
y_fit_peff = m * xaxis_n_d + b
# Plot the data and the fitted line
base_num, exp = find_exp_base(m)
label_txt = '$' + str(base_num) + '\\times 10^{'+ str(exp) + '}$' + ' $N_d + $' + r"$%0.2f$" % b


plt.figure(figsize=(8, 5), dpi=300)
plt.errorbar(num_center_list, p_eff, yerr = sigma_p, c = 'g', ls = 'none', marker = 'o', alpha = 1, ms = 10, label = '$p = 0.55$ $\kappa = 10^{-6}$ ')
plt.plot(xaxis_n_d, y_fit_peff, c = 'b', ls = 'dashed', lw = 4, label = label_txt)
plt.xlabel("$N_d$", fontsize = 15)
plt.ylabel("$p_{eff}$", fontsize = 15)
plt.xticks(fontsize = 15)
plt.yticks(fontsize = 15)

plt.subplots_adjust(right=0.96,left=0.18,top=0.94,bottom=0.12)
plt.legend(loc = 'best', fontsize = 15)
plt.savefig(base+plot_folder1+"peff_vs_N_radial_disp_radial_stress_" +pbond_string+"_"+tol_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+".png")
plt.clf()
plt.close()


#=====================================================================================================================
# Using the mu_m from the ratio of far field dipole moments with the sigma1/mu_m from the radial fits
#=====================================================================================================================
# mu_m from dfar ratios of p=0.55 and p=1 networks. this is taken from dipole_moment_means_srand_all_srand.py
mu_m_dfar = [3.72898213e-05, 3.01852967e-05, 3.64765343e-05, 4.19297665e-05,
       5.20026038e-05, 7.34576867e-05, 1.03803623e-04, 1.48035573e-04,
       2.22384919e-04, 3.61419353e-04]

# the error in mu_m; from the same code as mu_m
mu_m_dfar_sigma = [1.03575734e-05, 5.66122897e-06, 4.44677442e-06, 4.85623104e-06,
       6.85203332e-06, 1.31719066e-05, 2.91723175e-05, 5.54558167e-05,
       1.71678456e-04, 2.63994947e-04]


# sigma1_mu_m
# sigma1_mu_m_err

sigma1_dfar_rad_disp = sigma1_mu_m*mu_m_dfar
# sigma1_dfar_rad_disp_err = 

from matplotlib.ticker import FuncFormatter
# Define the formatter function
def scientific_format(x, pos):
    if x == 0:
        return "0"
    exponent = int(np.floor(np.log10(abs(x))))
    coefficient = x / 10**exponent
    return r"${:.1f} \times 10^{{{}}}$".format(coefficient, exponent)

# plotting sigma1 vs Nd
plt.figure(figsize=(8, 5), dpi=300)
# plt.errorbar(num_center_list, sigma1_dfar_rad_disp, yerr = sigma1_dfar_rad_disp_err, c = 'g', ls = 'none', marker = 'o', alpha = 1, ms = 10, label = '$p = 0.55$ $\kappa = 10^{-6}$ ')
plt.scatter(num_center_list, sigma1_dfar_rad_disp, c = 'g', marker = 'o', alpha = 1, s = 100, label = '$p = 0.55$ $\kappa = 10^{-6}$ ')
plt.xlabel("Number of dipoles, $N_d$", fontsize = 20)
plt.ylabel("$\Sigma_{1}$", fontsize = 20)
plt.xticks(fontsize = 15)
plt.yticks(fontsize = 15)

formatter = FuncFormatter(scientific_format)
plt.gca().yaxis.set_major_formatter(formatter)
plt.subplots_adjust(right=0.98,left=0.24,top=0.97,bottom=0.15)
plt.legend(loc = 'best', fontsize = 15)
plt.savefig(base+plot_folder1+"sigma1_vs_N_radial_disp_radial_stress_and_dfar_used_" +pbond_string+"_"+tol_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+".png")
plt.clf()
plt.close()


# plotting sigma1 vs Nd on a log scale
plt.figure(figsize=(8, 5), dpi=300)
plt.scatter(num_center_list, abs(sigma1_dfar_rad_disp), c = 'g', marker = 'o', alpha = 1, s = 100, label = '$p = 0.55$ $\kappa = 10^{-6}$ ')
plt.xlabel("Number of dipoles, $N_d$", fontsize = 20)
plt.ylabel("|$\Sigma_{1}$|", fontsize = 20)
plt.xticks(fontsize = 15)
plt.yticks(fontsize = 15)
plt.yscale('log')
plt.xscale('log')

formatter = FuncFormatter(scientific_format)
plt.gca().yaxis.set_major_formatter(formatter)
plt.subplots_adjust(right=0.98,left=0.24,top=0.97,bottom=0.15)
plt.legend(loc = 'best', fontsize = 15)
plt.savefig(base+plot_folder1+"sigma1_vs_N_radial_disp_radial_stress_and_dfar_used_log_scale_" +pbond_string+"_"+tol_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+".png")
plt.clf()
plt.close()

