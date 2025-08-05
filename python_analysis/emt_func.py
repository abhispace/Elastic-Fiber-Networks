#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May  4 08:13:26 2025

@author: abhinav
"""
import numpy as np
import matplotlib.pyplot as plt
import math
from scipy.interpolate import CubicSpline
import sys

#============================================================================================================================================================
# reading the data from EMT data thief file
#============================================================================================================================================================
def compute_p_eff(x_new, sigma_x_new=None):
    base = "/home/abhinav/david/"
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
    
    
    #============================================================================================================================================================
    # plotting all the shear data from emt
    #============================================================================================================================================================
    
    fig = plt.figure(figsize=(8, 4), dpi=300)
    plt.scatter(data_e_6[:,0], data_e_6[:,1], c = 'b', marker = 'o', alpha = 0.3, label = '$\kappa = 10^{-6}$')
    plt.scatter(data_e_5[:,0], data_e_5[:,1], c = 'r', marker = 's', alpha = 0.3, label = '$\kappa = 10^{-5}$')
    plt.scatter(data_e_4[:,0], data_e_4[:,1], c = 'g', marker = 'd', alpha = 0.3, label = '$\kappa = 10^{-4}$')
    plt.xlabel('$ p $', fontsize = 15)
    plt.ylabel('$ G $', fontsize = 15)
    plt.xticks(fontsize = 15)
    plt.yticks(fontsize = 15)
    plt.yscale('log')
    plt.subplots_adjust(top=0.92, left = 0.17, bottom = 0.15, right = 0.98)
    plt.legend()
    plt.savefig(base+emt_folder+"shear_modulus_log.png")  # just plotting the last step data
    plt.close()
    plt.clf()
    
    alpha_m_list_e_6_emt = 4*shear_e_6_emt/np.sqrt(3)
    alpha_m_list_e_5_emt = 4*shear_e_5_emt/np.sqrt(3)
    alpha_m_list_e_4_emt = 4*shear_e_4_emt/np.sqrt(3)
    
    bulk_m_list_e_6_emt = 0.5*alpha_m_list_e_6_emt*np.sqrt(3)
    bulk_m_list_e_5_emt = 0.5*alpha_m_list_e_5_emt*np.sqrt(3)
    bulk_m_list_e_4_emt = 0.5*alpha_m_list_e_4_emt*np.sqrt(3)
    
    # plotting all the alpha_m data from emt
    fig = plt.figure(figsize=(8, 4), dpi=300)
    plt.scatter(data_e_6[:,0], alpha_m_list_e_6_emt, c = 'b', marker = 'o', alpha = 0.3, label = '$\kappa = 10^{-6}$')
    plt.scatter(data_e_5[:,0], alpha_m_list_e_5_emt, c = 'r', marker = 's', alpha = 0.3, label = '$\kappa = 10^{-5}$')
    plt.scatter(data_e_4[:,0], alpha_m_list_e_4_emt, c = 'g', marker = 'd', alpha = 0.3, label = '$\kappa = 10^{-4}$')
    plt.xlabel('$ p $', fontsize = 15)
    plt.ylabel('$ \\alpha_{m} $', fontsize = 15)
    plt.xticks(fontsize = 15)
    plt.yticks(fontsize = 15)
    plt.yscale('log')
    plt.subplots_adjust(top=0.92, left = 0.17, bottom = 0.15, right = 0.98)
    plt.legend()
    plt.savefig(base+emt_folder+"alpha_m_log.png")  # just plotting the last step data
    plt.close()
    plt.clf()
    
    # plotting all the bulk moduli data from emt
    fig = plt.figure(figsize=(8, 4), dpi=300)
    plt.scatter(data_e_6[:,0], bulk_m_list_e_6_emt, c = 'b', marker = 'o', alpha = 0.3, label = '$\kappa = 10^{-6}$')
    plt.scatter(data_e_5[:,0], bulk_m_list_e_5_emt, c = 'r', marker = 's', alpha = 0.3, label = '$\kappa = 10^{-5}$')
    plt.scatter(data_e_4[:,0], bulk_m_list_e_4_emt, c = 'g', marker = 'd', alpha = 0.3, label = '$\kappa = 10^{-4}$')
    plt.xlabel('$ p $', fontsize = 15)
    plt.ylabel(' Bulk Modulus ', fontsize = 15)
    plt.xticks(fontsize = 15)
    plt.yticks(fontsize = 15)
    plt.yscale('log')
    plt.subplots_adjust(top=0.92, left = 0.17, bottom = 0.15, right = 0.98)
    plt.legend()
    plt.savefig(base+emt_folder+"bulk_modulus_log.png")  # just plotting the last step data
    plt.close()
    plt.clf()
    
    #============================================================================================================================================================
    # finding p_effective using cubic splines on alpha_m vs p data
    # this alpha_m is gotten by dfar calculations
    #============================================================================================================================================================
    # Sample data points
    x = alpha_m_list_e_6_emt    # x is the alpha_m values
    y = data_e_6[:,0]   # using only kappa e-6 for now; y is p value of the network
    
    # Create the cubic spline interpolator
    f = CubicSpline(x, y, bc_type='natural')  # 'natural' for natural boundary conditions
    
    #============================================================================================================================================================
    # using the function to get the p_eff values
    #============================================================================================================================================================
    
    # Interpolate y values using the spline
    y_new = f(x_new)
    if sigma_x_new is not None:
        dy_dx = f(x_new, 1)                    # the derivative of the function
        sigma_y = np.abs(dy_dx) * sigma_x_new  # the error in y
    else:
        sigma_y = None
    
    return y_new, sigma_y
