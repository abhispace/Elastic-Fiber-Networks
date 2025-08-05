#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 31 11:55:13 2025

@author: abhinav
"""
import subprocess
 
from math import log10, floor
def find_exp_base(number):
    exp = floor(log10(abs(number)))
    return round(number/10**exp, 2), exp 


# Arguments to be passed to the called script
# srand_list = ['112']
# for 10 seoarate networks with different dipole positions
# num_center_temp = ['1','2','5','10']
num_center_temp = ['1','2','5','10','15','20','25','30','35','40']
# srand_list = ['112','113', '114','116','117','118','119','120','121','122']                         # 115 is BAD!!!!!
srand_list = ['113', '114','116','117','118','119','120','121','122']                         # 115 is BAD!!!!!
# srand_list = ['112']


# for the L = 128 case
# num_center_temp = ['1','8','20','40','60']
# srand_list = ['116','117','118','119','120']                         # 115 is BAD!!!!!

# num_center_temp = ['80','100','120','140','160']
# # num_center_temp = ['100','120','140','160']
# srand_list = ['116','117']                         # 115 is BAD!!!!!

# Run the called script with arguments
for srand_val in srand_list:
        for num_center in num_center_temp:
            subprocess.run(['python', 'dipole_moment_plots.py', num_center, srand_val])
