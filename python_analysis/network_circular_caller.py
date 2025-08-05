#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 17 17:01:03 2024
it calls network_circular.py that is primarily used to plot the network
but in this instance, we use it to output the displacements of each node to a file.

This code also fits the displacements!!!!!
It also produces plots for EMT data such as alpha_m and p_eff, shear modulus and bulk modulus!!


NOTE: srand 114, nd = 30 has inf. need to come back and solve this later.

It also makes p_eff vs N_d plots.

@author: abhinav
"""
import subprocess
 
from math import log10, floor
def find_exp_base(number):
    exp = floor(log10(abs(number)))
    return round(number/10**exp, 2), exp 


# Arguments to be passed to the called script
# srand_list = ['112']
# num_center_temp = ['1']
diff_network_list = ['111111,101111', '111111,111111', '111111,121111', '111111,131111', '111111,141111', '111111,151111', '111111,161111', '111111,171111', '111111,181111', '111111,191111',
                     '222222,202222', '222222,212222', '222222,222222', '222222,232222', '222222,242222', '222222,252222', '222222,262222', '222222,272222', '222222,282222', '222222,292222',
                     '333333,303333', '333333,313333', '333333,323333', '333333,333333', '333333,343333', '333333,353333', '333333,363333', '333333,373333', '333333,383333', '333333,393333',
                     '444444,404444', '444444,414444', '444444,424444', '444444,434444', '444444,444444', '444444,454444', '444444,464444', '444444,474444', '444444,484444', '444444,494444',
                     '555555,505555', '555555,515555', '555555,525555', '555555,535555', '555555,545555', '555555,555555', '555555,565555', '555555,575555', '555555,585555', '555555,595555',
                     '666666,606666', '666666,616666', '666666,626666', '666666,636666', '666666,646666', '666666,656666', '666666,666666', '666666,676666', '666666,686666', '666666,696666',
                     '777777,707777', '777777,717777', '777777,727777', '777777,737777', '777777,747777', '777777,757777', '777777,767777', '777777,777777', '777777,787777', '777777,797777',
                     '888888,808888', '888888,818888', '888888,828888', '888888,838888', '888888,848888', '888888,858888', '888888,868888', '888888,878888', '888888,888888', '888888,898888',
                     '999999,909999', '999999,919999', '999999,929999', '999999,939999', '999999,949999', '999999,959999', '999999,969999', '999999,979999', '999999,989999', '999999,999999',
                     '900000,900000', '900000,900001', '900000,900002', '900000,900003', '900000,900004', '900000,900005', '900000,900006', '900000,900007', '900000,900008', '900000,900009'] 

# for 10 seoarate networks with different dipole positions
diff_network_list = ['111111,101111', '111111,111111', '111111,121111', '111111,131111', '111111,141111', '111111,151111', '111111,161111', '111111,171111', '111111,181111', '111111,191111']
# num_center_temp = ['1','2','5','10']
num_center_temp = ['1','2','5','10','15','20','25','30','35','40']
# num_center_temp = ['1','15','20','25','30','35','40']

# for L = 128: N_d = 1, 8, 20, 40, 60
# diff_network_list = ['111111,101111', '111111,111111']
# num_center_temp = ['1', '8', '20', '40', '60']
# num_center_temp = ['4']
# srand_list = ['116','117','118','119','120']                         # 115 is BAD!!!!!

# num_center_temp = ['80', '100', '120', '140', '160']
# srand_list = ['116','117']                         # 115 is BAD!!!!!

# num_center_temp = ['1']
# srand_list = ['113', '114','116','117','118','119','120','121','122']                         # 115 is BAD!!!!!
srand_list = ['112']
# diff_network_list = ['111111,191111']

# Run the called script with arguments
calling_flag = 1
p1_calling_flag = 0

# num_center = 5
if calling_flag == 1:
    for srand_val in srand_list:
        for diff_network in diff_network_list:
            for num_center in num_center_temp:
                subprocess.run(['python', 'network_circular.py', diff_network, num_center, srand_val])
        # subprocess.run(['python', 'network_circular.py', diff_network, str(num_center)])
        
# for p=1::
if p1_calling_flag == 1:
    for srand_val in srand_list:
        for num_center in num_center_temp:
            subprocess.run(['python', 'network_circular.py', num_center, srand_val])

''' 
#==============================================================================
# analysis below
import numpy as np
import matplotlib.pyplot as plt
import math

num_center = 40
num_dip = num_center*6
num = 11
pbond = 0.55
# pbond = 0.6

mu = 1
mu_c = 1
tol = 1.0e-7

kappa = 1e-6
kappa_str = "%.2e" % kappa # to write the filename

rlen = 0.9
rlen_txt = "%.4f" % rlen
rlen_txt_ring = "%.9f" % 0.1

srand = 112

full_flag = 1  # this is set to 1 to consider also the N_d = 1,2,5 values in the interpolation of p_eff values; otherwise set to zero

# automating this code lol. called by ...caller_caller lol
auto_flag = 0
if auto_flag == 1:
    import sys
    num_center = int(sys.argv[1])    
    num_dip = num_center*6
    srand = int(sys.argv[2])


pbond_string = "%.2f" % pbond # to write the filename
tol_str = "%.2e" % tol # to write the filename
    
mu_str = "%.4f" % mu # to write the filename
mu_c_str = "%.4f" % mu_c # to write the filename

base = "/home/abhinav/david/"

cluster_new_flag = 1 
hex_flag = 0
hex_rand_flag = 0
cluster_radial_flag = 1
    
ranseed = '667720,601210'    # -- network 1

#=========================================================================================================================================    
#=========================================================================================================================================    

if num_center == 2 and srand == 112:
    diff_network_list = ['111111,101111', '111111,111111', '111111,121111', '111111,131111', '111111,141111', '111111,151111', '111111,161111', '111111,171111', '111111,181111', '111111,191111',
                         '222222,202222', '222222,212222', '222222,222222', '222222,232222', '222222,242222', '222222,252222', '222222,262222', '222222,272222', '222222,282222', '222222,292222',
                         '333333,303333', '333333,313333', '333333,323333', '333333,333333', '333333,343333', '333333,353333', '333333,363333', '333333,373333', '333333,383333', '333333,393333',
                         '444444,404444', '444444,414444', '444444,424444', '444444,434444', '444444,444444', '444444,464444', '444444,474444', '444444,484444', '444444,494444',
                         '555555,505555', '555555,515555', '555555,525555', '555555,535555', '555555,545555', '555555,555555', '555555,565555', '555555,575555', '555555,585555', '555555,595555',
                         '666666,606666', '666666,616666', '666666,626666', '666666,636666', '666666,646666', '666666,656666', '666666,666666', '666666,676666', '666666,686666', '666666,696666',
                         '777777,707777', '777777,717777', '777777,727777', '777777,737777', '777777,747777', '777777,757777', '777777,767777', '777777,777777', '777777,787777', '777777,797777',
                         '888888,808888', '888888,818888', '888888,828888', '888888,838888', '888888,848888', '888888,858888', '888888,868888', '888888,878888', '888888,888888', '888888,898888',
                         '999999,909999', '999999,919999', '999999,929999', '999999,939999', '999999,949999', '999999,959999', '999999,969999', '999999,979999', '999999,989999', '999999,999999',
                         '900000,900000', '900000,900001', '900000,900002', '900000,900003', '900000,900004', '900000,900005', '900000,900006', '900000,900007', '900000,900008', '900000,900009']

if num_center == 5 and srand == 112:
    diff_network_list = ['111111,101111', '111111,111111', '111111,121111', '111111,131111', '111111,141111', '111111,161111', '111111,171111', '111111,181111', '111111,191111',
                          '222222,202222', '222222,212222', '222222,222222', '222222,232222', '222222,242222', '222222,252222', '222222,262222', '222222,272222', '222222,282222', '222222,292222',
                          '333333,303333', '333333,313333', '333333,323333', '333333,333333', '333333,343333', '333333,353333', '333333,363333', '333333,373333', '333333,383333', '333333,393333',
                          '444444,404444', '444444,414444', '444444,424444', '444444,434444', '444444,444444', '444444,464444', '444444,474444', '444444,484444', '444444,494444',
                          '555555,505555', '555555,515555', '555555,525555', '555555,535555', '555555,545555', '555555,555555', '555555,565555', '555555,575555', '555555,585555', '555555,595555',
                          '666666,606666', '666666,616666', '666666,626666', '666666,636666', '666666,646666', '666666,656666', '666666,666666', '666666,676666', '666666,686666', '666666,696666',
                          '777777,707777', '777777,717777', '777777,727777', '777777,737777', '777777,747777', '777777,757777', '777777,767777', '777777,777777', '777777,787777', '777777,797777',
                          '888888,818888', '888888,828888', '888888,838888', '888888,848888', '888888,858888', '888888,868888', '888888,878888', '888888,888888', '888888,898888',
                          '999999,909999', '999999,919999', '999999,929999', '999999,939999', '999999,949999', '999999,959999', '999999,969999', '999999,979999', '999999,989999', '999999,999999',
                          '900000,900000', '900000,900001', '900000,900002', '900000,900003', '900000,900005', '900000,900006', '900000,900007', '900000,900008', '900000,900009']

if (num_center == 10 or num_center == 15) and srand == 112:
    diff_network_list = ['111111,101111', '111111,111111', '111111,121111', '111111,131111', '111111,141111', '111111,161111', '111111,171111', '111111,181111', '111111,191111',
                          '222222,202222', '222222,212222', '222222,222222', '222222,232222', '222222,242222', '222222,252222', '222222,262222', '222222,272222', '222222,282222', '222222,292222',
                          '333333,303333', '333333,313333', '333333,323333', '333333,333333', '333333,343333', '333333,353333', '333333,363333', '333333,373333', '333333,383333', '333333,393333',
                          '444444,404444', '444444,414444', '444444,424444', '444444,434444', '444444,444444', '444444,464444', '444444,474444', '444444,484444', '444444,494444',
                          '555555,505555', '555555,515555', '555555,525555', '555555,535555', '555555,545555', '555555,555555', '555555,565555', '555555,575555', '555555,585555', '555555,595555',
                          '666666,606666', '666666,616666', '666666,626666', '666666,636666', '666666,646666', '666666,656666', '666666,666666', '666666,676666', '666666,686666', '666666,696666',
                          '777777,707777', '777777,717777', '777777,727777', '777777,737777', '777777,747777', '777777,757777', '777777,767777', '777777,777777', '777777,787777', '777777,797777',
                          '888888,818888', '888888,828888', '888888,838888', '888888,848888', '888888,858888', '888888,868888', '888888,888888', '888888,898888',
                          '999999,909999', '999999,919999', '999999,929999', '999999,939999', '999999,949999', '999999,959999', '999999,969999', '999999,979999', '999999,989999', '999999,999999',
                          '900000,900000', '900000,900001', '900000,900002', '900000,900003', '900000,900005', '900000,900006', '900000,900007', '900000,900008', '900000,900009']
    

if num_center == 20 and srand == 112:
    diff_network_list = ['111111,101111', '111111,111111', '111111,121111', '111111,131111', '111111,141111', '111111,161111', '111111,171111', '111111,181111', '111111,191111',
                          '222222,202222', '222222,212222', '222222,222222', '222222,232222', '222222,242222', '222222,252222', '222222,262222', '222222,272222', '222222,282222', '222222,292222',
                          '333333,303333', '333333,313333', '333333,323333', '333333,333333', '333333,343333', '333333,353333', '333333,363333', '333333,373333', '333333,383333', '333333,393333',
                          '444444,404444', '444444,414444', '444444,424444', '444444,434444', '444444,464444', '444444,484444', '444444,494444',
                          '555555,505555', '555555,515555', '555555,525555', '555555,545555', '555555,555555', '555555,565555', '555555,575555', '555555,585555', '555555,595555',
                          '666666,606666', '666666,616666', '666666,626666', '666666,636666', '666666,646666', '666666,656666', '666666,666666', '666666,676666', '666666,686666', '666666,696666',
                          '777777,707777', '777777,717777', '777777,727777', '777777,737777', '777777,747777', '777777,757777', '777777,767777', '777777,777777', '777777,787777', '777777,797777',
                          '888888,818888', '888888,828888', '888888,838888', '888888,848888', '888888,858888', '888888,868888', '888888,878888', '888888,888888', '888888,898888',
                          '999999,909999', '999999,919999', '999999,929999', '999999,939999', '999999,949999', '999999,959999', '999999,969999', '999999,979999', '999999,989999', '999999,999999',
                          '900000,900000', '900000,900001', '900000,900002', '900000,900003', '900000,900005', '900000,900006', '900000,900007', '900000,900008', '900000,900009']

if num_center == 25 and srand == 112:
    diff_network_list = ['111111,101111', '111111,111111', '111111,121111', '111111,131111', '111111,141111', '111111,161111', '111111,171111', '111111,181111', '111111,191111',
                          '222222,202222', '222222,212222', '222222,222222', '222222,232222', '222222,242222', '222222,252222', '222222,262222', '222222,272222', '222222,282222', '222222,292222',
                          '333333,303333', '333333,313333', '333333,323333', '333333,333333', '333333,343333', '333333,353333', '333333,363333', '333333,373333', '333333,383333', '333333,393333',
                          '444444,404444', '444444,414444', '444444,424444', '444444,464444', '444444,484444', '444444,494444',
                          '555555,505555', '555555,515555', '555555,525555', '555555,545555', '555555,565555', '555555,575555', '555555,585555', '555555,595555',
                          '666666,606666', '666666,616666', '666666,626666', '666666,636666', '666666,646666', '666666,656666', '666666,666666', '666666,676666', '666666,686666', '666666,696666',
                          '777777,707777', '777777,717777', '777777,727777', '777777,737777', '777777,747777', '777777,757777', '777777,767777', '777777,777777', '777777,787777', '777777,797777',
                          '888888,818888', '888888,828888', '888888,838888', '888888,848888', '888888,858888', '888888,868888', '888888,878888', '888888,888888', '888888,898888',
                          '999999,909999', '999999,919999', '999999,929999', '999999,939999', '999999,949999', '999999,959999', '999999,979999', '999999,989999', '999999,999999',
                          '900000,900000', '900000,900001', '900000,900002', '900000,900003', '900000,900005', '900000,900006', '900000,900007', '900000,900008', '900000,900009']

if num_center == 30 and srand == 112:
    diff_network_list = ['111111,101111', '111111,121111', '111111,131111', '111111,141111', '111111,161111', '111111,171111', '111111,181111', '111111,191111',
                         '222222,202222', '222222,212222', '222222,222222', '222222,232222', '222222,242222', '222222,252222', '222222,262222', '222222,282222', '222222,292222',
                         '333333,303333', '333333,313333', '333333,323333', '333333,333333', '333333,363333', '333333,393333',
                         '444444,404444', '444444,414444', '444444,424444', '444444,444444', '444444,484444', '444444,494444',
                         '555555,505555', '555555,515555', '555555,545555', '555555,565555', '555555,575555', '555555,585555', '555555,595555',
                         '666666,606666', '666666,616666', '666666,626666', '666666,636666', '666666,646666', '666666,656666', '666666,666666', '666666,676666', '666666,686666',
                         '777777,707777', '777777,717777', '777777,727777', '777777,737777', '777777,747777', '777777,757777', '777777,767777', '777777,777777', '777777,787777', '777777,797777',
                         '888888,818888', '888888,828888', '888888,838888', '888888,848888', '888888,858888', '888888,868888', '888888,878888', '888888,888888', '888888,898888',
                         '999999,909999', '999999,929999', '999999,929999', '999999,949999', '999999,959999', '999999,979999', '999999,989999', '999999,999999',
                         '900000,900000', '900000,900002', '900000,900003', '900000,900005', '900000,900006', '900000,900007', '900000,900008', '900000,900009']

if num_center == 35 and srand == 112:
    diff_network_list = ['111111,101111', '111111,121111', '111111,131111', '111111,141111', '111111,161111', '111111,171111', '111111,181111', '111111,191111',
                         '222222,212222', '222222,222222', '222222,232222', '222222,242222', '222222,252222', '222222,262222', '222222,282222', '222222,292222',
                         '333333,303333', '333333,313333', '333333,323333', '333333,333333', '333333,363333', '333333,393333',
                         '444444,404444', '444444,414444', '444444,424444', '444444,444444', '444444,484444', '444444,494444',
                         '555555,505555', '555555,515555', '555555,545555', '555555,565555', '555555,575555', '555555,585555', '555555,595555',
                         '666666,606666', '666666,616666', '666666,636666', '666666,646666', '666666,656666', '666666,666666', '666666,676666', '666666,686666',
                         '777777,707777', '777777,717777', '777777,727777', '777777,737777', '777777,757777', '777777,767777', '777777,777777', '777777,787777', '777777,797777',
                         '888888,818888', '888888,828888', '888888,838888', '888888,848888', '888888,858888', '888888,868888', '888888,878888', '888888,888888', '888888,898888',
                         '999999,909999', '999999,929999', '999999,929999', '999999,949999', '999999,959999', '999999,979999', '999999,989999', '999999,999999',
                         '900000,900000', '900000,900002', '900000,900005', '900000,900006', '900000,900007', '900000,900008', '900000,900009']


if num_center == 40 and srand == 112:
    diff_network_list = ['111111,101111', '111111,121111', '111111,131111', '111111,141111', '111111,161111', '111111,171111', '111111,181111', '111111,191111',
                         '222222,212222', '222222,222222', '222222,232222', '222222,242222', '222222,252222', '222222,282222', '222222,292222',
                         '333333,313333', '333333,323333', '333333,333333', '333333,363333', '333333,393333',
                         '444444,404444', '444444,414444', '444444,424444', '444444,444444', '444444,484444',
                         '555555,505555', '555555,515555', '555555,545555', '555555,565555', '555555,575555', '555555,585555', '555555,595555',
                         '666666,606666', '666666,616666', '666666,636666', '666666,646666', '666666,656666', '666666,666666', '666666,686666',
                         '777777,707777', '777777,717777', '777777,727777', '777777,737777', '777777,757777', '777777,767777', '777777,777777', '777777,797777',
                         '888888,828888', '888888,838888', '888888,848888', '888888,858888', '888888,868888', '888888,878888', '888888,888888', '888888,898888',
                         '999999,909999', '999999,929999', '999999,929999', '999999,949999', '999999,959999', '999999,979999', '999999,989999', '999999,999999',
                         '900000,900000', '900000,900002', '900000,900005', '900000,900006', '900000,900007', '900000,900008', '900000,900009']


if num_center == 40 and kappa == 1e-5 and srand == 112:
    diff_network_list = ['111111,101111', '111111,121111', '111111,131111', '111111,141111', '111111,171111', '111111,181111', '111111,191111',
                         '222222,212222', '222222,222222', '222222,232222', '222222,242222', '222222,252222', '222222,282222', '222222,292222',
                         '333333,313333', '333333,323333', '333333,333333', '333333,363333', '333333,393333',
                         '444444,404444', '444444,414444', '444444,424444', '444444,444444', '444444,484444',
                         '555555,505555', '555555,515555', '555555,545555', '555555,565555', '555555,575555', '555555,585555', '555555,595555',
                         '666666,606666', '666666,616666', '666666,636666', '666666,646666', '666666,656666', '666666,666666', '666666,686666',
                         '777777,707777', '777777,717777', '777777,727777', '777777,737777', '777777,757777', '777777,767777', '777777,777777', '777777,797777',
                         '888888,828888', '888888,838888', '888888,848888', '888888,858888', '888888,868888', '888888,878888', '888888,888888', '888888,898888',
                         '999999,909999', '999999,929999', '999999,929999', '999999,949999', '999999,959999', '999999,979999', '999999,989999', '999999,999999',
                         '900000,900000', '900000,900002', '900000,900005', '900000,900006', '900000,900007', '900000,900008', '900000,900009']


if srand != 112:
    diff_network_list = ['111111,101111', '111111,111111', '111111,121111', '111111,131111', '111111,141111', '111111,151111', '111111,161111', '111111,171111', '111111,181111', '111111,191111']
    # if (srand == 114 and num_center == 10) or (srand == 120 and num_center == 10):
    #     diff_network_list = ['111111,101111', '111111,111111', '111111,121111', '111111,131111', '111111,141111', '111111,151111', '111111,161111', '111111,171111', '111111,191111']
    if (srand == 114 and num_center == 30):
        diff_network_list = ['111111,101111', '111111,111111', '111111,121111', '111111,131111', '111111,141111', '111111,151111', '111111,161111', '111111,171111', '111111,191111']
    elif (srand == 113 and num_center == 30):
        diff_network_list = ['111111,101111', '111111,111111', '111111,121111', '111111,131111', '111111,151111', '111111,161111', '111111,171111', '111111,181111', '111111,191111']
    elif (srand == 117 and num_center == 40):
        diff_network_list = ['111111,101111', '111111,111111', '111111,121111', '111111,131111', '111111,141111', '111111,151111', '111111,171111', '111111,181111', '111111,191111']
    # elif (srand == 122 and num_center == 10):
    #     diff_network_list = ['111111,101111', '111111,121111', '111111,131111', '111111,141111', '111111,151111', '111111,161111', '111111,171111', '111111,181111', '111111,191111']

#***********************************************************************************************************************************************
# writing the bending dominated cases to file
#***********************************************************************************************************************************************
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

fout_folder = "lattice_nopbd_bndry_all_clamp_restlength_circle_inner_fixed_radial_arp_cluster/"
if cluster_new_flag == 1 and hex_flag == 0 and hex_rand_flag == 0 and cluster_radial_flag == 0:
    fout_folder = "cluster_download/lattice_nopbd_bndry_all_clamp_restlength_circle_inner_fixed_radial_arp_bash/"
elif cluster_new_flag == 1 and hex_flag == 1 and hex_rand_flag == 0 and cluster_radial_flag == 0:
    fout_folder = "cluster_download/lattice_nopbd_bndry_all_clamp_restlength_circle_inner_fixed_radial_arp_bash_hex/"
elif cluster_new_flag == 1 and hex_flag == 0 and hex_rand_flag == 1 and cluster_radial_flag == 0:
    fout_folder = "cluster_download/lattice_nopbd_bndry_all_clamp_restlength_circle_inner_fixed_radial_arp_bash_random_dipole_hex_bonds/"
elif cluster_new_flag == 1 and hex_flag == 0 and hex_rand_flag == 0 and cluster_radial_flag == 1:
    fout_folder = "cluster_download/lattice_nopbd_bndry_all_clamp_restlength_circle_inner_fixed_radial_arp_bash_radial/"

# bend_outfname = base+fout_folder + kappa_fname +ranseed+"/"+str(srand)+"/means/"+"bend_dom_"+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+".txt"
# bend_dom_cases = np.loadtxt(bend_outfname, dtype=str)
# print('Bending dominated cases for set p, N and kappa and srand read from file: ', bend_outfname)

#=========================================================================================================================================    
#=========================================================================================================================================    
    
    

# reading mean displacements
data_all = []
for diff_network in diff_network_list:
# for diff_network in bend_dom_cases:
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
    
    if cluster_new_flag == 1 and hex_flag == 1 and hex_rand_flag == 0 and cluster_radial_flag == 0:
        folder = 'cluster_download/lattice_nopbd_bndry_all_clamp_restlength_circle_inner_fixed_radial_arp_bash_hex/'+ kappa_fname +ranseed+"/"+str(srand)+"/"+diff_network+"/txt/displacement/"
        plot_folder = base + 'cluster_download/lattice_nopbd_bndry_all_clamp_restlength_circle_inner_fixed_radial_arp_bash_hex/'+ kappa_fname +ranseed+"/"+str(srand)
        disp_folder = base + 'cluster_download/lattice_nopbd_bndry_all_clamp_restlength_circle_inner_fixed_radial_arp_bash_hex/'+ kappa_fname +ranseed+"/"+str(srand)+"/means/"
        ranseed_folder = base + 'cluster_download/lattice_nopbd_bndry_all_clamp_restlength_circle_inner_fixed_radial_arp_bash_hex/'+ kappa_fname +ranseed+"/"
    elif cluster_new_flag == 1 and hex_flag == 0 and hex_rand_flag == 1 and cluster_radial_flag == 0:
        folder = 'cluster_download/lattice_nopbd_bndry_all_clamp_restlength_circle_inner_fixed_radial_arp_bash_random_dipole_hex_bonds/'+ kappa_fname +ranseed+"/"+str(srand)+"/"+diff_network+"/txt/displacement/"
        plot_folder = base + 'cluster_download/lattice_nopbd_bndry_all_clamp_restlength_circle_inner_fixed_radial_arp_bash_random_dipole_hex_bonds/'+ kappa_fname +ranseed+"/"+str(srand)
        disp_folder = base + 'cluster_download/lattice_nopbd_bndry_all_clamp_restlength_circle_inner_fixed_radial_arp_bash_random_dipole_hex_bonds/'+ kappa_fname +ranseed+"/"+str(srand)+"/means/"
        ranseed_folder = base + 'cluster_download/lattice_nopbd_bndry_all_clamp_restlength_circle_inner_fixed_radial_arp_bash_random_dipole_hex_bonds/'+ kappa_fname +ranseed+"/"
    elif cluster_new_flag == 1 and hex_flag == 0 and hex_rand_flag == 0 and cluster_radial_flag == 1:
        folder = 'cluster_download/lattice_nopbd_bndry_all_clamp_restlength_circle_inner_fixed_radial_arp_bash_radial/'+ kappa_fname +ranseed+"/"+str(srand)+"/"+diff_network+"/txt/displacement/"
        plot_folder = base + 'cluster_download/lattice_nopbd_bndry_all_clamp_restlength_circle_inner_fixed_radial_arp_bash_radial/'+ kappa_fname +ranseed+"/"+str(srand)
        disp_folder = base + 'cluster_download/lattice_nopbd_bndry_all_clamp_restlength_circle_inner_fixed_radial_arp_bash_radial/'+ kappa_fname +ranseed+"/"+str(srand)+"/means/"
        ranseed_folder = base + 'cluster_download/lattice_nopbd_bndry_all_clamp_restlength_circle_inner_fixed_radial_arp_bash_radial/'+ kappa_fname +ranseed+"/"

#     fname = base+folder+"bndry_node_radial_disp_mean_"+str(srand)+"_"+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+".txt"    
#     data = np.loadtxt(fname)

#     data_all.append(data)
#     # i += 1
    
# data_all = np.array(data_all)
# bins = data_all[0,:,0]
# # finding the mean of all the networks
# disp = data_all[:,:,1]
# mean_disp = np.nanmean(disp,axis=0)

disp_mean_fname = disp_folder +"bndry_node_radial_disp_mean_"+str(srand)+"_"+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+".txt"
print('reading the mean of (after binning) radial displacements on all nodes from file : ',disp_mean_fname)    
data_all = np.loadtxt(disp_mean_fname)
bins = data_all[:,0]
mean_disp = data_all[:,1]
std_disp = data_all[:,2]

# scatter plot - outer region :: only works for kappa = 1e-6 for now
r1 = 14
r2 = 25
bulk_m = 0.5*np.sqrt(3)*mu
shear_m = 0.25*np.sqrt(3)*mu

# xaxis of plot
# xdata = bins[12:]
# ydata = mean_disp[12:]  # yaxis only for outer region
# std_disp_outer = std_disp[12:]

xdata = bins[12:]
ydata = mean_disp[12:]  # yaxis only for outer region
std_disp_outer = std_disp[12:]

# first reading the standard deviation of every network then averaging them, then finding standard error of mean
# std_disp = np.mean(data_all[:,:,3], axis = 0)/np.sqrt(len(data_all))
func_prefactor = 1. / ( ((r1**2) / (2*((np.sqrt(3)/2) * r1**2 + ((np.sqrt(3)/4) * r2**2))))  * ( (r2**2 - xdata**2)/xdata ) )
fit_std = np.mean(std_disp_outer*func_prefactor)



if pbond == 1:
    sigma1 = -5.8e-4
    sigma1_outer = -5.8e-4

elif pbond == 0.9:
    sigma1 = -6.2e-4
    sigma1_outer = -7.0e-4

    shear_m = 3.04e-1
    bulk_m = 6.08e-01

elif pbond == 0.8:
    sigma1 = -6.6e-4
    sigma1_outer = -9.2e-4

    shear_m = 1.87e-1
    bulk_m = 3.74e-01

elif pbond == 0.55 and num_center == 1:
    sigma1 = -3.0e-8
    sigma1_outer = -7e-8
    shear_m = 1.2e-5
    bulk_m = 2.4e-05
elif pbond == 0.55 and num_center == 5:
    sigma1 = -3.0e-8
    sigma1_outer = -20e-8
    shear_m = 1.2e-5
    bulk_m = 2.4e-05
elif pbond == 0.55 and num_center == 10:
    sigma1 = -3.0e-8
    sigma1_outer = -34e-8
    shear_m = 1.2e-5
    bulk_m = 2.4e-05
    
fit_x = np.linspace(bins[0],bins[-1],100)
fit_x_outer = np.linspace(r1+1, r2, 100)     # outer region x range
# u_r = ( (sigma1 * r1**2) / (2*(bulk_m * r1**2+ shear_m * r2**2)) ) * ( (r2**2 - fit_x**2)/fit_x )       #  theory ar + b/r fit for whole range 
# u_r_outer = ( (sigma1_outer * r1**2) / (2*(bulk_m * r1**2+ shear_m * r2**2)) ) * ( (r2**2 - fit_x_outer**2)/fit_x_outer )       #  theory ar + b/r fit for whole range 

# fitting data
# import scipy.optimize
from scipy.optimize import curve_fit
# def func(x):
#     return ( (x * r1**2) / (2*(bulk_m * r1**2+ shear_m * r2**2)) ) * ( (r2**2 - fit_x_outer**2)/fit_x_outer )
def func(x,a):
    return (a) * (r1**2) / (2*((np.sqrt(3)/2) * r1**2 + ((np.sqrt(3)/4) * r2**2)))  * ( (r2**2 - x**2)/x )

# f2str = '{:.2f}'
popt, pcov = curve_fit(func, xdata, ydata, maxfev=10000)
# label = '$%s x^{-%s}$' % tuple(f2str.format(t) for t in [popt[0], popt[1]])
print('Fit parameters: ', popt)


# scatter plot of binned data (wrt initial radial distance) and also fits
fig = plt.figure(figsize=(8, 4), dpi=300)
plt.scatter(bins, mean_disp)
# plt.plot(fit_x, u_r, c = 'r')
# plt.plot(fit_x, y_fit, c = 'g')
plt.xlabel('Distance from center', fontsize = 15)
plt.ylabel('Radial Displacement', fontsize = 15)
plt.xticks(fontsize = 15)
plt.yticks(fontsize = 15)
# plt.ylim([-0.004,0])   # for p = 1; N = 1
# plt.ylim([-0.02,0.01])   # for p = 0.9; N = 1
# plt.ylim([-0.025,0.015])   # for p = 0.8; N = 1

# plt.xlim([12,radius+dr])
plt.subplots_adjust(top=0.92, left = 0.17, bottom = 0.15, right = 0.98)

plt.savefig(plot_folder+"/displacement_scatter_mean_"+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+mu_str+"_"+mu_c_str+"_"+kappa_str+"_"+tol_str+"_"+rlen_txt+"_"+str(num-1)+".png")  # just plotting the last step data
plt.close()
plt.clf()


z = np.polyfit(bins[12:], mean_disp[12:], 2)
f = np.poly1d(z)
# fit_label = ' %5.3e $x^{2}$ + %5.3e x + %5.3e ' % (f[2], f[1], f[0])

f2str = '{:.2e}'
fit_label = '$ %s x^{2} + %s x + %s $' % tuple(f2str.format(t) for t in [f[2], f[1], f[0]])

# calculate new x's and y's
# x_new = np.linspace(x[0], x[-1], 50)
y_fit = f(fit_x_outer)

fit_theory_outer_xdata = np.linspace(13, 25)
fit_theory_outer_ydata = func(fit_theory_outer_xdata, *popt)
base_num, exp = find_exp_base(popt[0])
base_num_err, exp_err = find_exp_base(np.sqrt(np.diag(pcov))[0])
label_txt = '$\\Sigma_{1}/\\mu_{m} =' + str(base_num) + ' \\times 10^{'+ str(exp) + '}$' + ' $ \\pm '+ str(base_num_err) + '\\times 10^{'+ str(exp_err) + '}$'
# print(base_num, )

fig = plt.figure(figsize=(8, 4), dpi=300)
# plt.scatter(bins, mean_disp)
plt.errorbar(bins, mean_disp, yerr = std_disp, marker = 'o', ls = 'none')
# plt.plot(fit_x_outer, u_r_outer, c = 'r')
# plt.plot(fit_x_outer, y_fit, c = 'g', label = fit_label)
# plt.plot(xdata, func(xdata, *popt), 'b-', label = '$ \\Sigma_{1}/\\mu_{m} = %0.2e \\pm %0.2e $' % (popt[0], np.sqrt(np.diag(pcov))))    
# plt.plot(xdata, func(xdata, *popt), 'b-', label = label_txt)    
plt.plot(fit_theory_outer_xdata, fit_theory_outer_ydata, 'b-', label = label_txt)    
plt.xlabel('Distance from center', fontsize = 15)
plt.ylabel('Radial Displacement', fontsize = 15)
plt.xticks(fontsize = 15)
plt.yticks(fontsize = 15)
plt.ylim([np.min(ydata),0])   # for p = 1; N = 1

if pbond == 1:
    plt.ylim([-0.004,0.0])   # for p = 0.9; N = 1
if pbond == 0.9:
    plt.ylim([-0.005,0.0])   #   for p = 0.9; N = 1
if pbond == 0.8:
    plt.ylim([-0.007,0.0])   # for p = 0.8; N = 1
if pbond == 0.55 and num_center == 1:
    plt.ylim([-0.020,0.0])   # for p = 0.8; N = 1
if pbond == 0.55 and num_center == 5:
    plt.ylim([-0.05,0.01])   # for p = 0.8; N = 1
if pbond == 0.55 and num_center == 10:
    plt.ylim([-0.08,0.0])   # for p = 0.8; N = 1

plt.xlim([12+0.5,bins[-1]+0.5])
plt.subplots_adjust(top=0.92, left = 0.17, bottom = 0.15, right = 0.98)
plt.legend(loc = 'best')

plt.savefig(plot_folder+"/displacement_scatter_mean_outer_"+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+mu_str+"_"+mu_c_str+"_"+kappa_str+"_"+tol_str+"_"+rlen_txt+"_"+str(num-1)+".png")  # just plotting the last step data
plt.close()
plt.clf()

#==============================================================================
# commented out because it is already being done in dipole_moment_plots.py
# # writing binned data means!!
# disp_mean_outfname = disp_folder+"bndry_node_radial_disp_mean_"+str(srand)+"_"+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+".txt"
# print('here all mean (after binning) radial displacements on all nodes are written to file : ',disp_mean_outfname)    
# heading = 'bin mid       Mean Radial dist        Meadian Radial dist        Std Radial dist'
# fmt = '%10.2f', '%15.7e', '%15.7e'
# np.savetxt(disp_mean_outfname, np.column_stack((bins, mean_disp, std_disp)), header = heading, fmt = fmt)

#==============================================================================
# alpha_m vs p plot below
#==============================================================================
p_list = np.array([0.55, 0.6, 0.7, 0.8, 0.9, 1])
shear_m_list = np.array([1.2e-5, 3.1e-5, 0.05, 0.187, 0.304, 0.433])
alpha_m_list = 4*shear_m_list/np.sqrt(3)
bulk_m_list = 0.5*alpha_m_list*np.sqrt(3)

fig = plt.figure(figsize=(8, 4), dpi=300)
plt.scatter(p_list, alpha_m_list)
plt.xlabel('$ p $', fontsize = 15)
plt.ylabel('$ \\alpha_{m} $', fontsize = 15)
plt.xticks(fontsize = 15)
plt.yticks(fontsize = 15)

# plt.xlim([12+0.5,bins[-1]+0.5])
plt.subplots_adjust(top=0.92, left = 0.17, bottom = 0.15, right = 0.98)
plt.legend()

plt.savefig(base+"cluster_download/alpha_m_vs_p_passive_network.png")  # just plotting the last step data
plt.close()
plt.clf()


fig = plt.figure(figsize=(8, 4), dpi=300)
# plt.scatter(p_list, np.log(shear_m_list))
plt.scatter(p_list, alpha_m_list)
plt.xlabel('$ p $', fontsize = 15)
plt.ylabel('$ \\alpha_{m} $', fontsize = 15)
plt.xticks(fontsize = 15)
plt.yticks(fontsize = 15)
plt.yscale('log')
plt.subplots_adjust(top=0.92, left = 0.17, bottom = 0.15, right = 0.98)
# plt.legend()
plt.savefig(base+"cluster_download/alpha_m_vs_p_passive_network_log.png")  # just plotting the last step data
plt.close()
plt.clf()



#==============================================================================
# a* and b* for kappa = 1e-6
#==============================================================================
alpha = 1
a_star = (alpha*p_list - alpha_m_list)/(alpha - alpha_m_list)
b_star = 2/3 - a_star


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

# plotting all the shear data from emt
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
fig = plt.figure(figsize=(8, 5), dpi=300)
plt.scatter(data_e_6[:,0], alpha_m_list_e_6_emt, c = 'b', marker = 'o', s = 100, alpha = 0.3, label = '$\widetilde\kappa = 10^{-6}$')
plt.scatter(data_e_5[:,0], alpha_m_list_e_5_emt, c = 'r', marker = 's', s = 100, alpha = 0.3, label = '$\widetilde\kappa = 10^{-5}$')
plt.scatter(data_e_4[:,0], alpha_m_list_e_4_emt, c = 'g', marker = 'd', s = 100, alpha = 0.3, label = '$\widetilde\kappa = 10^{-4}$')
plt.xlabel('$ p $', fontsize = 20)
plt.ylabel('$ \\mu_{m} $', fontsize = 20)
plt.xticks(fontsize = 20)
plt.yticks(fontsize = 20)
plt.yscale('log')
ax = plt.gca()
# Increase tick width and length for major and minor ticks
ax.tick_params(axis='both', which='major', width=3, length=10)  # Major ticks
ax.tick_params(axis='both', which='minor', width=2, length=5)   # Minor ticks
plt.subplots_adjust(top=0.96, left = 0.15, bottom = 0.16, right = 0.98)
plt.legend(loc = 'best', fontsize = 14)
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

#==============================================================================
# finding p_effective using cubic splines on alpha_m vs p data
# this alpha_m is gotten by dfar calculations
#==============================================================================

from scipy.interpolate import CubicSpline

# Sample data points
x = alpha_m_list_e_6_emt    # x is the alpha_m values
y = data_e_6[:,0]   # using only kappa e-6 for now; y is p value of the network

# Create the cubic spline interpolator
f = CubicSpline(x, y, bc_type='natural')  # 'natural' for natural boundary conditions

# the following are for model 1 btw
# Generate new x values for interpolation
# the first value has changed so we are automatically reading it from the file
# if srand == 112:
#     x_new = np.array([4.30712478e-05, 2.49063480e-05, 3.45629444e-05, 2.97510556e-05,
#            3.41823708e-05, 5.47801744e-05, 6.97469465e-05, 8.90492913e-05,
#            2.05679724e-04, 3.14077173e-04])         # for kappa = 1e-6 case; these are alpha_m values for N = 1,2,5,10....40; also for Model 1
    
#     sigma_x_new = [1.56358896e-06, 1.92828857e-06, 1.38911084e-06, 1.13312057e-06,
#            1.19314532e-06, 2.31284766e-06, 3.52200806e-06, 4.90733186e-06,
#            4.98969141e-05, 6.54395199e-05]         # for kappa = 1e-6 case; these are errors in alpha_m values for N = 1,2,5,10....40; also for Model 1


# following data is taken from code: dipole_moment_means_srand.py.  At the end of that code, there are two variables:
# alpha_m_kappae_6: this is the x_new for this code 
# yerr_alpha_m_kappae_6: this is the sigma_x_new for this code

# following is the data when N = 1 is chosen as per the srand value at a random location
# reading the files
# the values corresponding to 112 has changed so we are automatically reading it from the file
# if srand != 112:
#     num_center_val = 10    # this was set before, when Nd went uto 10 for srand = 113...122 cases
#     num_center_val = 40    # now we can go upto 40 although srand 114 and 122 do not work anymore
#     pbond_string_val = '1.00'
#     fname_alpha = disp_folder +  "alpha_m_values_"+str(srand)+"_"+str(num_center_val)+"_"+str(num_center_val*6)+"_"+pbond_string_val+"_"+tol_str+"_"+kappa_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+".txt"
#     alpha_all_data = np.loadtxt(fname_alpha)
#     x_new = alpha_all_data[:,1]
#     sigma_x_new = alpha_all_data[:,2]

num_center_val = 40    # now we can go upto 40 although srand 114 and 122 do not work anymore
pbond_string_val = '1.00'
fname_alpha = disp_folder +  "alpha_m_values_"+str(srand)+"_"+str(num_center_val)+"_"+str(num_center_val*6)+"_"+pbond_string_val+"_"+tol_str+"_"+kappa_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+".txt"
alpha_all_data = np.loadtxt(fname_alpha)
x_new = alpha_all_data[:,1]
sigma_x_new = alpha_all_data[:,2]

# elif srand == 113:
#     x_new = [3.96048118e-05, 2.71653913e-05, 3.12904084e-05, 3.85422404e-05]
#     sigma_x_new = [3.45073106e-06, 3.43455982e-06, 3.42872923e-06, 3.82909106e-06]
# # # for srand 114:
# elif srand == 114:
#     x_new = [3.96048118e-05, 2.75433704e-05, 3.06398932e-05, 5.01094239e-05]
#     sigma_x_new = [3.45073106e-06, 2.61614237e-06, 3.09724538e-06, 5.02831748e-06]
# # # for srand 116:
# elif srand == 116:
#     x_new = [3.96048118e-05, 3.07689897e-05, 4.05965908e-05, 4.15942661e-05]
#     sigma_x_new = [3.45073106e-06, 2.87268589e-06, 3.22259298e-06, 3.41580409e-06]
# # # for srand 117:
# elif srand == 117:
#     x_new = [3.96048118e-05, 4.17988621e-05, 4.24313151e-05, 3.82888003e-05]
#     sigma_x_new = [3.45073106e-06, 8.61640787e-06, 4.74267423e-06, 4.86217713e-06]
# # # for srand 118:
# elif srand == 118:
#     x_new = [3.96048118e-05, 2.59485331e-05, 2.59911370e-05, 4.29861792e-05]
#     sigma_x_new = [3.45073106e-06, 5.81455293e-06, 3.14624400e-06, 4.86612349e-06]
# # # for srand 119:
# elif srand == 119:
#     x_new = [3.96048118e-05, 2.75541968e-05, 3.60437600e-05, 4.27530424e-05]
#     sigma_x_new = [3.45073106e-06, 2.79952537e-06, 3.81173998e-06, 3.80823584e-06]
# # # for srand 120:
# elif srand == 120:
#     x_new = [3.96048118e-05, 1.48423526e-05, 3.31455994e-05, 3.45750861e-05]
#     sigma_x_new = [3.45073106e-06, 1.51057395e-06, 3.33009449e-06, 3.29166370e-06]
# # # for srand 121:
# elif srand == 121:
#     x_new = [3.96048118e-05, 6.52695462e-05, 6.21996410e-05, 6.17526039e-05]
#     sigma_x_new = [3.45073106e-06, 1.66564710e-05, 9.35761554e-06, 9.73821432e-06]
# # # for srand 122:
# elif srand == 122:
#     x_new = [3.96048118e-05, 2.41081235e-05, 3.07034178e-05, 4.19991037e-05]
#     sigma_x_new = [3.45073106e-06, 1.76479546e-06, 2.58831115e-06, 3.84342166e-06]

# 
# # following data is when in the N = 1 case, the dipole is at the center of the simulation box at 2080
# # for srand 113:
# elif srand == 113:
#     x_new = [3.96048118e-05, 2.71653913e-05, 3.12904084e-05, 3.85422404e-05]
#     sigma_x_new = [3.45073106e-06, 3.43455982e-06, 3.42872923e-06, 3.82909106e-06]
# # # for srand 114:
# elif srand == 114:
#     x_new = [3.96048118e-05, 2.75433704e-05, 3.06398932e-05, 5.01094239e-05]
#     sigma_x_new = [3.45073106e-06, 2.61614237e-06, 3.09724538e-06, 5.02831748e-06]
# # # for srand 116:
# elif srand == 116:
#     x_new = [3.96048118e-05, 3.07689897e-05, 4.05965908e-05, 4.15942661e-05]
#     sigma_x_new = [3.45073106e-06, 2.87268589e-06, 3.22259298e-06, 3.41580409e-06]
# # # for srand 117:
# elif srand == 117:
#     x_new = [3.96048118e-05, 4.17988621e-05, 4.24313151e-05, 3.82888003e-05]
#     sigma_x_new = [3.45073106e-06, 8.61640787e-06, 4.74267423e-06, 4.86217713e-06]
# # # for srand 118:
# elif srand == 118:
#     x_new = [3.96048118e-05, 2.59485331e-05, 2.59911370e-05, 4.29861792e-05]
#     sigma_x_new = [3.45073106e-06, 5.81455293e-06, 3.14624400e-06, 4.86612349e-06]
# # # for srand 119:
# elif srand == 119:
#     x_new = [3.96048118e-05, 2.75541968e-05, 3.60437600e-05, 4.27530424e-05]
#     sigma_x_new = [3.45073106e-06, 2.79952537e-06, 3.81173998e-06, 3.80823584e-06]
# # # for srand 120:
# elif srand == 120:
#     x_new = [3.96048118e-05, 1.48423526e-05, 3.31455994e-05, 3.45750861e-05]
#     sigma_x_new = [3.45073106e-06, 1.51057395e-06, 3.33009449e-06, 3.29166370e-06]
# # # for srand 121:
# elif srand == 121:
#     x_new = [3.96048118e-05, 6.52695462e-05, 6.21996410e-05, 6.17526039e-05]
#     sigma_x_new = [3.45073106e-06, 1.66564710e-05, 9.35761554e-06, 9.73821432e-06]
# # # for srand 122:
# elif srand == 122:
#     x_new = [3.96048118e-05, 2.41081235e-05, 3.07034178e-05, 4.19991037e-05]
#     sigma_x_new = [3.45073106e-06, 1.76479546e-06, 2.58831115e-06, 3.84342166e-06]
# 





x_plot = np.array([1,2,5,10,15,20,25,30,35,40])    
# if srand != 112:
    # x_plot = np.array([1,2,5,10])

# Interpolate y values using the spline
y_new = f(x_new)
dy_dx = f(x_new, 1)                    # the derivative of the function
sigma_y = np.abs(dy_dx) * sigma_x_new  # the error in y

# fitting a line to the higher N values
if srand == 112:
    coefficients = np.polyfit(x_plot[3:], y_new[3:], 1, w = 1./sigma_y[3:])
    if full_flag == 1:
        coefficients = np.polyfit(x_plot, y_new, 1, w = 1./sigma_y)
        
else:
    coefficients = np.polyfit(x_plot, y_new, 1, w = 1./sigma_y)
m, b = coefficients
# Create a line based on the fitted coefficients

if srand == 112:
    y_fit = m * x_plot[3:] + b
    if full_flag == 1:
        y_fit = m * x_plot + b        
else:
    y_fit = m * x_plot + b
    

# Plot the original points and the interpolated curve
# plt.plot(x, y, 'o', label='Original data')
# plt.plot(x_new, y_new, label='Cubic spline')
# plt.legend()
# plt.show()


# plotting p_eff vs N
plt.figure(figsize=(8, 5), dpi=300)
# plt.scatter(x_plot, y_new, marker = '.', s = 500, alpha=0.5, color = 'r')#, label = pbond_string)# Network 1')    
plt.errorbar(x_plot, y_new, yerr = sigma_y, linestyle = 'none', marker = '.', ms = 20, alpha=1, color = 'g')#, label = pbond_string)# Network 1')    
# plt.scatter(x_new, y_new, marker = '.', s = 500, alpha=0.5, color = 'r')#, label = pbond_string)# Network 1')    

# Plot the data and the fitted line
base_num, exp = find_exp_base(m)
# label_txt = r"$p_{eff} = %.2e N_{d} + %.2f $" % (m, b)
label_txt = '$' + str(base_num) + '\\times 10^{'+ str(exp) + '}$' + ' $N_d + $' + r"$%0.2f$" % b
if srand == 112 and full_flag == 0:
    plt.plot(x_plot[3:], y_fit, label=label_txt)
elif srand == 112 and full_flag == 1:
    plt.plot(x_plot, y_fit, label=label_txt)
else:
    plt.plot(x_plot, y_fit, label=label_txt)
# theoretical line
y_theory = 0.54 + (2./517)*x_plot    # 517: inner nodes; 2440: all nodes; 0.54: for N = 0 or 1
# plt.plot(x_plot, y_theory, label='theory')

plt.xlabel("Number of Dipoles", fontsize = 15)
plt.ylabel("$ p_{eff} $", fontsize = 15)
plt.yticks(fontsize = 15)
plt.xticks(fontsize = 15)
plt.subplots_adjust(right=0.96,left=0.18,top=0.94,bottom=0.12)
plt.legend(loc = 'best')
if full_flag == 0:
    plt.savefig(base+emt_folder+"p_eff_vs_N_kappa_e-6_"+ 'srand_'+ str(srand) +"_" +str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+kappa_str+"_"+tol_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+".png")
elif full_flag == 1:
    plt.savefig(base+emt_folder+"p_eff_vs_N_kappa_e-6_all_data_"+ 'srand_'+ str(srand) +"_" +str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+kappa_str+"_"+tol_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+".png")

plt.clf()
plt.close()


from matplotlib.ticker import FuncFormatter
# Define the formatter function
def scientific_format(x, pos):
    if x == 0:
        return "0"
    exponent = int(np.floor(np.log10(abs(x))))
    coefficient = x / 10**exponent
    return r"${:.1f} \times 10^{{{}}}$".format(coefficient, exponent)
# plt.ticklabel_format(axis='y', style='sci', scilimits=(0,0))

# plotting alpha_m vs N
plt.figure(figsize=(8, 5), dpi=300)
plt.errorbar(x_plot, x_new, yerr = sigma_x_new, linestyle = 'none', marker = '.', ms = 20, alpha=1, color = 'g')#, label = pbond_string)# Network 1')    
plt.xlabel("Number of Dipoles", fontsize = 20)
plt.ylabel("$ D_{far,p=0.55}/D_{far, p=1} $", fontsize = 20)
plt.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
plt.yticks(fontsize = 20)
plt.xticks(fontsize = 20)
# plt.yscale('log')
# plt.xscale('log')
formatter = FuncFormatter(scientific_format)
plt.gca().yaxis.set_major_formatter(formatter)
plt.subplots_adjust(right=0.98,left=0.25,top=0.98,bottom=0.15)
# plt.legend(loc = 'best')
plt.savefig(base+emt_folder+"alpha_m_vs_N_kappa_e-6_"+ 'srand_'+ str(srand) +"_" +str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+kappa_str+"_"+tol_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+".png")

plt.clf()
plt.close()

#==============================================================================
# dloc values for p=1 network for srand 112
if srand == 112:
    dloc_list_p1 = np.array([5.7661642e-01, 1.1533254e+00, 2.8838063e+00, 5.7692234e+00, 8.6572705e+00, 1.1546441e+01, 1.4438562e+01, 1.7332349e+01, 2.0227529e+01, 2.3124529e+01])
    dloc_list_p55 = np.array([5.4000025e-01, 1.0800009e+00, 2.7000028e+00, 5.4000048e+00, 8.1000132e+00, 1.0800024e+01, 1.3500034e+01, 1.6200052e+01, 1.8900142e+01, 2.1600217e+01])
    dloc_ratio = dloc_list_p55/dloc_list_p1
    
    plt.figure(figsize=(8, 5), dpi=300)
    plt.scatter(x_plot, dloc_ratio, marker = '.', s = 1200, alpha=1, color = 'g')#, label = pbond_string)# Network 1')    
    plt.xlabel("Number of Dipoles", fontsize = 25)
    plt.ylabel("$ D_{loc, p=0.55}/D_{loc, p=1} $", fontsize = 25)
    plt.yticks(fontsize = 25)
    plt.xticks(fontsize = 25)
    plt.ylim([0.9,1.0])
    # plt.yscale('log')
    # plt.xscale('log')
    plt.subplots_adjust(right=0.96,left=0.18,top=0.94,bottom=0.168)
    # plt.legend(loc = 'best')
    plt.savefig(base+emt_folder+"dloc_ratio_vs_N_kappa_e-6_"+ 'srand_'+ str(srand) +"_" +str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+kappa_str+"_"+tol_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+".png")
    
    plt.clf()
    plt.close()


#==============================================================================
x_test = np.logspace(-6, -4,100)
y_test = f(x_test)
# plt.scatter(x_test, y_test)
plt.scatter(y_test, x_test)
plt.yscale('log')


#==============================================================================
# writing fit to sigma1 or sigma1/alpha_m to file
#==============================================================================
disp_outfname = disp_folder+"fit_radial_disp_"+str(srand)+"_"+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+".txt"
print('here the fits are written to file : ',disp_outfname)    
heading = 'pbond         kappa             N        Sigma1/alpha_m'
fmt = '%12d', '%10.2e', '%12d', '%15.7e'
np.savetxt(disp_outfname, np.column_stack((pbond, kappa, num_center, popt[0])), header = heading, fmt = fmt)


#===================================================================================================================================================
# writing the N_d and corresponding p_eff to file
#===================================================================================================================================================
p_eff_outfname = disp_folder+"peff_vs_Nd_"+str(srand)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+".txt"
print('here effitive p values are written to file : ',p_eff_outfname)    
heading = '# Dipoles         mean p_eff          std p_eff'
fmt = '%12d', '%15.7e', '%15.7e'
np.savetxt(p_eff_outfname, np.column_stack((x_plot, y_new, sigma_y)), header = heading, fmt = fmt)








'''










'''
#==============================================================================
# temp plotting the fit paramters and calculated sigma1 or stress
# also reading the mean inner stress data and plotting
#==============================================================================
# reading stress values from file for the radial dipole in p = 1 network
base = '/home/abhinav/david/'

# pbond = 1.0
pbond = 0.55
pbond_string = "%.2f" % pbond # to write the filename

kappa = 1e-6
kappa_str = "%.2e" % kappa # to write the filename

rlen = 0.9
rlen_txt = "%.4f" % rlen
rlen_txt_ring = "%.9f" % 0.1

srand = 112
num_center_list = [1,2,3,4,5,10,15,20,25,30,35,40]

radial_flag = 1
hex_flag = 0
# for true p = 1 case
if pbond == 1:
    hex_flag = 1
    radial_flag = 0
    
if radial_flag == 1:
    folder = 'lattice_nopbd_bndry_all_clamp_restlength_circle_inner_fixed_only_radial_arp/'
elif hex_flag == 1:
    folder = 'lattice_nopbd_bndry_all_clamp_restlength_circle_inner_fixed_radial_arp/'
        

stress_data = []
data_all_inner_stress = []
for i in range(0,len(num_center_list)):

    if pbond == 1.0:
        num_dip = 6*num_center_list[i]
        stress_bond_fname = base+folder+"txt/force/"+"inner_bndry_bond_radial_force_total_"+str(srand)+"_"+str(num_center_list[i])+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+".txt"
            
        print('here all individual radial forces on all bonds that cross the inner boundary are read from file : ',stress_bond_fname)    
        stress_data.append(np.loadtxt(stress_bond_fname))
        
    if pbond != 1.0:
        stress_data.append([])
        for diff_network in diff_network_list:
            if hex_flag == 1:
                folder = 'cluster_download/lattice_nopbd_bndry_all_clamp_restlength_circle_inner_fixed_radial_arp_bash_hex/'+ kappa_fname +ranseed+"/"+str(srand)+"/"+diff_network+"/"
                plot_folder = base + 'cluster_download/lattice_nopbd_bndry_all_clamp_restlength_circle_inner_fixed_radial_arp_bash_hex/'+ kappa_fname +ranseed+"/"+str(srand)
            elif radial_flag == 1:
                folder = 'cluster_download/lattice_nopbd_bndry_all_clamp_restlength_circle_inner_fixed_radial_arp_bash_radial/'+ kappa_fname +ranseed+"/"+str(srand)+"/"+diff_network+"/"
                plot_folder = base + 'cluster_download/lattice_nopbd_bndry_all_clamp_restlength_circle_inner_fixed_radial_arp_bash_radial/'+ kappa_fname +ranseed+"/"+str(srand)
    
            num_dip = 6*num_center_list[i]
            stress_bond_fname = base+folder+"txt/force/"+"inner_bndry_bond_radial_force_total_"+str(srand)+"_"+str(num_center_list[i])+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+".txt"
                
            # print('here all individual radial forces on all bonds that cross the inner boundary are read from file : ',stress_bond_fname)    
            stress_data[i].append(np.loadtxt(stress_bond_fname))


inner_stress_num_center_list = [1,2,5,10,15,20,25,30, 35, 40]
for i in range(0,len(inner_stress_num_center_list)):
    num_dip = 6*inner_stress_num_center_list[i]
    if pbond == 1.0:
        # inner stress data reading
        fname_inner_stress = base+folder+"txt/area/"+"total_inner_stress_"+str(srand)+"_"+str(inner_stress_num_center_list[i])+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+".txt"    

    if pbond != 1.0:
        # inner stress data reading
        fname_inner_stress = plot_folder+"/means/"+"mean_bend_dom_inner_stress_"+str(inner_stress_num_center_list[i])+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+".txt"    
            
    print('reading inner stress data : ',fname_inner_stress)    
    # inner stress mean data
    data_inner_stress = np.loadtxt(fname_inner_stress)        
    data_all_inner_stress.append(data_inner_stress)

data_all_inner_stress = np.array(data_all_inner_stress)
inner_stress_means = data_all_inner_stress[:,1]    
if pbond == 0.55:
    inner_stress_std = data_all_inner_stress[:,2]    

stress_data = np.array(stress_data)
xaxis = np.array([1,2,3,4,5,10,15,20,25,30,35,40])
if pbond == 1.0:
    sigma1_calc = stress_data[:,1]

    # multiplying by -1 below so that we kep the convention of inner boundary force as positive and vice versa
    # for radial only bonds:
    if radial_flag == 1:
        yaxis = -1 * np.array([-5.2509722e-04, -1.0406139e-03, -1.5573606e-03, -2.0817818e-03, -2.6030781e-03, -5.2585336e-03, -7.7890617e-03, -1.0367808e-02, -1.2904824e-02, -1.5424978e-02, -1.7937837e-02, -2.0466067e-02])
    # for radial and HEX bonds:
    elif hex_flag == 1:
        # yaxis = -1 * np.array([-3.19e-04, -6.35e-04, -9.51e-04, -1.27e-03, -1.59e-03, -3.18e-03, -4.78e-03, -6.37e-03, -7.97e-03, -9.55e-03, -1.11e-02, -1.27e-02])
        yaxis = -1 * np.array([-2.92e-04, -5.81e-04, -8.71e-04, -1.16e-03, -1.46e-03, -2.91e-03, -4.37e-03, -5.83e-03, -7.29e-03, -8.74e-03, -1.02e-02, -1.17e-02])
    fig = plt.figure(figsize=(8, 4), dpi=300)
    plt.scatter(xaxis, yaxis, c = 'g', marker = 'd', label = '$p = 1$')
    # plt.scatter(xaxis, yaxis, c = 'g', marker = 'd', label = 'fitted')
    # plt.scatter(xaxis, sigma1_calc, c = 'b', marker = 'd', label = 'measured')
    plt.xlabel('$ N $', fontsize = 15)
    plt.ylabel(' $ \\Sigma1 / \\alpha_{m}$ ', fontsize = 15)
    plt.xticks(fontsize = 15)
    plt.yticks(fontsize = 15)
    # plt.yscale('log')
    plt.subplots_adjust(top=0.92, left = 0.17, bottom = 0.15, right = 0.98)
    plt.legend()
    if radial_flag == 1:
        plt.savefig(base+emt_folder+"sigma1_alpha_m_fit_p1_radial.png")  # just plotting the last step data
    elif hex_flag == 1:
        plt.savefig(base+emt_folder+"sigma1_alpha_m_fit_p1_hex.png")  # just plotting the last step data
    plt.close()
    plt.clf()

# # p = 0.55
if pbond == 0.55:
    sigma1_calc = stress_data[:,:,1]    
    mask = np.zeros((np.shape(sigma1_calc)[0], np.shape(sigma1_calc)[1]))
    mask.fill(False)
    mask[-1,[16, 20, 71]] = True         # masking the bad values for N = 40 case
    mask[-2,93] = True
    masked_sigma1_calc = np.ma.masked_array(sigma1_calc, mask)
    masked_sigma1_calc_mean = np.ma.mean(masked_sigma1_calc, axis = 1)
    masked_sigma1_calc_std = np.ma.std(masked_sigma1_calc, axis = 1)
    sigma1_calc_mean = np.mean(sigma1_calc, axis = 1)
    sigma1_calc_std = np.std(sigma1_calc, axis = 1)
    
    sigma1_calc_plot = [6.3678183873e-09, 1.1141360774299997e-08,
                       # 2.2077923490100006e-08, 3.0624675618799996e-08,
                       4.362253635499999e-08, 4.933839134230001e-08,
                       8.640829362300002e-08, 1.5212972245700003e-07,
                       2.3328033657399999e-07, 3.6801390054e-07]
                       # 2.7528307650707067e-06, 1.6886466444432985e-05]
    
    xaxis = np.array([1,2,5,10,15,20,25,30])
    # multiplying by -1 below so that we kep the convention of inner boundary force as positive and vice versa
    yaxis = -1 * np.array([-2.2027909e-03, -1.8331579e-03, -7.1262771e-03, -1.2000501e-02, -1.9262847e-02, -3.2567470e-02, -4.3929266e-02, -5.4770616e-02])
    fig = plt.figure(figsize=(8, 4), dpi=300)
    plt.scatter(xaxis, yaxis, c = 'g', marker = 'd', alpha = 1, label = '$p = 0.55$ $\kappa = 10^{-6}$ ')
    # plt.scatter(xaxis, yaxis, c = 'g', marker = 'd', alpha = 0.3, label = '$p = 0.55$ $\kappa = 10^{-6}$ Fitted')
    # plt.scatter(xaxis, sigma1_calc_plot, c = 'r', marker = 'd', alpha = 0.3, label = '$p = 0.55$ $\kappa = 10^{-6}$ Calculated')
    plt.xlabel('$ N $', fontsize = 15)
    plt.ylabel(' $ \\Sigma1 / \\alpha_{m}$ ', fontsize = 15)
    plt.xticks(fontsize = 15)
    plt.yticks(fontsize = 15)
    # plt.yscale('log')
    plt.subplots_adjust(top=0.92, left = 0.17, bottom = 0.15, right = 0.98)
    plt.legend()
    if radial_flag == 1:
        plt.savefig(base+emt_folder+"sigma1_alpha_m_fit_p_55_radial.png")  # just plotting the last step data
    elif hex_flag == 1:
        plt.savefig(base+emt_folder+"sigma1_alpha_m_fit_p_55_hex.png")  # just plotting the last step data
    plt.close()
    plt.clf()

'''

'''

# plotting p = 1 (hex) and p = 0.55 (radial) cases together
xaxis = [1,2,5,10,15,20,25,30, 35, 40]
yaxis_hex_p1 = np.array([0.000292, 0.000581, 0.00146 , 0.00291 , 0.00437 , 0.00583 , 0.00729 , 0.00874, 0.0102  , 0.0117])
yaxis_hex_p1_norm = yaxis_hex_p1/yaxis_hex_p1[0]
# yaxis_radial_p55 = np.array([0.00220279, 0.00183316, 0.00712628, 0.0120005 , 0.01926285, 0.03256747, 0.04392927, 0.05477062, ])
# yaxis_radial_p55 = np.array([0.00200697, 0.00166728, 0.00649298, 0.01091102 , 0.01753511, 0.02980966, 0.04023195, 0.05011304, 0.06349771, 0.07446804])
# yaxis_radial_p55_err = np.array([8.41e-5, 5.01e-5, 2.12e-4, 3.42e-4, 5.73e-4, 8.67e-4, 1.07e-3, 1.28e-3, 1.54e-3, 1.73e-3, ])

yaxis_radial_p55 = np.array([0.00180373, 0.00149844, 0.00583546, 0.00980611, 0.01575941, 0.02679096, 0.03615783, 0.04503831, 0.05706757, 0.06692698])
yaxis_radial_p55_err = np.array([7.55e-5, 4.5e-5, 1.91e-4, 3.08e-4, 5.15e-4, 7.79e-4, 9.64e-4, 1.15e-3, 1.38e-3, 1.55e-3])

yaxis_radial_p55_norm = yaxis_radial_p55/yaxis_radial_p55[0]
yaxis_radial_p55_err_norm = yaxis_radial_p55_err/yaxis_radial_p55[0]

fig = plt.figure(figsize=(8, 6), dpi=300)
# plt.scatter(xaxis, yaxis_radial_p55, c = 'g', marker = 'o', alpha = 1, s = 150, label = '$p = 0.55$ $\kappa = 10^{-6}$ ')
# plt.scatter(xaxis, yaxis_hex_p1, c = 'b', marker = 'd', alpha = 1, s = 150, label = '$p = 1$')
# plt.scatter(xaxis, yaxis_radial_p55_norm, c = 'g', marker = 'o', alpha = 1, s = 150, label = '$p = 0.55$ $\kappa = 10^{-6}$ ')
plt.errorbar(xaxis, yaxis_radial_p55_norm, yerr = yaxis_radial_p55_err_norm, c = 'g', ls = 'none', marker = 'o', alpha = 1, ms = 10, label = '$p = 0.55$ $\kappa = 10^{-6}$ ')
plt.scatter(xaxis, yaxis_hex_p1_norm, c = 'b', marker = 'd', alpha = 1, s = 150, label = '$p = 1$')
plt.xlabel('$ N $', fontsize = 25)
# plt.ylabel(' $ \\Sigma1 / \\alpha_{m}$ ', fontsize = 25)
plt.ylabel(' Normalized $ \\Sigma_1 / \\mu_{m}$ ', fontsize = 25)
plt.xticks(fontsize = 25)
plt.yticks(fontsize = 25)
# plt.yscale('log')
# plt.xscale('log')
plt.subplots_adjust(top=0.96, left = 0.17, bottom = 0.15, right = 0.98)
plt.legend(loc = 'best', fontsize = 15)
plt.savefig(base+emt_folder+"sigma1_alpha_m_fit_p_55_radial_p1_hex.png")  # just plotting the last step data
plt.close()
plt.clf()

fig = plt.figure(figsize=(8, 6), dpi=300)
# plt.scatter(xaxis, yaxis_radial_p55, c = 'g', marker = 'o', alpha = 1, s = 150, label = '$p = 0.55$ $\kappa = 10^{-6}$ ')
# plt.scatter(xaxis, yaxis_hex_p1, c = 'b', marker = 'd', alpha = 1, s = 150, label = '$p = 1$')
# plt.scatter(xaxis, yaxis_radial_p55_norm, c = 'g', marker = 'o', alpha = 1, s = 150, label = '$p = 0.55$ $\kappa = 10^{-6}$ ')
plt.errorbar(xaxis, yaxis_radial_p55_norm, yerr = yaxis_radial_p55_err_norm, c = 'g', ls = 'none', marker = 'o', alpha = 1, ms = 10, label = '$p = 0.55$ $\kappa = 10^{-6}$ ')
plt.scatter(xaxis, yaxis_hex_p1_norm, c = 'b', marker = 'd', alpha = 1, s = 150, label = '$p = 1$')
plt.xlabel('$ N $', fontsize = 25)
# plt.ylabel(' $ \\Sigma1 / \\alpha_{m}$ ', fontsize = 25)
plt.ylabel(' Normalized $ \\Sigma_1 / \\mu_{m}$ ', fontsize = 25)
plt.xticks(fontsize = 25)
plt.yticks(fontsize = 25)
plt.yscale('log')
plt.xscale('log')
plt.subplots_adjust(top=0.96, left = 0.17, bottom = 0.15, right = 0.98)
plt.legend(loc = 'best', fontsize = 15)
plt.savefig(base+emt_folder+"sigma1_alpha_m_fit_logscale_p_55_radial_p1_hex.png")  # just plotting the last step data
plt.close()
plt.clf()


# plotting the inner stress::
fig = plt.figure(figsize=(8, 6), dpi=300)
if pbond == 0.55:
    plt.errorbar(xaxis, inner_stress_means, yerr = inner_stress_std, c = 'g', ls = 'none', marker = 'o', alpha = 1, ms = 10, label = '$p = 0.55$ $\kappa = 10^{-6}$ ')
elif pbond == 1:
    plt.scatter(xaxis, inner_stress_means, c = 'g', marker = 'o', alpha = 1, s = 100, label = '$p = 1$')
plt.xlabel('$ N_d $', fontsize = 25)
# plt.ylabel(' $ \\Sigma1 / \\alpha_{m}$ ', fontsize = 25)
plt.ylabel(' Inner Stress ', fontsize = 25)
plt.xticks(fontsize = 25)
plt.yticks(fontsize = 25)
plt.yscale('log')
plt.xscale('log')
plt.subplots_adjust(top=0.96, left = 0.17, bottom = 0.15, right = 0.98)
plt.legend(loc = 'best', fontsize = 15)
if pbond == 0.55:
    plt.savefig(base+emt_folder+"calculated_inner_stress_p55_radial.png")  # just plotting the last step data
elif pbond == 1:
    plt.savefig(base+emt_folder+"calculated_inner_stress_p1_hex.png")  # just plotting the last step data
plt.close()
plt.clf()

# '''

'''
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# y1 = np.array([y_new1, y_new2, y_new3, y_new4, y_new5, y_new6, y_new7, y_new8, y_new9])

# new_srand_list = ['113', '114', '116','117','118','119','120','121', '122']                         # 115 is BAD!!!!!
new_srand_list = ['112','113','114', '116','117','118','119','120','121', '122']                         # 115 is BAD!!!!!
all_peff = []

for srand_val in new_srand_list:    
    p_eff_outfname = ranseed_folder + srand_val +"/means/" + "peff_vs_Nd_"+srand_val+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+".txt"
    data_peff = np.loadtxt(p_eff_outfname)
    all_peff.append(data_peff)

all_peff = np.array(all_peff)

# y1 = np.column_stack((y_new1, y_new2, y_new3, y_new4, y_new5, y_new6, y_new7, y_new8, y_new9))
# yerr_all = np.column_stack((yerr1, yerr2, yerr3, yerr4, yerr5, yerr6, yerr7, yerr8, yerr9))

y1 = all_peff[:,:,1]   # this is the peff value
yerr_all = all_peff[:,:,2]
nd_vals = all_peff[0,:,0]

peff_plot = np.mean(y1, axis = 0)
peff_err = np.mean(yerr_all, axis = 0)

coefficients = np.polyfit(nd_vals, peff_plot, 1, w = 1./peff_err)
m, b = coefficients
y_fit_peff = m * nd_vals + b
# Plot the data and the fitted line
base_num, exp = find_exp_base(m)
label_txt = '$' + str(base_num) + '\\times 10^{'+ str(exp) + '}$' + ' $N_d + $' + r"$%0.2f$" % b


fig = plt.figure(figsize=(8, 6), dpi=300)
plt.errorbar(nd_vals, peff_plot, yerr = peff_err, c = 'g', ls = 'none', marker = 'o', alpha = 1, ms = 10, label = '$p = 0.55$ $\kappa = 10^{-6}$ ')
plt.plot(nd_vals, y_fit_peff, label=label_txt)
plt.xlabel('$ N_d $', fontsize = 25)
# plt.ylabel(' $ \\Sigma1 / \\alpha_{m}$ ', fontsize = 25)
plt.ylabel(' $P_{eff}$ ', fontsize = 25)
plt.xticks(fontsize = 25)
plt.yticks(fontsize = 25)
# plt.yscale('log')
# plt.xscale('log')
plt.subplots_adjust(top=0.96, left = 0.17, bottom = 0.15, right = 0.98)
plt.legend(loc = 'best', fontsize = 15)
plt.savefig(base+emt_folder+"peff_srand_113_through_122_pooled.png")  # just plotting the last step data
plt.close()
plt.clf()



'''
# gotta run the mu_m code (sigma1... code) first for the following!

# y_plot = f(mu_m)

# fig = plt.figure(figsize=(8, 6), dpi=300)
# plt.scatter(x_plot, y_plot, s = 150, label='$\widetilde\kappa = 10^{-6}$')
# plt.xlabel('$ N_d $', fontsize = 25)
# plt.ylabel(' $P_{eff}$ ', fontsize = 25)
# plt.xticks(fontsize = 25)
# plt.yticks(fontsize = 25)
# # plt.yscale('log')
# # plt.xscale('log')
# plt.subplots_adjust(top=0.96, left = 0.17, bottom = 0.15, right = 0.98)
# plt.legend(loc = 'best', fontsize = 15)
# plt.savefig(base+emt_folder+"peff_srand112_from_stress.png")  # just plotting the last step data
# plt.close()
# plt.clf()


'''
y_plot2 = f(mu_m2)

fig = plt.figure(figsize=(8, 6), dpi=300)
plt.scatter(x_plot, y_plot2, s = 150, label='$\widetilde\kappa = 10^{-6}$')
plt.xlabel('$ N_d $', fontsize = 25)
plt.ylabel(' $P_{eff}$ ', fontsize = 25)
plt.xticks(fontsize = 25)
plt.yticks(fontsize = 25)
# plt.yscale('log')
# plt.xscale('log')
plt.subplots_adjust(top=0.96, left = 0.17, bottom = 0.15, right = 0.98)
plt.legend(loc = 'best', fontsize = 15)
plt.savefig(base+emt_folder+"peff_srand113_122_from_stress_rings.png")  # just plotting the last step data
plt.close()
plt.clf()


y_plot3 = f(mu_m_stress_fits)   # this is p_eff found from mu_m3 from the simga1_calc... file

fig = plt.figure(figsize=(8, 6), dpi=300)
plt.scatter(x_plot, y_plot3, s = 150, label='$\widetilde\kappa = 10^{-6}$')
plt.xlabel('$ N_d $', fontsize = 25)
plt.ylabel(' $P_{eff}$ ', fontsize = 25)
plt.xticks(fontsize = 25)
plt.yticks(fontsize = 25)
# plt.yscale('log')
# plt.xscale('log')
plt.subplots_adjust(top=0.96, left = 0.17, bottom = 0.15, right = 0.98)
plt.legend(loc = 'best', fontsize = 15)
plt.savefig(base+emt_folder+"peff_srand113_122_from_stress_fits.png")  # just plotting the last step data
plt.close()
plt.clf()

'''