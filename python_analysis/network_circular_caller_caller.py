#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 31 13:20:04 2025

@author: abhinav
"""
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
srand_list = ['112','113', '114','116','117','118','119','120','121','122']                         # 115 is BAD!!!!!
# srand_list = ['112']
# Run the called script with arguments
for srand_val in srand_list:
        for num_center in num_center_temp:
            subprocess.run(['python', 'network_circular_caller.py', num_center, srand_val])
