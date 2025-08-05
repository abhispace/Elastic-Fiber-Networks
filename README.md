# Elastic-Fiber-Networks
This repository contains the codes to generate a disordered triangular network with active internal force dipoles embedded in it. This can be used to model biological fiber networks such as the acto-myosin networks and the cell-in-ECM networks. This code was first written by Dr. David Quint during his PhD to study shear behavior of actin networks. In its current form, it has been adapted by Abhinav Kumar (under the guidance of Dr. Kinjal Dasbiswas) to simulate active internal forces. The rights to the code belong to Abhinav Kumar, Dr. David Quint and Dr. Kinjal Dasbiswas.

The codes use a randomly generated triangular network to create fiber networks, and conjugate gradient for energy minimization after internal forces have been applied using isotropic force dipoles.

------------------------------------------------------------------------------------------------------------------------------------------------------------------------
C Codes
------------------------------------------------------------------------------------------------------------------------------------------------------------------------

1. Movies-Working_2D_Lamel_network_force_cont_dip_bndry_clampled_restlength_circle_inner_fixed_radial_arp_bash_radial.c: all the transverse bonds (bonds that connect outer dipole nodes) around every force dipole are absent. This is model 1.
2. Movies-Movies-Working_2D_Lamel_network_force_cont_dip_bndry_clampled_restlength_circle_inner_fixed_radial_arp_bash_rand.c: the transverse bonds (bonds that connect outer dipole nodes) around every force dipole may or may not be present (decided randomly, based on global 'p'). This is model 2.
3. Movies-Working_2D_Lamel_network_force_cont_dip_bndry_clampled_restlength_circle_inner_fixed_radial_arp_bash_hex.c: all the transverse bonds (bonds that connect outer dipole nodes) around every force dipole are present. This is model 3.

------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Compilation and execution
------------------------------------------------------------------------------------------------------------------------------------------------------------------------
gcc c_codename.c -Ofast -o lamel -m64 nrutil_jen.c -lm -std=c99

./lamel > output_log.txt &

------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Bash script
------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Example file: arp_radial_dip_cluster_dip_rand1.sub

Modify this file to suit your system settings and directory structure. It was originally written to run on a cluster. Therefore, the line numbers 3 through 8 may need to be commented out.

------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Python analysis codes
------------------------------------------------------------------------------------------------------------------------------------------------------------------------
1. python_analysis/network_circular.py: this is the first code that should be run to start the analysis process. Edit the directory and folder names to suit your own system. Set different flags to 0 or 1 depending on whether you want plotting or not. This code reads in the simulation results of each time step and plots the network strain plot. It also calculates radial displacements of each node and outputs them.
2. python_analysis/network_circular_caller.py: this code automates the running of network_circular.py. Most of this code is commented out now and only the top few lines that are un-commented are good to use. The bottom commented part is left as a learning tool for the reader.
3. python_analysis/network_circular_caller_caller.py: this code automates the analysis in network_circular_caller.py which has been commented out and is not used. Left for educational purposes.
4. 

