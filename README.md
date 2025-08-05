# Elastic-Fiber-Networks
This repository contains the codes to generate a disordered triangular network with active internal force dipoles embedded in it. This can be used to model biological fiber networks such as the acto-myosin networks and the cell-in-ECM networks. This code was first written by Dr. David Quint during his PhD to study shear behavior of actin networks. In its current form, it has been adapted by Abhinav Kumar (under the guidance of Dr. Kinjal Dasbiswas) to simulate active internal forces. The rights to the code belong to Abhinav Kumar, Dr. David Quint and Dr. Kinjal Dasbiswas.

The codes use a randomly generated triangular network to create fiber networks, and conjugate gradient for energy minimization after internal forces have been applied using isotropic force dipoles.

C Codes
1. Movies-Working_2D_Lamel_network_force_cont_dip_bndry_clampled_restlength_circle_inner_fixed_radial_arp_bash_radial.c: all the transverse bonds (bonds that connect outer dipole nodes) around every force dipole are absent. This is model 1.
2. Movies-Movies-Working_2D_Lamel_network_force_cont_dip_bndry_clampled_restlength_circle_inner_fixed_radial_arp_bash_rand.c: the transverse bonds (bonds that connect outer dipole nodes) around every force dipole may or may not be present (decided randomly, based on global 'p'). This is model 2.
3. Movies-Working_2D_Lamel_network_force_cont_dip_bndry_clampled_restlength_circle_inner_fixed_radial_arp_bash_hex.c: all the transverse bonds (bonds that connect outer dipole nodes) around every force dipole are present. This is model 3.
  
