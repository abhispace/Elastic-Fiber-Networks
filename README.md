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
4. python_analysis/bndry_forces_circular.py: this code is as important as the first code above. It calculates the local and far field dipole moments and stores the results in output files. It also calculates stress according to LAMMPS documentation.
5. python_analysis/bndry_forces_circular_caller.py: it automates the calling and running of bndry_forces_circular.py.
6. python_analysis/dipole_moment_plots.py: it writes a file with mean dfar, dloc with their stds (for a given srand, which decides a random dipole placement). that file is read by dipole_moment_means.py.
7. python_analysis/dipole_moment_plots_caller.py: automates the running of dipole_moment_plots.py
8. python_analysis/dipole_moment_means.py: it reads the mean values (that have been written to file using dipole_moment_plots.py) for each kappa and N value and produces master plots for a given srand.
9. python_analysis/dipole_moment_means_srand.py: same as dipole_moment_means.py. It was created as a test code to test the results corresponding to different random dipole placements. For each specific dipole placement, it writes the mean far field dipole moment over all possible Nd values.
10. python_analysis/dipole_moment_means_srand_all_srand.py: finds the mean far field dipole moments for larger system size (L=128 or R2 = 50) and compares the results to that of the original system size. Also creates effective p plots.
11. python_analysis/dipole_moment_means_srand_all_models_srand.py: reads in the mean dfar and energies of all three dipole model results for kappa = 1e-6 and plots them as a function of dipole nmber. One limitation is that it reads the means corresponding to each srand first and then averages all ten srand means. This works out if all srand cases have same number of simulations for each Nd, but may not always be mathematically the best approach. all the results here are better created by all_srand_..._plotter.py.
12. python_analysis/all_srand_all_nd_diff_network.py: it reads in each simulation data and calculates mean dfar and radial displacements for all networks which are bending dominated. the results are written to file which is in turn used by all_srand..._plotter.py to create plots.
13. python_analysis/all_srand_all_nd_diff_network_plotter.py: This code reads the values of dfar and energies written out by all_srand_all_nd_diff_network.py code and reads them to make plots. this is where the plot for the paper is created for dfar vs Nd and en ratio vs Nd!! it also calculates the p_eff of all three models.
14. python_analysis/mean_radial_disp_p1.py: only for p=1 cases, for each Nd value, it reads all srand cases and finds the average radial displacements.
15. python_analysis/emt_func.py: it can be imported by other files to pass the values of mu_m (also called alpha_m) and returns p_eff values.

