%The number of simulation iterations used to calculate the transition matrix
transition_loops='2000'; 

%The number of simulation iterations used to calculate the Voronoi centers
voronoi_loops='50';  

%The file to store the temporary BNG files
temporary_save_location='/home/grad/Documents/MATLAB/multi_gene_MISA/temp_scratch'; 

%Number of Voronoi centers 
target_bin_number='100' ;

%Target number of replicas per voronoi center
Replicas_per_bins='100'; 

%Example output file name
output_file_name='transition_matrix_schlogl_021318'; 

%The simulation lagtime tau
tau='1'; 

%The number of chemical species to use in the Voronoi center movement calculation
species='1' ;

%location to store the output .mat file
output_save_file='/home/grad/Documents/MATLAB/multi_gene_MISA'; 

%The BioNetGen root
BNG_root='/home/grad/BioNetGen-2.2.6-stable/bin/run_network'; 

%The BioNetGen .net model file location
model_name='/home/grad/Documents/MATLAB/multi_gene_MISA/schlogl.net' ;

%Change the matlab call to your MATLABROOT
General_bin_movement_and_transition_local_linux(transition_loops,voronoi_loops,temporary_save_location,target_bin_number,Replicas_per_bins,output_file_name,tau,species,output_save_file,BNG_root,model_name)

%Outputs a matfile: output_file_name.mat
%Inside the matfile, the output contains:
%replicas: the current replica positions, tags, and replica weights
%transition_matrix_count: the transition matrix calculated by finding the counts of transitions from sampling region i to sampling region during every lagtime tau
%transition_matrix_count: the transition matrix calculated by summing the weight of transitions from sampling region i to sampling region during every lagtime tau
%iteration: the current 0-4 temporary iteration number for temporarily stored replica information
%toc: the current total computational time
%ijk: the current ijk iteration of transition matrix calculation

%Also ouputs a temporary voronoi.txt, transition_matrix_count.txt and transition_matrix.txt file for every 2nd ijk iteration under output_save_file/output_save_name/
%containing the voronoi centeres, current transition matrix calculated from i->j counts, and current transition matrix calcluated from i-> weights


vor_inf='/home/grad/Documents/MATLAB/multi_gene_MISA/schlogl_bins.txt'; 

%General_bin_movement_and_transition_local_linux_voronoi_inp(transition_loops,voronoi_loops,temporary_save_location,target_bin_number,Replicas_per_bins,output_file_name,tau,species,output_save_file,BNG_root,model_name,vor_inf)
