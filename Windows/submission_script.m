%The number of simulation iterations used to calculate the transition matrix
transition_loops='2000' 

%The number of simulation iterations used to calculate the Voronoi centers
voronoi_loops='1000'  

%The file to store the temporary BNG files
temporary_save_location='/home/usr/tempfiles' 

%Number of Voronoi centers 
target_bin_number='100' 

%Target number of replicas per voronoi center
Replicas_per_bins='50' 

%Example output file name
output_file_name='transition_matrix_schlogl_tau_10_bin_100_010418.mat' 

%The simulation lagtime tau
tau='10' 

%The number of chemical species to use in the Voronoi center movement calculation
species='1' 

%location to store the output .mat file
output_save_file='/home/usr/we_output_storage' 

%The BioNetGen root
BNG_root='/home/usr/BioNetGen_files/run_network' 

%The BioNetGen .net model file location
model_name='/home/usr/BioNetGen_files/models/schlogl.net' 

General_bin_movement_and_transition_local_windows(transition_loops,voronoi_loops,temporary_save_location,target_bin_number,Replicas_per_bins,output_file_name,tau,species,output_save_file,BNG_root,model_name)

