%!/bin/sh
%The number of simulation iterations used to calculate the transition matrix
transition_loops='2000';

%The number of simulation iterations used to calculate the Voronoi centers
voronoi_loops='1000';  

%The file to store the temporary BNG files
temporary_save_location='C:\Users\usr\tempfiles'; 

%Number of Voronoi centers 
target_bin_number='100'; 

%Target number of replicas per voronoi center
Replicas_per_bins='50'; 

%Example output file name
output_file_name='transition_matrix_schlogl_tau_10_bin_100_010418.mat'; 

%The simulation lagtime tau
tau='10'; 

%The number of chemical species to use in the Voronoi center movement calculation
species='1'; 

%location to store the output .mat file
output_save_file='C:\Users\usr\we_output_storage'; 

%The BioNetGen root
BNG_root='C:\Users\usr\BioNetGen_files\run_network'; 

%The BioNetGen .net model file location
model_name='C:\Users\usr\BioNetGen_files\models\schlogl.net'; 

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

%______________________________________________________________________________


%Runs the rate calculation with previously definited Voronoi centers and uses a region of interest definition where region1 and region 2 are selections of Voronoi centers

%The number of simulation iterations used to calculate the transition rate
transition_rate_loops='2000';

%The file location specifying the original voronoi centers
original_voronoi_center='voronoi_centers.txt';

%The file location that maps each voronoi center to a region e.g. from 0..10
mapping_centers='map.txt';

%The voronoi centers mapping to region1 are the 1st region of interest
region1='1'; 

%The voronoi centers mapping to region2 are the 2nd region of interest
region2='2'; 

%Change the matlab call to your MATLABROOT and uncomment if you want to start a rate calculation
%General_bin_movement_and_rates_voronoi_bins_local_linux(transition_rate_loops,temporary_save_location,Replicas_per_bins,output_file_name,tau,species,output_save_file,BNG_root,model_name,original_voronoi_center,mapping_centers,region1,region2)

%Outputs a matfile: output_file_name.mat
%Inside the matfile, the output contains:
%replicas: the current replica positions, tags, and replica weights
%rates: a matlab struct file containing the sub_variables: 
%     total_number_of_replicas_transferred: the total number of replicas transferred from 2 -> 1 (in the first index) and 1 -> 2 (in the 2nd index)
%     total_weight_of_replicas_transferred: the total number of replicas transferred from 2 -> 1 (in the first index) and 1 -> 2 (in the 2nd index)
%     total_probability_in_region_i
%     tau

%iteration: the current 0-4 temporary iteration number for temporarily stored replica information
%toc: the current total computational time
%ijk: the current ijk iteration of transition matrix calculation

%Also ouputs a temporary MFPT.txt and rate.txt file for under output_save_file/output_save_name/
%containing the MFPT and rates for every ijk iteration


%______________________________________________________________________________


%Runs the rate calculation with previously definited Voronoi centers and uses a region of interest definition where region1 and region 2 are hyperspheres of radius 'radius'
%around a state specified in 'region1' and 'region2'

%The number of simulation iterations used to calculate the transition rate
%transition_rate_loops='2000' 

%The file location specifying the original voronoi centers
%original_voronoi_center='voronoi_centers.txt'

%The radius of each hypersphere
%radius='3'

%The number of chemical species used to definite the hypersphere centers
%species='2' 

%The center of the hypershpere in ROI 1
%region1='[5,3]' 

%The center of the hypershpere in ROI 2
%region2='[3,5]' 

%Change the matlab call to your MATLABROOT and uncomment if you want to start a rate calculation
%General_bin_movement_and_rates_hypersphere_local_linux(transition_rate_loops,temporary_save_location,Replicas_per_bins,output_file_name,tau,species,output_save_file,BNG_root,model_name,original_voronoi_center,radius,region1,region2)

%Outputs a matfile: output_file_name.mat
%Inside the matfile, the output contains:
%replicas: the current replica positions, tags, and replica weights
%rates: a matlab struct file containing the sub_variables: 
%     total_number_of_replicas_transferred: the total number of replicas transferred from 2 -> 1 (in the first index) and 1 -> 2 (in the 2nd index)
%     total_weight_of_replicas_transferred: the total number of replicas transferred from 2 -> 1 (in the first index) and 1 -> 2 (in the 2nd index)
%     total_probability_in_region_i
%     tau

%iteration: the current 0-4 temporary iteration number for temporarily stored replica information
%toc: the current total computational time
%ijk: the current ijk iteration of transition matrix calculation

%Also ouputs a temporary MFPT.txt and rate.txt file for under output_save_file/output_save_name/
%containing the MFPT and rates for every ijk iteration
