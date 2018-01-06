# Adaptive_WE_public
This repository contains code for running the weighted ensemble (WE) sampling method with adaptive Voronoi-based discretization on biochemical networks.

There are three versions of the code, one for running the WE method to find a transition matrix of transitions between Voronoi-based sampling regions (General_bin_movement_and_transition_local_linux.m and General_bin_movement_and_transition_local_windows.m), and two for running the WE method to find the transition rate between two states of interest.

There are two versions of the code for the WE method in transition rate estimation mode, one which finds the rate between two hyperspheres (General_bin_movement_and_rates_hypersphere_local_linux.m), and one that finds the rate between two sets of Voronoi centers(General_bin_movement_and_rates_voronoi_bins_local_linux.m).

Both versions can be run on Windows or Linux, so long as permissions to read the BioNetGen run_network application (under BioNetGen_files) and .net files are properly set.

To run the code, use submission_script.sh or submission_script.m for linux or windows, respectively.
Inside submission_script.sh (or *.m), all parameters for running the WE method can be set. If you want to use one of the three versions of the code, comment/ uncomment the appropriate lines. The default is set to run the transition matrix calculation.
The WE method codes are currently set to open parallel pools for MATLAB R2012b and earlier versions. The parallel computing toolbox is necessary to use paralell computing. If you are using a more recent version of MATLAB, use the parpool option in comments in the main .m files.

The output of General_bin_movement_and_transition_local_linux.m returns a .mat file with the prefix specified in the input parameters (see submission_script.m or submission_script.sh).

Inside the .mat file are several variables:

VoronoiCenters: the current position of all Voronoi centers

replicas: the current replica positions, tags, and replica weights

transition_matrix_count: the transition matrix calculated by finding the counts of transitions from sampling region i to sampling region during every lagtime tau

transition_matrix_count: the transition matrix calculated by summing the weight of transitions from sampling region i to sampling region during every lagtime tau

iteration: the current 0-4 temporary iteration number for temporarily stored replica information

toc: the current total computational time

ijk: the current ijk iteration of transition matrix calculation
