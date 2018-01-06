# Adaptive_WE_public
This repository contains code for running the weighted ensemble (WE) sampling method with adaptive Voronoi-based discretization on biochemical networks.
There are three versions of the code, one for running the WE method to find a transition matrix of transitions between Voronoi-based sampling regions (General_bin_movement_and_transition_local_linux.m and General_bin_movement_and_transition_local_windows.m), and two for running the WE method to find the transition rate between two states of interest.
There are two versions of the code for the WE method in transition rate estimation mode, one which finds the rate between two hyperspheres (General_bin_movement_and_rates_hypersphere_local_linux.m), and one that finds the rate between two sets of Voronoi centers(General_bin_movement_and_rates_voronoi_bins_local_linux.m).
Both versions can be run on Windows or Linux, so long as permissions to read the BioNetGen run_network application (under BioNetGen_files) and .net files are properly set.

To run the code, use submission_script.sh or submission_script.m for linux or windows, respectively.
Inside submission_script.sh (or *.m), all parameters for running the WE method can be set. If you want to use one of the three versions of the code, comment/ uncomment the appropriate lines. The default is set to run the transition matrix calculation.
The WE method codes are currently set to open parallel pools for MATLAB R2012b and earlier versions. The parallel computing toolbox is necessary to use paralell computing. If you are using a more recent version of MATLAB, use the parpool option in comments in the main .m files.
