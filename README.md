# Rare-Event-Sampling-Gene-Networks
This repository contains code for running the weighted ensemble (WE) sampling method with adaptive Voronoi-based discretization on biochemical networks.

### Requirements
- Operating System Requirements/Recommendations
  - Linux (Recommended) or Mac OS X 
  - Windows is not fully implemented but a basic version is included as an example
- Matlab 2013a or newer

There are three versions of the code, one for running the WE method to find a transition matrix of transitions between Voronoi-based sampling regions (General_bin_movement_and_transition_local_linux.m and General_bin_movement_and_transition_local_windows.m), and two for running the WE method to find the transition rate between two states of interest.

### Included

There are three versions of the code (found in the Linux folder): 
1. Running the WE method to find a transition matrix of transitions between Voronoi-based sampling regions 
   - `General_bin_movement_and_transition_local.m` 
2. Running the WE method to find the transition rate between two hypershperes.
   - `General_bin_movement_and_rates_hypersphere_local.m`
3. Running the WE method to find the transition rate between two sets of Voronoi centers.
   - `General_bin_movement_and_rates_voronoi_bins_local.m`

### Instructions

Permissions to read the BioNetGen `run_network` application (under `/BioNetGen_files/`) and any `*.net` files must be properly set.

To run the code, use `submission_script.sh` or `submission_script.m` for linux or windows, respectively.

Inside `submission_script.sh` (or `submission_script.m`), all parameters for running the WE method can be set. 

Comment/uncomment the appropriate lines to select which version (1-3) to run. The default is the transition matrix calculation.

The output of `General_bin_movement_and_transition_local_linux.m` is a .mat file with the prefix specified in the input parameters (see `submission_script.m` or `submission_script.sh`).


### Inside the .mat file are the following variables:

* VoronoiCenters: the current position of all Voronoi centers

* replicas: the current replica positions, tags, and replica weights

* transition_matrix_count: the transition matrix calculated by finding the counts of transitions from sampling region i to sampling region during every lagtime tau

* transition_matrix_count: the transition matrix calculated by summing the weight of transitions from sampling region i to sampling region during every lagtime tau

* iteration: the current 0-4 temporary iteration number for temporarily stored replica information

* toc: the current total computational time

* ijk: the current ijk iteration of transition matrix calculation
