function General_bin_movement_and_rates_hypersphere_local_linux(transition_rate_loops,temporary_save_location,Replicas_per_bins,output_file_name,tau,species,output_save_file,BNG_root,model_name,original_voronoi_center,regionradius,region1,region2)

%%converts string inputs into doubles.
Replicas_per_bins = str2double(Replicas_per_bins);
tau = str2double(tau);
transition_rate_loops = str2double(transition_rate_loops);
species = str2double(species);

tic = cputime;

%%Initializes the temporary directory and the final saave directory.
save_location = [output_save_file '/' output_file_name '.mat'];
temp_dir_transition = [output_save_file '/' output_file_name '/'];
if not(exist(temp_dir_transition))
    mkdir(temp_dir_transition);
end
matobj_save_info = matfile(save_location,'Writable',true); %saving the updated weights of the regions

%%loads the information for input regions of interest for the rate
%%conditions
VoronoiCenters = load(original_voronoi_center); %The original voronoi centers
regionradius = load(regionradius); %Specifics the radius for hypersphere region 1
region1 = str2double(region1); %specifies the states for region 1
region2 = str2double(region2); %specifies the states for region 2
%%%%%%%%%%%%%%%%%%%%%%%%%%
%**************************************************************

%***************************************************************
timestep = tau/2;
number_simulation_steps = 2;
%%Sets the reactions to record via BioNetGen to 2 (to avoid input/output
%%issues with BNG).



%%examples of BNG_root
%BNG_root = '/home/grad/BioNetGen-2.2.6-stable/bin/run_network'; %BionetGen simulation root
%BNG_root = '/dfs2/elread/rxn-share/BioNetGen-2.2.6-stable/bin/run_network';
%BNG_root = 'C:\Users\Maggie\Documents\BioNetGen-2.2.6-stable\bin\run_network.exe';

BioNetGen_replica_folder = [temporary_save_location]; %Creates folder to store the temporary BioNetGen simulations
if not(exist([BioNetGen_replica_folder '/t0']))
    mkdir([BioNetGen_replica_folder '/t0']);
end
%%initializes a temporary subfolder t0 to store the first temporary BNG
%%output files


iteration = 1;
%%MODIFY THIS PORTION OF THE CODE IS YOU ARE USING A MATLAB VERSION MORE RECENT THAN MATLAB R2012b
%___________________________________________________________________________
%%Opens the parallel pool for MATLAB R2012b and earlier versions.
cores = 4; %%my computer only had 4 cores.
if matlabpool('size')==0 %opens 'cores' cores for parallel use if the pool isn't already open.
    matlabpool(cores)
end

%% Opens the parallel pool for versions more recent than MATLAB R2012b

% poolobj = gcp('nocreate'); % If no pool, do not create new one.
% if isempty(poolobj)
%    poolsize = 0;
% else
%    poolsize = poolobj.NumWorkers;
% end

% if poolsize == 0
%    parpool('local')
% end
%___________________________________________________________________________



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% First iteration
replicas = [];
initial_weight = 1/(Replicas_per_bins);

%%Step 1: Simulating Replicas_per_bins new replicas all initialized at the
%%starting location in the file specified by model_name
parfor ii = 1:Replicas_per_bins
    
    %%runs the command line for the BioNetGen ssa simulation
    %%outputs model_name.gdat (simulation data)
    %%model_name.cdat, temporary simulation data
    %%model_name_end.net, which ia a BioNetGen.net file used to simulate
    %%the next replica iteration
    [~, ~] = unix([BNG_root ' -o ' BioNetGen_replica_folder '/t0/newsim_' num2str(ii) ' -e -p ssa -h %RANDOM% --cdat 0 --fdat 0 -g ' ...
        model_name ' ' model_name ' ' num2str(timestep) ' ' num2str(number_simulation_steps) ]);
    current_replica_location = dlmread([BioNetGen_replica_folder '/t0/newsim_' num2str(ii) '.gdat'],'',2,1);
    
    data_keep = [current_replica_location(end,: )];
    %%stores the replicas as [replica positions, tag for resimulating the
    %%.net file, weight of the replica]
    replicas = vertcat(replicas,[data_keep,ii,initial_weight]);
    
end


%Specifies the temporary filename
previous_BNG_dir = [BioNetGen_replica_folder '/t0'];
%Normalizes the replica weights


%specifies if a given replica was previous in region1 (denoted as 1) or in
%region2 (denoted as -1)

distsfromROI1 = pdist2(region1,replicas(:,1:species));
distsfromROI2 = pdist2(region2,replicas(:,1:species));
inROI2=find(distsfromROI2<=regionradius);
InROI1=find(distsfromROI1<=regionradius);

temporary_region1or2 = zeros(length(replicas(:,end)),1);
temporary_region1or2(InROI1) = 1;
temporary_region1or2(inROI2) = -1;
replicas = [replicas(:,1:end-2),temporary_region1or2,replicas(:,end-1:end)];
replicas(:,end) = replicas(:,end)./sum(replicas(:,end));


%performs the WEstep on the replicas
replicas = WEstep_072116(replicas,VoronoiCenters,Replicas_per_bins,species);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Stores transition rate information in the struct variable 'rates'
for i = 1:2
    rates(i).total_number_of_replicas_transferred = zeros(transition_rate_loops,1);
    rates(i).total_weight_of_replicas_transferred = zeros(transition_rate_loops,1);
    rates(i).total_probability_in_region_i = zeros(transition_rate_loops,1);
    rates(i).tau = tau;
end

update_MFPT = zeros(transition_rate_loops,2);
%%simulation loop for calculating the transition matrix
for ijk = 1:transition_rate_loops
    
    
    %%definitions are the same here as in the Voronoi center calculation.
    
    newBNGdir = [BioNetGen_replica_folder '/t' num2str(iteration)];
    mkdir(newBNGdir);
    iteration = iteration + 1;
    %%%HARDCODED full replica data storage file iterations. Only four iterations are allowed at a time
    if iteration >=4
        iteration = 0;
    end
    
    %%Finds the starting location of the current replicas
    Previous_ROI = replicas(:,end-2);
    weights = replicas(:,end);
    parallel_BNG_name = replicas(:,end-1);
    replicas_new = [];
    parfor kk = 1:length(weights)
        [~, ~] = unix([BNG_root ' -o ' newBNGdir '/newsim_' num2str(kk) ' -e -p ssa -h $RANDOM --cdat 0 --fdat 0 -g ' ...
            previous_BNG_dir '/newsim_' num2str(parallel_BNG_name(kk)) '_end.net ' ...
            previous_BNG_dir '/newsim_' num2str(parallel_BNG_name(kk)) '_end.net ' ...
            num2str(timestep) ' ' num2str(number_simulation_steps)  ]);
        current_replica_location = dlmread([newBNGdir '/newsim_' num2str(kk) '.gdat'],'',2,1);
        
        data_keep = [current_replica_location(end,:)];
        
        replicas_new = vertcat(replicas_new,[data_keep,Previous_ROI,kk,weights(kk)]);
    end
    
    
    if ijk > 1
        rmdir(previous_BNG_dir,'s');
    end
    previous_BNG_dir = newBNGdir;
    pause(1);
    
    
    
    %%calculates the rates
    
    %%finds which replicas have changed from region 1 to 2 and vice-versa.
    
    distsfromROI1 = pdist2(region1,replicas_new(:,1:species));
    distsfromROI2 = pdist2(region2,replicas_new(:,1:species));
    inROI2=find(distsfromROI2<=regionradius);
    InROI1=find(distsfromROI1<=regionradius);
    Current_ROI = zeros(length(replicas_new(:,1)),1);
    Current_ROI(InROI1') = 1;
    Current_ROI(inROI2') = -1;
    switch1to2 = find(Current_ROI-Previous_ROI==2);
    switch2to1 = find(Current_ROI-Previous_ROI==-2);
    
    %%updates the replicas with their new regions
    replicas_new(InROI1',end-2) = 1;
    replicas_new(InROI2',end-2) = -1;
    
    %%updates the stored rate information informaiton
    rates(1).total_probability_in_region_i(ijk,1) = sum(replicas_new(replicas_new(:,end-2)==1,end));
    rates(2).total_probability_in_region_i(ijk,1) = sum(replicas_new(replicas_new(:,end-2)==-1,end));
    
    [sa,~] = size(switch1to2);
    [sb,~] = size(switch2to1);
    
    prob_rep_A = replicas_new(switch2to1,end);
    rates(1).total_number_of_replicas_transferred(ijk,1) = sb;
    rates(1).total_weight_of_replicas_transferred(ijk,1) = sum(prob_rep_A);
    
    
    prob_rep_B = replicas_new(switch1to2,end);
    rates(2).total_number_of_replicas_transferred(ijk,1) = sa;
    rates(2).total_weight_of_replicas_transferred(ijk,1) =sum(prob_rep_B);
    spec_indA = find(replicas_new(:,end-2)== 1);
    spec_indB = find(replicas_new(:,end-2)== -1);
    spec_indC = find(replicas_new(:,end-2)==0);
    
    
    
    replicas1 = WEstep_072116(replicas_new(spec_indA,:),VoronoiCenters,Replicas_per_bins,species);
    replicas2 = WEstep_072116(replicas_new(spec_indB,:),VoronoiCenters,Replicas_per_bins,species);
    replicas3 = WEstep_072116(replicas_new(spec_indC,:),VoronoiCenters,Replicas_per_bins,species);
    
    replicas = [replicas1;replicas2;replicas3];
    
    
    
    %%keeping the last replica positions in case a restart is necessary
    dlmwrite([BioNetGen_replica_folder '/replicas_' num2str(iteration) '.txt'],replicas);
    
    
    
    matobj_save_info.replicas = replicas;
    matobj_save_info.VoronoiCenters = VoronoiCenters;
    matobj_save_info.rates = rates;
    if rem(ijk,2)==0
        save([save_location ],'VoronoiCenters','replicas','rates','iteration','toc','ijk','-v7.3')
    end
    
    
    for i = 1:2
        update_MFPT(ijk,i) = 1/(rates(i).total_weight_of_replicas_transferred(ijk,1)./(rates(i).total_probability_in_region_i(ijk,1))./(tau));
        
    end
    
    
    %creates a file updated with the current MFPT calculatoin
    if ijk == 1
        dlmwrite([temp_dir_transition 'MFPT_per_ijk_iteration.txt'],update_MFPT(ijk,:) );
        dlmwrite([temp_dir_transition 'rate_per_ijk_iteration.txt'],1./update_MFPT(ijk,:) );
        
    else
        dlmwrite([temp_dir_transition 'MFPT_per_ijk_iteration.txt'],update_MFPT(ijk,:) ,'-append');
        
        dlmwrite([temp_dir_transition 'rate_per_ijk_iteration.txt'],1./update_MFPT(ijk,:) ,'-append');
    end
    
    
    
    
    toc = cputime-tic;
    %disp(toc)
    %disp(ijk)
    
end



toc = cputime-tic;

save([save_location ],'VoronoiCenters','replicas','rates','iteration','toc','ijk','-v7.3')

end


function [newreps] = WEstep_072116(replicas,Voronoi_List,numreps,species)
newreps = [];
bins = knnsearch(Voronoi_List,replicas(:,1:species));
%This code assumes the columns of the replicas matrix are the species, and that the weights are in the last column
parfor i = 1:length(Voronoi_List(:,1))
    replicas_in_bini = replicas(bins==i,:);
    if not(isempty(replicas_in_bini))
        newreps_temp = WEstep(replicas_in_bini,numreps);
        if  abs(sum(replicas_in_bini(:,end)) -   sum(newreps_temp(:,end)) > 10^(-12))
            sum(newreps_temp(:,end))
            sum(replicas_in_bini(:,end))
            error('WEstep failure');
        end
        newreps = [newreps; newreps_temp];
        
    end
end


end



function newreps = WEstep(curr_reps,numreps)  %%%%%error in weight calc

blah = length(curr_reps(:,1));
if isempty(blah)
    error('reps_empty!')
end
%%%%%number of repilcas incorred incremented by 1
newreps = curr_reps;

curr_weights = newreps(:,end);
total_weight = sum(curr_weights);
idx_remove = [];
[SortedWeights, SortedInd] = sort(curr_weights);
min_idx = 2;
max_idx = length(SortedInd);
smallest_weight = SortedWeights(1);

%%removes smallest weight outlier
while (smallest_weight < total_weight/(3*numreps)) && (min_idx < max_idx)
    smallest_weight = sum(SortedWeights(1:min_idx));
    idx_remove = [ SortedInd(1:min_idx)];
    
    min_idx = min_idx+1;
end
if isempty(idx_remove) == 0
    new_rep_ind = randsample(SortedInd(1:min_idx),1,true,SortedWeights(1:min_idx));
    
    new_reps_temp = newreps(new_rep_ind,:);
    new_reps_temp(:,end) = smallest_weight;
    newreps = [newreps; new_reps_temp];
    newreps(idx_remove,:) = [];
end

w2 = sum(newreps(:,end));

curr_weights =  newreps(:,end);
%%removes largest weight outlier
[largest_weight, max_idx] = max(curr_weights);
tmpweight = largest_weight;
mxx = 0;
max_factor = 1;

if largest_weight > total_weight/numreps*3
    
    while largest_weight > total_weight/numreps*3
        largest_weight = tmpweight/max_factor;
        max_factor = max_factor+1;
    end
    
    new_reps_temp = repmat(newreps(max_idx,:),max_factor+mxx,1);
    new_reps_temp(:,end) = repmat(tmpweight/(max_factor+mxx),max_factor+mxx,1);
    newreps = [newreps; new_reps_temp];
    newreps(max_idx,:) = [];
end


curr_weights = newreps(:,end);
curr_len = length(curr_weights);
w3 = sum(newreps(:,end));

%%Merges weights until the number of replicas = numreps

[SortedWeights, SortedInd] = sort(curr_weights);
min_idx = 2;
max_idx = length(SortedInd);
smallest_weight = SortedWeights(1);
idx_remove = [];

while curr_len  > numreps && (min_idx < max_idx)
    smallest_weight = sum(SortedWeights(1:min_idx));
    idx_remove = [ SortedInd(1:min_idx)];
    min_idx = min_idx+1;
    curr_len = curr_len-1;
end
if isempty(idx_remove) == 0
    new_rep_ind = randsample(SortedInd(1:min_idx),1,true,SortedWeights(1:min_idx));
    
    new_reps_temp = newreps(new_rep_ind,:);
    new_reps_temp(:,end) = smallest_weight;
    newreps = [newreps; new_reps_temp];
    newreps(sort(idx_remove),:) = [];
    
end



w4 = sum(newreps(:,end));

%%splits weights until the number of replicas = numreps

curr_weights =  newreps(:,end);
curr_len = length(curr_weights);

[largest_weight, max_idx] = max(curr_weights);
tmpweight = largest_weight;
mxx_fact2 = 0;

max_factor = 1;


if curr_len < numreps
    while curr_len < numreps
        max_factor = max_factor+1;
        curr_len = curr_len+1;
    end
    new_reps_temp = repmat(newreps(max_idx,:),max_factor+mxx,1);
    a = length(new_reps_temp(:,1));
    
    new_reps_temp(:,end) = tmpweight/(max_factor+mxx_fact2).*ones(a,1);
    newreps = [newreps; new_reps_temp];
    shouldbezero = (length(newreps(:,end))-1-numreps);
    
    
    newreps(max_idx,:) = [];
    
    
    shouldbezero2 = (length(newreps(:,end))-numreps);
    
    
    if shouldbezero2~=0
        error(' replica number from WEstep after splitting is wrong');
    end
    
    
end

w5 = sum(newreps(:,end));


%%Checking for errors in the weight calculation.

if (all(abs([total_weight-w2,w2-w3,w3-w4,w4-w5]))) > 10^(-12)
    abs([total_weight-w2,w2-w3,w3-w4,w4-w5])
    error('Error is in one of the three subsections of WE step')
elseif length(newreps(:,end)) ~= numreps
    error('Final replica length is wrong');
end


end



