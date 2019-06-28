function General_bin_movement_and_transition_local_linux(transition_loops,voronoi_loops,temporary_save_location,target_bin_number,Replicas_per_bins,output_file_name,tau,species,output_save_file,BNG_root,model_name)
    %% Initializing timing, variables and fixing types
    transition_loops = str2double(transition_loops);
    target_bin_number = str2double(target_bin_number);
    Replicas_per_bins = str2double(Replicas_per_bins);
    tau = str2double(tau);
    voronoi_loops = str2double(voronoi_loops);
    species = str2double(species);
    
    % Sets BioNetGen (BNG) output steps to 2 to avoid BNG I/O errors
    timestep = tau/2;
    number_simulation_steps = 2;
    
    % Starting timer
    tic = cputime;

    %% Setting up folders 
    % Initializing temporary and the final save directory.
    save_location = [output_save_file '/' output_file_name '.mat'];
    temp_dir_transition = [output_save_file '/' output_file_name '/'];
    if not(exist(temp_dir_transition, 'dir'))
        mkdir(temp_dir_transition);
    end
    matobj_save_info = matfile(save_location,'Writable',true); 

    BioNetGen_replica_folder = [temporary_save_location]; 
    if not(exist([BioNetGen_replica_folder '/t0'], 'dir'))
        mkdir([BioNetGen_replica_folder '/t0']);
    end

    %% Initializing parallel pool depending on matlab version
    if verLessThan('matlab','8.1')
        % Runs on R2012b and earlier
        cores = 4; 
        if matlabpool('size')==0   
            matlabpool(cores)
        end
    else
        % Runs on R2013a and later
        poolobj = gcp('nocreate'); % If no pool, do not create new one.
        if isempty(poolobj)
            poolsize = 0;
        else
            poolsize = poolobj.NumWorkers;
        end
        if poolsize == 0
            parpool('local');
        end
    end
    
    %% Step 0: Initializing first iteration
    iteration = 1;
    replicas = [];
    initial_weight = 1/(Replicas_per_bins);

    %% Step 1: Simulating Replicas_per_bins: new replicas all initialized at the starting location in the file specified by model_name
    parfor ii = 1:Replicas_per_bins
        % runs the command line for the BNG ssa simulation outputs: model_name.gdat
        % (simulation data), model_name.cdat, temporary simulation data
        % model_name_end.net, which is a BioNetGen.net file used to simulate the next
        % replica iteration
        [~, ~] = unix([BNG_root ' -o ' BioNetGen_replica_folder '/t0/newsim_' num2str(ii) ' -e -p ssa -h $RANDOM --cdat 0 --fdat 0 -g ' ...
            model_name ' ' model_name ' ' num2str(timestep) ' ' num2str(number_simulation_steps) ]);
        current_replica_location = dlmread([BioNetGen_replica_folder '/t0/newsim_' num2str(ii) '.gdat'],'',2,1);
        data_keep = [current_replica_location(end,: )];
        
        % storing replicas as:
        % [replica positions, tag for resimulating the .net file, weight of the replica]
        replicas = vertcat(replicas,[data_keep,ii,initial_weight]);
    end

    % Specify the temporary filename
    previous_BNG_dir = [BioNetGen_replica_folder '/t0'];
    % Normalize the replica weights
    replicas(:,end) = replicas(:,end)./sum(replicas(:,end));

    %% Step 2: create the first iteration of Voronoi centers
    % Randomly picks one Voronoi center from the list of replicas as the starting Voronoi
    % center, then iteratively picks the furthest replica from the chosen Voronoi centers
    % until target_bin_number is reached.
    newVoronoi = zeros(target_bin_number,species);
    temp_reps = replicas;
    new_ind_voronoi = randsample(length(temp_reps(:,1)),1);
    newVoronoi(1,:) = temp_reps(new_ind_voronoi,1:species);
    for i = 2:target_bin_number
        [~,Distances] = knnsearch(newVoronoi(1:i-1,:),temp_reps(:,1:species));
        [~,idxmax] = max(Distances);
        newVoronoi(i,:) = temp_reps(idxmax,1:species);
    end
    VoronoiCenters = newVoronoi;
    
    %% Step 3: Perform WEstep on the replicas
    replicas = WEstep_072116(replicas,VoronoiCenters,Replicas_per_bins,species);

    %% Step 4: Determine voronoi positions
    % Loop for moving the voronoi positions. Transition matrix calculation is disabled for this step
    for ijk = 1:voronoi_loops
        display(['Voronoi iteration ' num2str(ijk) ]);
        % Creates a new temporary file, limited to a max of 4.
        newBNGdir = [BioNetGen_replica_folder '/t' num2str(iteration)];
        mkdir(newBNGdir);
        iteration = iteration + 1;

        % HARDCODED: Cap of 4 full replica data storage file iterations at once.
        if iteration >=4
            iteration = 0;
        end

        % Using matlab's parallel environment to simulate a subsection of tagged replicas
        % from the previous temporary file.
        weights = replicas(:,end);
        parallel_BNG_name = replicas(:,end-1);
        replicas_new = [];
        parfor kk = 1:length(weights)
            % Running BNG simulation
            [~,~] = unix([BNG_root ' -o ' newBNGdir '/newsim_' num2str(kk) ' -e -p ssa -h $RANDOM --cdat 0 --fdat 0 -g ' ...
                previous_BNG_dir '/newsim_' num2str(parallel_BNG_name(kk)) '_end.net ' ...
                previous_BNG_dir '/newsim_' num2str(parallel_BNG_name(kk)) '_end.net ' ...
                num2str(timestep) ' ' num2str(number_simulation_steps)   ]);
            current_replica_location = dlmread([newBNGdir '/newsim_' num2str(kk) '.gdat'],'',2,1);
            data_keep = [current_replica_location(end,:)];
            replicas_new = vertcat(replicas_new,[data_keep,kk,weights(kk)]);
        end

        if ijk > 2
            rmdir(previous_BNG_dir,'s');
        end
        previous_BNG_dir = newBNGdir;
        % The pause statement is needed to keep the matlab script from erroring while
        % folders are deleted.
        pause(1);

        % Moving voronoi centers
        newVoronoi = zeros(target_bin_number,species);
        new_ind_voronoi = randsample(length(replicas_new(:,1)),1);
        newVoronoi(1,:) = replicas_new(new_ind_voronoi,1:species);
        for i = 2:target_bin_number
            [~,Distances] = knnsearch(newVoronoi(1:i-1,:),replicas_new(:,1:species));
            [~,idxmax] = max(Distances);
            newVoronoi(i,:) = replicas_new(idxmax,1:species);
        end
        VoronoiCenters = newVoronoi;
        replicas = WEstep_072116(replicas_new,VoronoiCenters,Replicas_per_bins,species);
        
        % Saving last replica positions in case a restart is needed
        dlmwrite([BioNetGen_replica_folder '/replicas_' num2str(iteration) '.txt'],replicas);
        matobj_save_info.replicas = replicas;
        matobj_save_info.VoronoiCenters = VoronoiCenters;

        % Stores the replica information and current voronoi positions every 2 iterations
        if rem(ijk,2)==0
            % save([save_location ],'VoronoiCenters','replicas','-v7.3')
            save(save_location,'VoronoiCenters','replicas','-v7.3')
            dlmwrite([temp_dir_transition 'voronoi' num2str(ijk) '.txt'],VoronoiCenters );
        end
        toc = cputime-tic; % Stopping timer
        
        % Debug options
        % display(toc);
        % display(ijk)
    end
    % Initializes the transition matrix where each i->j movement is counted
    transition_matrix_count = zeros(target_bin_number); 
    
    % Initializes the transition matrix where the replica weight of each i->j movement is counted
    transition_matrix = zeros(target_bin_number);

    %% Step 5: Calculating the transition matrix
    for ijk = 1:transition_loops
        display(['Iteration ' num2str(ijk) ]);
        newBNGdir = [BioNetGen_replica_folder '/t' num2str(iteration)];
        mkdir(newBNGdir);
        iteration = iteration + 1;
        
        % HARDCODED: Cap of 4 full replica data storage file iterations at once.
        if iteration >=4
            iteration = 0;
        end

        % Finding the starting location of the current replicas
        Previous_bin_locations = knnsearch(VoronoiCenters,replicas(:,1:species));
        weights = replicas(:,end);
        parallel_BNG_name = replicas(:,end-1);
        replicas_new = [];
        parfor kk = 1:length(weights)
            % Running BNG simulation
            [~, ~] = unix([BNG_root ' -o ' newBNGdir '/newsim_' num2str(kk) ' -e -p ssa -h $RANDOM --cdat 0 --fdat 0 -g ' ...
                previous_BNG_dir '/newsim_' num2str(parallel_BNG_name(kk)) '_end.net ' ...
                previous_BNG_dir '/newsim_' num2str(parallel_BNG_name(kk)) '_end.net ' ...
                num2str(timestep) ' ' num2str(number_simulation_steps)  ]);
            current_replica_location = dlmread([newBNGdir '/newsim_' num2str(kk) '.gdat'],'',2,1);
            data_keep = [current_replica_location(end,:)];
            replicas_new = vertcat(replicas_new,[data_keep,kk,weights(kk)]);
        end

        if ijk > 1
            rmdir(previous_BNG_dir,'s');
        end
        previous_BNG_dir = newBNGdir;
        pause(1);

        % Calculates the transition matrix only after simulating the first iteration ijk
        if ijk > 1
            Bins= knnsearch(VoronoiCenters,replicas_new(:,1:species));
            count_rows = [Previous_bin_locations,Bins];
            temp_urows = unique(count_rows,'rows');
            for irows = 1:length(temp_urows(:,1))
                outof = temp_urows(irows,1);
                into = temp_urows(irows,2);
                times_visited = (ismember(count_rows,temp_urows(irows,:),'rows'));
                transition_matrix(sub2ind(size(transition_matrix_count),outof,into)) = transition_matrix(sub2ind(size(transition_matrix_count),outof,into))+sum(replicas_new(times_visited,end));
                transition_matrix_count(sub2ind(size(transition_matrix_count),outof,into)) = transition_matrix_count(sub2ind(size(transition_matrix_count),outof,into))+sum(times_visited);
            end

        end
        replicas = WEstep_072116(replicas_new,VoronoiCenters,Replicas_per_bins,species);
        
        % Saving last replica positions in case a restart is needed
        dlmwrite([BioNetGen_replica_folder '/replicas_' num2str(iteration) '.txt'],replicas);
        matobj_save_info.replicas = replicas;
        matobj_save_info.VoronoiCenters = VoronoiCenters;
        matobj_save_info.transition_matrix_count = transition_matrix_count ;
        matobj_save_info.transition_matrix = transition_matrix ;
        
        % Stores the replica information, current voronoi positions, and transition matrices every 2 iterations
        if rem(ijk,2)==0
            save([save_location ],'VoronoiCenters','replicas','transition_matrix_count','transition_matrix','iteration','toc','ijk','-v7.3')
            dlmwrite([temp_dir_transition 'transition_matrix_count' num2str(ijk) '.txt'],transition_matrix_count );
            dlmwrite([temp_dir_transition 'transition_matrix' num2str(ijk) '.txt'],transition_matrix );
        end
        toc = cputime-tic; % Stopping timer
        
        % Debug options
        %disp(toc)
        %disp(ijk)
    end
    toc = cputime-tic; % Stopping timer

    %% Step 6: Saving final results and calculating shortest lifetimes
    % Saves final replica replica information, current voronoi positions, and transition matrices 
    save([save_location ],'VoronoiCenters','replicas','transition_matrix_count','transition_matrix','iteration','toc','ijk','-v7.3')

    % Calculating slowest lifetimes and printing results
    trans_matrix_temp_var_raw = transition_matrix;
    row_summation_trans_matrix = sum(trans_matrix_temp_var_raw,2);
    idx_zeros_in_trans_matrix = find(row_summation_trans_matrix==0);
    trans_matrix_temp_var_raw(idx_zeros_in_trans_matrix,:) = [];
    trans_matrix_temp_var_raw(:,idx_zeros_in_trans_matrix) = [];
    row_summation_trans_matrix(idx_zeros_in_trans_matrix) = [];
    trans_matrix_tot_row=repmat(row_summation_trans_matrix,1,size(trans_matrix_temp_var_raw,2));
    T=trans_matrix_temp_var_raw./trans_matrix_tot_row; %this is the transition probability 
    [Eigenvectors,Eigenvalues]=eig(transpose(T));
    [sorted_Eigevnalues,ik]=sort((real(diag(Eigenvalues))));
    sorted_Eigenvalues = flipud(sorted_Eigevnalues(1:end));
    lifetimes = -tau./log(flipud(sorted_Eigevnalues(1:end-1)));

    display(['The three slowest lifetimes are: ' num2str(lifetimes(1)) ', ' num2str(lifetimes(2)) ', and ' num2str(lifetimes(3))])
end

%% WEstep wrapper
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


%% WEstep algorithm 
function newreps = WEstep(curr_reps,numreps)  %%%%%error in weight calc
    blah = length(curr_reps(:,1));
    if isempty(blah)
        error('reps_empty!')
    end
    % number of repilcas in corred incremented by 1
    newreps = curr_reps;

    curr_weights = newreps(:,end);
    total_weight = sum(curr_weights);
    idx_remove = [];
    [SortedWeights, SortedInd] = sort(curr_weights);
    min_idx = 2;
    max_idx = length(SortedInd);
    smallest_weight = SortedWeights(1);

    % removes smallest weight outlier
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
    
    % removes largest weight outlier
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

    % Merges weights until the number of replicas = numreps
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
    
    % splits weights until the number of replicas = numreps
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
        newreps(max_idx,:) = [];
        
        shouldbezero2 = (length(newreps(:,end))-numreps);
        if shouldbezero2~=0
            error(' replica number from WEstep after splitting is wrong');
        end
    end

    w5 = sum(newreps(:,end));

    % Checking for errors in the weight calculation.
    if (all(abs([total_weight-w2,w2-w3,w3-w4,w4-w5]))) > 10^(-12)
        abs([total_weight-w2,w2-w3,w3-w4,w4-w5])
        error('Error is in one of the three subsections of WE step')
    elseif length(newreps(:,end)) ~= numreps
        error('Final replica length is wrong');
    end
end
