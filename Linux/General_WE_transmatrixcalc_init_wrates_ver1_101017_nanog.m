function General_WE_transmatrixcalc_restart_wrates_ver1_092317(vorloops,loops2,temp_save_location,num_nodes,ReplicasRequired,name,tau,tau_name,final_save_location,curr_iteration,species,bash_clump,rstart,full_restart)

loops2 = str2double(loops2);
num_nodes = str2double(num_nodes);
ReplicasRequired = str2double(ReplicasRequired);
tau = str2double(tau);
species = str2double(species);
bash_clump = str2double(bash_clump);
vorloops = str2double(vorloops);
tic = cputime;
curr_iteration = str2double(curr_iteration);
rstart = str2double(rstart);
full_restart = str2double(full_restart);

save_location = [final_save_location '/trans_' name '_tau' tau_name '.mat'];
temp_dir_transition = [final_save_location '/trans_' name '_tau' tau_name '/'];
if not(exist(temp_dir_transition))
    mkdir(temp_dir_transition);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%

VoronoiLocs = load('state.txt');
MAP = load('MAP.txt');


timestep = tau/2;

matobj_save_info = matfile(save_location,'Writable',true); %saving the updated weights of the regions
BioNetGen_replica_folder = temp_save_location; %Creates folder to store the temporary BioNetGen simulations
        newBNGdir = [BioNetGen_replica_folder '/t' num2str(0)];
        if ~exist(newBNGdir,'dir')
            mkdir(newBNGdir);
        end

cd(temp_save_location);


  

    dlmwrite([temp_save_location '/curr_rep_ind.txt'],[1:1:(ReplicasRequired*100)]','precision','%i');
    dlmwrite([temp_save_location '/out_rep_data_ind.txt'],[1:1:(ReplicasRequired*100)]','precision','%i');
    dlmwrite([temp_save_location '/out_rep_data_weight.txt'],ones(ReplicasRequired*100,1),'precision','%10.6e');


unix(['rm -rf ' temp_save_location '/tempout_*.txt']);

    unix([ ' > ' temp_save_location '/curr_out_rep_' num2str(1) '.txt']);
unix([ ' > ' temp_save_location '/finished.txt']);

    iters = floor(ReplicasRequired*100/bash_clump);

    temp= regexp(fileread('sub_BNG_initialization.sh'),'\n','split');
    temp{17} = sprintf('OUT_NUM="%d"',1);
    temp{18} = sprintf('TMP="%s"',temp_save_location);
    temp{19} = sprintf('TSTEP="%f"',timestep);
    temp{21} = sprintf('PREV_ID="%d"',0);
    temp{22} = sprintf('NEW_ID="%d"',0);
    temp{24} = sprintf('STOP="%d"',bash_clump);
    temp{25} = sprintf('PREV_FILE="BNG_net_files/%s"','nanog.net');

    temp{32} = sprintf('#$ -t 1-%d',iters);

    
    
    
    
    
    fid = fopen(['sub_BNG_scripts.sh'],'w');
    fprintf(fid,'%s \n',temp{:});
    fclose(fid);
    
    unix(['qsub sub_BNG_scripts.sh']);
     [~,Is_done] = unix(['wc -l < finished.txt']);

  while str2num(Is_done) < iters
     pause(10)
     [~,Is_done] = unix(['wc -l < finished.txt']);
     disp(str2num(Is_done));
  end




unix(['cat ' temp_save_location '/tempout_*.txt >> ' temp_save_location '/curr_out_rep_' num2str(1) '.txt']);
unix(['rm -rf ' temp_save_location '/finished.txt']);


   

    replicas_forwards = load([temp_save_location '/curr_out_rep_' num2str(1) '.txt']);
    [~,sortind,~] = unique(replicas_forwards(:,end-1));
    replicas_forwards = [replicas_forwards(sortind,2:end)];

replicas_forwards(:,end) = replicas_forwards(:,end)./sum(replicas_forwards(:,end));
    

bins = knnsearch(VoronoiLocs, replicas_forwards(1:species));
inA = find(MAP(bins)==2);
inB = find(MAP(bins)==3);
temploc = zeros(length(replicas_forwards(:,end)),1);
temploc(inA) = 1;
temploc(inB) = -1;
replicas_forwards = [replicas_forwards(:,1:end-2),temploc,replicas_forwards(:,end-1:end)];

    dlmwrite([temp_save_location '/curr_iteration.txt'],0,'precision','%i');
        dlmwrite([temp_save_location '/replica_reserve' num2str(0) '.txt'],replicas_forwards,'precision','%10.10e');




if rstart+1 < vorloops
    for ijk = rstart+1:vorloops
        
        
        
        
        newBNGdir = [BioNetGen_replica_folder '/t' num2str(iteration)];
        if ~exist(newBNGdir,'dir')
            mkdir(newBNGdir);
        end
        
        Binsprevx = knnsearch(VoronoiLocs(:,1:species),replicas_forwards(:,1:species));
        
        dlmwrite([temp_save_location '/curr_rep_ind.txt'],replicas_forwards(:,end-1),'precision','%i');
        dlmwrite([temp_save_location '/out_rep_data_ind.txt'],[1:1:length(Binsprevx)]','precision','%i');
        dlmwrite([temp_save_location '/out_rep_data_weight.txt'],replicas_forwards(:,end),'precision','%10.6e');
        
        weights = replicas_forwards(:,end);
        unix(['rm -rf ' temp_save_location '/tempout_*.txt']);
        
        unix([ ' > ' temp_save_location '/curr_out_rep_' num2str(ijk) '.txt']);
        unix([ ' > ' temp_save_location '/finished.txt']);
        
        iters = floor(length(Binsprevx)/bash_clump);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%
        temp= regexp(fileread('sub_BNG_batch_template.sh'),'\n','split');
        temp{17} = sprintf('OUT_NUM="%d"',ijk);
        temp{18} = sprintf('TMP="%s"',temp_save_location);
        temp{19} = sprintf('TSTEP="%f"',timestep);
        temp{21} = sprintf('PREV_ID="%d"',curr_iteration);
        temp{22} = sprintf('NEW_ID="%d"',iteration);
        temp{24} = sprintf('STOP="%d"',bash_clump);
        temp{32} = sprintf('#$ -t 1-%d',iters);
        
        
        
        
        
        fid = fopen(['sub_BNG_scripts.sh'],'w');
        fprintf(fid,'%s \n',temp{:});
        fclose(fid);
        
        
        
        
        if rem(length(Binsprevx),bash_clump) ~=0
            temp= regexp(fileread('sub_BNG_batch_template.sh'),'\n','split');
            temp{4}= sprintf('#$ -N temp_fin');
            temp{17} = sprintf('OUT_NUM="%d"',ijk);
            temp{18} = sprintf('TMP="%s"',temp_save_location);
            temp{19} = sprintf('TSTEP="%f"',timestep);
            temp{21} = sprintf('PREV_ID="%d"',curr_iteration);
            temp{22} = sprintf('NEW_ID="%d"',iteration);
            temp{24} = sprintf('STOP="%d"',bash_clump);
            
            temp{32} = sprintf('#temp');
            temp{34} = sprintf('for i in {%d..%d}',iters*bash_clump+1,length(Binsprevx));
            temp{36} = sprintf('	CURRLINE="${i}"');
            temp{46} = sprintf('echo 1 > ${TMP}/finished2.txt');
            
            
            fid = fopen(['sub_BNG_scripts_fin.sh'],'w');
            fprintf(fid,'%s \n',temp{:});
            fclose(fid);
            unix(['/opt/gridengine/bin/lx-amd64/qsub ' temp_save_location '/sub_BNG_scripts_fin.sh']);
            
        end
        
        unix(['/opt/gridengine/bin/lx-amd64/qsub ' temp_save_location '/sub_BNG_scripts.sh']);
        [~,Is_done] = unix(['wc -l < finished.txt']);
        
        while str2num(Is_done) < iters
            pause(10)
            [~,Is_done] = unix(['wc -l < finished.txt']);
            disp(str2num(Is_done));
        end
        
        
        unix(['rm -rf ' temp_save_location '/finished.txt']);
        
        if  rem(length(Binsprevx),bash_clump) ~=0
            
            while ~exist([temp_save_location '/finished2.txt']);
                pause(10)
                [~,Is_done] = unix(['wc -l < curr_out_rep_' num2str(ijk) '.txt']);
                disp(str2num(Is_done));
            end
            unix(['rm -rf ' temp_save_location '/finished2.txt']);
            
        end
        
        
        unix(['cat ' temp_save_location '/tempout_*.txt >> ' temp_save_location '/curr_out_rep_' num2str(ijk) '.txt']);
        unix(['rm -rf ' temp_save_location '/finished.txt']);
        
        
        replicas_new = load(['curr_out_rep_' num2str(ijk) '.txt']);
        
        [~,sortind,~] = unique(replicas_new(:,end-1));
        replicas_new = [replicas_new(sortind,2:end-1),weights];
        curr_iteration = iteration;
        iteration = iteration+1;
        if iteration > 3
            iteration = 0;
            unix('rm -rf temp.o*');
            unix('rm -rf temp_fin.o*');
            
        end
        
        newVoronoi = zeros(num_nodes,species);
        disp(size(replicas_new))
        disp(size(newVoronoi))
        new_ind_voronoi = randsample(length(replicas_new(:,1)),1);
        disp(size(new_ind_voronoi))
        newVoronoi(1,:) = replicas_new(new_ind_voronoi,1:species);
        for i = 2:num_nodes
            [~,Distances] = knnsearch(newVoronoi(1:i-1,:),replicas_new(:,1:species));
            [~,idxmax] = max(Distances);
            newVoronoi(i,:) = replicas_new(idxmax,1:species);
        end
        
        %
        VoronoiLocs = newVoronoi;
        
        Binsprev = knnsearch(VoronoiLocs,replicas_new(:,1:species));
        
        
        replicas_forwards = WEstep_051117(replicas_new,VoronoiLocs,ReplicasRequired,Binsprev);
        
        toc = cputime-tic;
        display([ijk,toc]);
        
        matobj_save_info.replicas_forwards = replicas_forwards;
        matobj_save_info.VoronoiLocs = VoronoiLocs;
        matobj_save_info.iteration = iteration;
        save([save_location ],'VoronoiLocs','replicas_forwards','iteration','toc','-v7.3')
        
        
        if rem(ijk,2)==0
            dlmwrite([temp_dir_transition 'voronoi' num2str(ijk) '.txt'],VoronoiLocs );
            
        end
        
        
        
        
        
        
        
        
        
        
    end
end
MFPT_restart = load([temp_save_location '/restart.txt']);

if full_restart == 1



    for i = 1:2
        MFPT(i).counts = zeros(2000,1);
        MFPT(i).prob = zeros(2000,1);
        MFPT(i).weights = zeros(2000,1);
        MFPT(i).replica_number = zeros(2000,1);
        
    end

    MFPT_restart = 1;
end


for ijk = MFPT_restart:MFPT_restart+loops2
    
    if ijk == MFPT_restart
        curr_iteration = load([temp_save_location '/curr_iteration.txt']);
        iteration = curr_iteration + 1;
        replicas_forwards = load([temp_save_location '/replica_reserve' num2str(curr_iteration) '.txt']);
        newBNGdir = [BioNetGen_replica_folder '/t' num2str(iteration)];
        if ~exist(newBNGdir,'dir')
            mkdir(newBNGdir);
        end
        
        binsr = replicas_forwards(:,end-2);
        
    else
        
        
        newBNGdir = [BioNetGen_replica_folder '/t' num2str(iteration)];
        if ~exist(newBNGdir,'dir')
            mkdir(newBNGdir);
        end
        dlmwrite([temp_save_location '/replica_reserve' num2str(curr_iteration) '.txt'],replicas_forwards,'precision','%10.8e');
        

    end
    
    binsr = replicas_forwards(:,end-2);

    
    dlmwrite([temp_save_location '/restart.txt'],ijk+1,'precision','%i');

    dlmwrite([temp_save_location '/out_rep_AtoB.txt'],replicas_forwards(:,end-2),'precision','%i');
    dlmwrite([temp_save_location '/curr_iteration.txt'],curr_iteration,'precision','%i');
    
    dlmwrite([temp_save_location '/curr_rep_ind.txt'],replicas_forwards(:,end-1),'precision','%i');
    dlmwrite([temp_save_location '/out_rep_data_ind.txt'],[1:1:length(replicas_forwards(:,1))]','precision','%i');
    dlmwrite([temp_save_location '/out_rep_data_weight.txt'],replicas_forwards(:,end),'precision','%10.6e');
    weights = replicas_forwards(:,end);
    unix([ ' > ' temp_save_location '/curr_out_rep_' num2str(ijk) '.txt']);
    unix([ ' > ' temp_save_location '/finished.txt']);
    unix([ ' > ' temp_save_location '/finished2.txt']);
    
    iters = floor(length(replicas_forwards(:,1))/bash_clump);
    unix(['rm -rf ' temp_save_location '/tempout_*.txt']);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%
    temp= regexp(fileread('sub_BNG_batch_template_AtoB.sh'),'\n','split');
    temp{17} = sprintf('OUT_NUM="%d"',ijk);
    temp{18} = sprintf('TMP="%s"',temp_save_location);
    temp{19} = sprintf('TSTEP="%f"',timestep);
    temp{21} = sprintf('PREV_ID="%d"',curr_iteration);
    temp{22} = sprintf('NEW_ID="%d"',iteration);
    temp{24} = sprintf('STOP="%d"',bash_clump);
    temp{32} = sprintf('#$ -t 1-%d',iters);
    
    
    
    
    
    fid = fopen(['sub_BNG_scripts.sh'],'w');
    fprintf(fid,'%s \n',temp{:});
    fclose(fid);
    
    
    if rem(length(replicas_forwards(:,1)),bash_clump) ~=0
        temp= regexp(fileread('sub_BNG_batch_template_AtoB.sh'),'\n','split');
        temp{4}= sprintf('#$ -N temp_fin');
        temp{17} = sprintf('OUT_NUM="%d"',ijk);
        temp{18} = sprintf('TMP="%s"',temp_save_location);
        temp{19} = sprintf('TSTEP="%f"',timestep);
        temp{21} = sprintf('PREV_ID="%d"',curr_iteration);
        temp{22} = sprintf('NEW_ID="%d"',iteration);
        temp{24} = sprintf('STOP="%d"',bash_clump);
        
        temp{32} = sprintf('#temp');
        temp{34} = sprintf('for i in {%d..%d}',iters*bash_clump+1,length(replicas_forwards(:,1)));
        temp{36} = sprintf('	CURRLINE="${i}"');
        temp{46} = sprintf('echo 1 > ${TMP}/finished2.txt');
        
        
        fid = fopen(['sub_BNG_scripts_fin.sh'],'w');
        fprintf(fid,'%s \n',temp{:});
        fclose(fid);
        
        pause(2)
        
        unix([temp_save_location '/sub_temp_fin.sh']);
        
    end
    pause(2)
    
    unix([temp_save_location '/sub_temp.sh']);
    
    
    pause(2)
    [~,Is_done] = unix(['wc -l < finished.txt']);
    tempa = str2num(Is_done);
    count = 0;
    while str2num(Is_done) < iters
        pause(10)
        [~,Is_done] = unix(['wc -l < finished.txt']);
        tempb = str2num(Is_done);
        disp(str2num(Is_done));
        if str2num(Is_done)==0 || tempa == tempb
            count = count+1;
            if count > 500
                unix([ ' > ' temp_save_location '/finished.txt']);
                unix(['rm -rf ' temp_save_location '/tempout_*.txt']);

                unix([temp_save_location '/sub_temp.sh']);
                 count = 0;
            end
        end
        tempa = tempb;
        
    end
    
    
    unix(['rm -rf ' temp_save_location '/finished.txt']);
    if  rem(length(replicas_forwards(:,1)),bash_clump) ~=0
        
        while ~exist([temp_save_location '/finished2.txt']);
            pause(10)
            [~,Is_done] = unix(['wc -l < curr_out_rep_' num2str(ijk) '.txt']);
            disp(str2num(Is_done));
        end
        unix(['rm -rf ' temp_save_location '/finished2.txt']);
        
    end
    
    unix(['cat ' temp_save_location '/tempout_*.txt >> ' temp_save_location '/curr_out_rep_' num2str(ijk) '.txt']);
    
    replicas_new = load(['curr_out_rep_' num2str(ijk) '.txt']);
    unix(['rm -rf ' temp_save_location '/finished.txt']);
    
    
    [~,sortind,~] = unique(replicas_new(:,end-1));
    replicas_new = [replicas_new(sortind,2:end-1),weights];
    curr_iteration = iteration;
    iteration = iteration+1;
    if iteration > 3
        iteration = 0;
        
        
    end
    unix('rm -rf init.o*');
    unix('rm -rf tempAtoB.o*');
    unix('rm -rf temp_fin.o*');

    bins = knnsearch(VoronoiLocs,replicas_new(:,1:species));

    b=find(MAP(bins)==2);
    a=find(MAP(bins)==3);
    
    
    temp_switch = zeros(length(replicas_new(:,1)),1);
    temp_switch(a') = 1;
    temp_switch(b') = -1;
    switchA = find(temp_switch-binsr==2);
    switchB = find(temp_switch-binsr==-2);
    
    
    replicas_new(a',end-2) = 1;
    replicas_new(b',end-2) = -1;
    
    
    MFPT(1).weights(ijk,1) = sum(replicas_new(replicas_new(:,end-2)==1,end));
    MFPT(2).weights(ijk,1) = sum(replicas_new(replicas_new(:,end-2)==-1,end));
    MFPT(1).replica_number(ijk,1) = numel(find(replicas_new(:,end-2)==1));
    MFPT(2).replica_number(ijk,1) = numel(find(replicas_new(:,end-2)==-1));
    
    
    
    
    
    [sa,~] = size(switchA);
    [sb,~] = size(switchB);
    
    prob_rep_A = replicas_new(switchB,end);
    MFPT(1).counts(ijk,1) = sb;
    MFPT(1).prob(ijk,1) = sum(prob_rep_A);
    
    
    prob_rep_B = replicas_new(switchA,end);
    MFPT(2).counts(ijk,1) = sa;
    MFPT(2).prob(ijk,1) =sum(prob_rep_B);
    spec_indA = find(replicas_new(:,end-2)== 1);
    spec_indB = find(replicas_new(:,end-2)== -1);
    spec_indC = find(replicas_new(:,end-2)==0);




    BinsprevA = knnsearch(VoronoiLocs,replicas_new(spec_indA,1:species));
    BinsprevB = knnsearch(VoronoiLocs,replicas_new(spec_indB,1:species));
    BinsprevC = knnsearch(VoronoiLocs,replicas_new(spec_indC,1:species));
    
    if ~isempty(spec_indB)
        replicas_forwardsa = WEstep_051117(replicas_new(spec_indA,:),VoronoiLocs,ReplicasRequired,BinsprevA );
        replicas_forwardsb = WEstep_051117(replicas_new(spec_indB,:),VoronoiLocs,ReplicasRequired,BinsprevB );
        replicas_forwardsc = WEstep_051117(replicas_new(spec_indC,:),VoronoiLocs,ReplicasRequired,BinsprevC );
        
        replicas_forwards = [replicas_forwardsa; replicas_forwardsb;replicas_forwardsc];
    else
        replicas_forwards = [WEstep_051117(replicas_new(spec_indA,:),VoronoiLocs,ReplicasRequired,BinsprevA );WEstep_051117(replicas_new(spec_indC,:),VoronoiLocs,ReplicasRequired,BinsprevC )];
    end
    
    
    
    
    
    
    
    matobj_save_info.replicas_forwards = replicas_forwards;
    matobj_save_info.VoronoiLocs = VoronoiLocs;


    matobj_save_info.MFPT = MFPT ;
    matobj_save_info.iteration = iteration;
    
    if rem(ijk,2)==0
        save([save_location ],'VoronoiLocs','replicas_forwards','MFPT','iteration','-v7.3')
        
    end
    
    for i = 1:2
        currentMFPTprobraw(ijk,i) = 1/(MFPT(i).prob(ijk,1)./(MFPT(i).weights(ijk,1))./(tau));
        
    end
    
    
    
    display([currentMFPTprobraw(ijk,1:2)]);    %
    display(ijk);
    display([MFPT(1).weights(ijk,1),MFPT(2).weights(ijk,1)])
    
    if ijk == 1
        
        dlmwrite([temp_dir_transition 'MFPT_prob.txt'],currentMFPTprobraw(ijk,:) );
        
    else
        dlmwrite([temp_dir_transition 'MFPT_prob.txt'],currentMFPTprobraw(ijk,:) ,'-append');
    end
    
    
    
    
    
    
    
    toc = cputime-tic;
    display(toc);
    
    
    
    
    
    
    
end


toc = cputime-tic;
display(toc);
save([save_location ],'VoronoiLocs','replicas_forwards','MFPT','iteration','-v7.3')


end


function [newreps] = WEstep_051117(replicas,Voronoi_List,numreps,bins)
newreps = [];
%This code assumes the columns of the replicas matrix are the species, and that the weights are in the last column
for i = 1:length(Voronoi_List(:,1))
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





