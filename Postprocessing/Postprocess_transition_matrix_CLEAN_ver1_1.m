clear all
close all


%%Finds the eigenvalues and eigenvectors of the transition matrix, and
%%formats the data for use in MSMBuilder and PyEMMA.


%%%%%%%%%******************************************************
images=1; %Set this value to 1 if you want a plot of the eigenvalues/
%lifetimes and a scatter plot of the voronoi centeres.


tau = 10;
bins = 300;
%File name, without .mat root, of the output of adpative WE sampling
file_name = 'trans_nanog_tau10_021218_f10_4';
load(file_name);

%%%For older versions, the variable trans_matrix_prob needs to be converted
%%%to trans_matrix;

trans_matrix = trans_matrix_prob;

%creates an output directory with the same name as the file_name
out_dir = [file_name '/'];
mkdir(out_dir);

%%Counts  number of unfilled Voronoi regions.
unfilled = find(trans_matrix<0);

%%sets unfilled to 0 if unfilled is empty, meaning that there are no empty
%%Voronoi regions.
if isempty(unfilled)
    unfilled = 0;
end
disp(['There are ' num2str(unfilled) ' unfilled Voronoi bins.']);

%%Normalizes the transition matrix such that rows sum to 1, and removes
%%any rows (and their associated columns) that are entirely zero.
trans_matrix_temp_var_raw = trans_matrix;
row_summation_trans_matrix = sum(trans_matrix_temp_var_raw,2);
idx_zeros_in_trans_matrix = find(row_summation_trans_matrix==0);
trans_matrix_temp_var_raw(idx_zeros_in_trans_matrix,:) = [];
trans_matrix_temp_var_raw(:,idx_zeros_in_trans_matrix) = [];
row_summation_trans_matrix(idx_zeros_in_trans_matrix) = [];
trans_matrix_tot_row=repmat(row_summation_trans_matrix,1,size(trans_matrix_temp_var_raw,2));
T=trans_matrix_temp_var_raw./trans_matrix_tot_row; %this is the transition probability that you should give to Brian

%%removes bins from the Voronoi centers that are not visited.
StatesList = VoronoiLocs;
StatesList(idx_zeros_in_trans_matrix,:) = [];


%%Finds the eigevnalues, eigenvectors, lifetimes, and probability
%%distribution of the transition matrix T.
[Eigenvectors,Eigenvalues]=eig(transpose(T));
[sorted_Eigevnalues,ik]=sort((real(diag(Eigenvalues))));
sorted_Eigenvalues = flipud(sorted_Eigevnalues(1:end));
lifetimes = -tau./log(flipud(sorted_Eigevnalues(1:end-1)));
Eigenvectors = Eigenvectors(:,flipud(ik));
Probability_Surface = ((Eigenvectors(:,1))./sum(Eigenvectors(:,1)));



save([out_dir 'data.mat'],'T','trans_matrix','Probability_Surface','lifetimes','StatesList','VoronoiLocs','Eigenvectors','sorted_Eigenvalues','-v7.3');

%Changing the naming convention for use in MSMbuilder


save([out_dir 'StatesList.mat'],'StatesList');
TransitionMatrix = T;
save([out_dir 'TransitionMatrix.mat'],'TransitionMatrix');
DDs = sorted_Eigenvalues;
VVs = Eigenvectors;
save([out_dir 'DDs.mat'],'DDs');
save([out_dir 'VVs.mat'],'VVs');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Plots the eigenvalues and lifetimes

fontsize = 20;
linewidth = 2;
markersize = 25;
species1 = 2;
species2 = 3;
eigshow = 20;
prob_scale = 50000;

x1 = 1:1:eigshow;
x2 = 2:1:eigshow;


figure(1)
set(1,'PaperUnits','Inches','PaperSize',[7.5,5],'PaperPosition',[0 0 7.5,5],'Position',[100 100 750,500])
set(1,'Color','w');


plot(x1',sorted_Eigenvalues(1:eigshow),'- . k','LineWidth',linewidth,'MarkerSize',markersize);
set(gca,'FontSize',fontsize);
axis(gca,'square');

set(gca,'Ylim',[0 1])
set(gca,'YTick',[0 0.2 0.4 0.6 0.8 1]);
ylabel('Eigenvalues')
xlabel('Eigenvalue Index')

title(['\tau = ' num2str(tau)]);
eig_curr = sorted_Eigenvalues(1:eigshow);
indicies = [0.2,0.4,0.6,0.8,1];
tscale = lifetimes(1:eigshow-1);

yind = -tau./log(indicies);
ax1 = gca;
set(ax1,'FontSize',fontsize)
axis(ax1,'square');

set(gca,'Xlim',[0 eigshow]);

ax1_pos = ax1.Position;

ax2 = axes( 'Position',ax1_pos,   'YAxisLocation','right',...
    'Color','none');
set(ax2,'FontSize',fontsize)
set(ax2,'Ylim',[0 1])
set(ax2,'YTick',[0 0.2 0.4 0.6 0.8 1]);
set(ax2,'YTickLabel',{0 yind(1) yind(2) yind(3) yind(4) '\infty'})
set(ax2,'Xlim',[0 eigshow]);
ylabel(ax2,'Lifetimes (t)')
axis(ax2,'square');

if images == 1
    print([out_dir 'eigenvalues.png'],'-dpng')
end



%%Scatter plots species1 vs. species2
figure(2)
scatter(StatesList(:,species1),StatesList(:,species2),(Probability_Surface).*prob_scale,(Probability_Surface),'filled')
hold on
voronoi(StatesList(:,species1),StatesList(:,species2),'k');
set(gca,'FontSize',fontsize)
axis('square')

box on


if images == 1
    print([out_dir 'scatter.png'],'-dpng')
end



