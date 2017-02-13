% load data
clear;
close all;
clc;

load('/Users/mattperich/Data/TrialDataFiles/Chewie_CO_FF_2016-10-07.mat');

%%
n_bins = 3;

[~,trial_data] = getTDidx(trial_data,'epoch','BL');
% get preparatory activity
td = truncateAndBin(trial_data,n_bins,{'idx_movement_on',-50},{'idx_movement_on',70});
td_prep = truncateAndBin(trial_data,n_bins,{'idx_go_cue',-30},{'idx_go_cue',0});
td_move = truncateAndBin(trial_data,n_bins,{'idx_movement_on',10},{'idx_movement_on',40});

%%
% get subspaces
arrays = {'M1'};
pca_params = struct('array',{arrays},'do_smoothing',true,'bin_size',0.01*n_bins,'kernel_SD',0.05,'trial_avg',true,'trial_avg_cond','target_direction');
[w,mu] = getPCA(td,pca_params);
[w_prep,mu_prep] = getPCA(td_prep,pca_params);
[w_move,mu_move] = getPCA(td_move,pca_params);

%% Project activity into different subspaces
td_td_prep   = getPCA(td,w_prep,mu,pca_params);
td_td_move   = getPCA(td,w_move,mu,pca_params);

%%
figure;
subplot(2,1,1); hold all;
for i = 1:8
    temp = td_td_prep(i).([[arrays{:}] '_pca']);
    plot(temp(:,1),'LineWidth',2);
end
axis('tight');
set(gca,'Box','off','TickDir','out','FontSize',14,'YLim',[-1.5 1.5]);
ylabel('First PC Proj','FontSize',14)
title('Preparatory Subspace','FontSize',14);
subplot(2,1,2); hold all;
for i = 1:8
    temp = td_td_move(i).([[arrays{:}] '_pca']);
    plot(temp(:,1),'LineWidth',2)
end
axis('tight');
set(gca,'Box','off','TickDir','out','FontSize',14,'YLim',[-1.5 1.5]);
ylabel('First PC Proj','FontSize',14);
title('Movement Subspace','FontSize',14);



