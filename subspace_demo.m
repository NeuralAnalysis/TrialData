% load data
clear;
close all;
clc;

% put your filename here to load an existing one
%   or call trial_data = parseFileByTrial(cds);
load('/Users/mattperich/Data/TrialDataFiles/Chewie_CO_FF_2016-10-07.mat');

%%
n_bins = 3;

[~,trial_data] = getTDidx(trial_data,'epoch','BL');

% get struct that has preparatory and movement activity
td = truncateAndBin(trial_data,n_bins,{'idx_movement_on',-50},{'idx_movement_on',70});
% get struct that has only preparatory activity
td_prep = truncateAndBin(trial_data,n_bins,{'idx_go_cue',-30},{'idx_go_cue',0});
% get struct that has only movement activity
td_move = truncateAndBin(trial_data,n_bins,{'idx_movement_on',10},{'idx_movement_on',40});

%%
% get subspaces
arrays = {'M1','PMd'}; % uses neurons from both arrays
pca_params = struct( ...
    'array',{arrays}, ...
    'do_smoothing',true, ...
    'bin_size',0.01*n_bins, ...
    'kernel_SD',0.05, ...
    'trial_avg',true, ...
    'trial_avg_cond','target_direction');
% we want the means for each neurons when we project the data later
[~,temp] = getPCA(td,pca_params);
mu = temp.mu; clear temp;
% get covariance matrices and means
[~,temp] = getPCA(td_prep,pca_params);
w_prep = temp.w;
[~,temp] = getPCA(td_move,pca_params);
w_move = temp.w;

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



