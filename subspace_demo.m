% load data
clear;
close all;
clc;

load('/Users/mattperich/Data/TrialDataFiles/Chewie_CO_VR_2016-10-06.mat');

%%
n_bins = 1;

[~,trial_data] = getTDidx(trial_data,'epoch','BL');
% get preparatory activity
td = truncateAndBin(trial_data,n_bins,{'idx_movement_on',-50},{'idx_movement_on',70});
td_prep = truncateAndBin(trial_data,n_bins,{'idx_go_cue',-30},{'idx_go_cue',0});
td_move = truncateAndBin(trial_data,n_bins,{'idx_movement_on',10},{'idx_movement_on',40});

%%
% get subspaces
arrays = {'M1','PMd'};
pca_params = struct('array',{arrays},'do_smoothing',true,'bin_size',0.01*n_bins,'kernel_SD',0.05,'trial_avg',true,'trial_avg_cond','target_direction');
[w,mu] = getPCA(td,pca_params);
[w_prep,mu_prep] = getPCA(td_prep,pca_params);
[w_move,mu_move] = getPCA(td_move,pca_params);

%%
figure;
subplot(121); imagesc(w_prep,[-1,1]); axis('square');
subplot(122); imagesc(w_move,[-1,1]); axis('square');

%% Project activity into different subspaces
td_td_prep   = getPCA(td,w_prep,mu,pca_params);
td_td_move   = getPCA(td,w_move,mu,pca_params);
td_move_prep = getPCA(td_move,w_prep,mu_move,pca_params);
td_move_move = getPCA(td_move,w_move,mu_move,pca_params);
td_prep_move = getPCA(td_prep,w_move,mu_prep,pca_params);
td_prep_prep = getPCA(td_prep,w_prep,mu_prep,pca_params);

%% plot all movement activity in preparatory subspace
trial_idx = 1:length(td_move_prep);
plot_dims = 1:3;
figure;
temp = cat(1,td_move_prep(trial_idx).([[arrays{:}] '_pca']));
subplot(221); plot(temp(:,plot_dims)); set(gca,'YLim',[-1 1]);
temp = cat(1,td_move_move(trial_idx).([[arrays{:}] '_pca']));
subplot(222); plot(temp(:,plot_dims)); set(gca,'YLim',[-1 1]);

temp = cat(1,td_prep_prep(trial_idx).([[arrays{:}] '_pca']));
subplot(223); plot(temp(:,plot_dims)); set(gca,'YLim',[-1 1]);
temp = cat(1,td_prep_move(trial_idx).([[arrays{:}] '_pca']));
subplot(224); plot(temp(:,plot_dims)); set(gca,'YLim',[-1 1]);

%%
figure;
subplot(2,1,1); hold all;
for i = 1:8
    temp = td_td_prep(i).([[arrays{:}] '_pca']);
    plot(temp(:,1),'LineWidth',2);
end
axis('tight');
set(gca,'Box','off','TickDir','out','FontSize',14,'YLim',[-1 1]);
ylabel('First PC Proj','FontSize',14)
title('Preparatory Subspace','FontSize',14);
subplot(2,1,2); hold all;
for i = 1:8
    temp = td_td_move(i).([[arrays{:}] '_pca']);
    plot(temp(:,1),'LineWidth',2)
end
axis('tight');
set(gca,'Box','off','TickDir','out','FontSize',14,'YLim',[-1 1]);
ylabel('First PC Proj','FontSize',14);
title('Movement Subspace','FontSize',14);



