% load data
clear;
close all;
clc;

% put your filename here to load an existing one
%   or call trial_data = parseFileByTrial(cds);
load('/Users/mattperich/Data/TrialDataFiles/Chewie_CO_FF_2016-10-07.mat');

%%
n_bins = 3;
spiking_inputs = {'M1_spikes'};
% note: to do M1 and PMd arrays:
%       spiking_inputs = {'M1_spikes','PMd_spikes'};

% this is a curl field experiment, but I just want baseline reward trials
[~,trial_data] = getTDidx(trial_data,'epoch','BL','result','R');
% I want to smooth the data to get ready for PCA
trial_data = smoothSignals(trial_data,struct( ...
    'signals',{spiking_inputs}, ...
    'sqrt_transform',true, ...
    'do_smoothing', true, ...
    'kernel_SD', 0.1));

% get struct that has preparatory and movement activity
td = truncateAndBin(trial_data,n_bins,{'idx_movement_on',-50},{'idx_movement_on',70});
% get struct that has only preparatory activity
td_prep = truncateAndBin(trial_data,n_bins,{'idx_go_cue',-30},{'idx_go_cue',0});
% get struct that has only movement activity
td_move = truncateAndBin(trial_data,n_bins,{'idx_movement_on',10},{'idx_movement_on',40});

% Now all of the trials will be uniform lenght for each data section, so
% it's easy to trial average, which will clean up our results
td = trialAverage(td,'target_direction');
td_prep = trialAverage(td_prep,'target_direction');
td_move = trialAverage(td_move,'target_direction');

%%
% get subspaces
arrays = strrep(spiking_inputs,'_spikes','');
pca_params = struct( ...
    'signals',{spiking_inputs});
% get covariance matrices and means
[~,pca_prep] = getPCA(td_prep,pca_params);
[~,pca_move] = getPCA(td_move,pca_params);

%% Project activity for whole trial into the different subspaces
td_td_prep   = getPCA(td,setfield(pca_params,'w',pca_prep.w));
td_td_move   = getPCA(td,setfield(pca_params,'w',pca_move.w));

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



