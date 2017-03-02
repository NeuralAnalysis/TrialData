%% Load data
clear;
clc;

filename = '/Users/mattperich/Data/TrialDataFiles/Chewie_CO_FF_2016-10-07.mat';
% Load data
load(filename);

[~,td] = getTDidx(trial_data,'epoch','BL','result','R');
td = truncateAndBin(td,{'idx_target_on',0},{'idx_trial_end',-20});

td = removeBadTrials(td);
td = removeBadNeurons(td,struct('min_fr',3));

%% Visualize

params = struct( ...
    'trials',       1, ...
    'plot_signals', {{'vel'}}, ...
    'plot_pca',false, ...
    'pca_array','M1',...
    'pos_offset',[0 -30]);
visData(td,params)

