%% Load data
clear;
clc;

filename = '/Users/mattperich/Data/TrialDataFiles/Chewie_CO_FF_2016-10-07.mat';
% Load data
load(filename);
[~,td] = getTDidx(trial_data,'epoch','BL');
td = truncateAndBin(td,1,{'idx_movement_on',-20},{'idx_movement_on',50});

%% Process
td = getPCA(td,struct('array','M1','do_smoothing',true,'sqrt_transform',true,'trial_avg',true,'trial_avg_cond','target_direction'));


%% Visualize

params = struct( ...
    'trials',       1, ...
    'plot_signals', {{'vel','force'}}, ...
    'plot_pca',false, ...
    'pca_array','M1',...
    'pos_offset',[0 -30]);
visData(td,params)

