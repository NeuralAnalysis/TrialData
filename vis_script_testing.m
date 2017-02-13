%% Load data
clear;
clc;


% Pat's computer
% data_folder = 'C:\Users\pnlawlor\Box Sync\PatAndMattData';

% Matt's computer
data_folder = '/Users/mattperich/Data/TrialDataFiles';

data_fname = 'Chewie_CO_FF_2016-10-07.mat';
% Load data
load([data_folder '/' data_fname])

%% Load scripts

% addpath(genpath('C:\Users\pnlawlor\GoogleDrive\Research\Projects\PatAndMatt\Scripts_PatAndMatt'))

%% Choose visualization parameters

param_struct.trials = 10;
param_struct.plot_signals = {'pos','vel','force'};
param_struct.pos_location = 'left';


%% Visualize

vis_data(trial_data,param_struct)

