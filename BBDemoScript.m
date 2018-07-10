%% TD Introduction
clear all
load('Han_20171206_TD.mat')
trial_data

%% Trial Data Manipulation
%% getTDFields: returns column headers, with options for which to pick
trial_data_events = getTDfields(trial_data, 'idx');

%% getTDidx: returns indices and trial data structure with the appropriate
%characteristics
[~, trial_data] = getTDidx(trial_data, 'result', 'r', 'epoch', 'BL');
trial_data = trial_data(~isnan([trial_data.target_direction]));

%% binTD: Rebin the data at a new sampling rate
%   inputs: trial_data, bins_to_combine
bins_to_combine = 5;
trial_data = binTD(trial_data,5);

%% removeBadTrials: removes any trials that have NA for any idx field
trial_data_noBumps = removeBadTrials(trial_data);

%% removeBadNeurons: removes neurons that are low firing
trial_data = removeBadNeurons(trial_data);

%% getMoveOnsetAndPeak: adds new idx heading for movement-onset and peak
% speed
trial_data = getMoveOnsetAndPeak(trial_data, params);


%% getSpeed: adds new 'speed' column

trial_data = getSpeed(trial_data);


%% trimTD: cut the data from your trial_data into sections
trial_data_movement = trimTD(trial_data, {'idx_movement_on', 0}, {'idx_movement_on', 6});
trial_data_premovement = trimTD(trial_data, {'idx_movement_on', -5}, {'idx_movement_on', -2});

%% getVars: get the data out of all trials of trial_data structure
velocity = get_vars(trial_data_movement, {'vel',1:2});
% alternatively:
velocity2 = cat(1, trial_data_movement.vel);
speed1 = cat(2, trial_data_movement.speed);


%% Trial data analysis functions:
%% getModel: fits a model between chosen inputs. This model is general, and
% can be neural nets, glm, linear regression. Can be expanded to fit new
% models

params.model_type    =  'nn';
params.model_name    =  'movementNN';
params.in_signals    =  {'S1_spikes'};% {'name',idx; 'name',idx};
params.out_signals   =  {'vel'};% {'name',idx};
params.train_idx     =  1:length(trial_data_movement);
% GLM-specific parameters
params.do_lasso      =  false;
params.lasso_lambda  =  0;
params.lasso_alpha   =  0;
params.nn_params = [10,20];

params.glm_distribution     =  'normal';
[trial_data_glm, model_info] = getModel(trial_data_movement, params);
%%
params.eval_metric      =  'vaf';
params.num_boots        =  1000;
r2_glm = evalModel(trial_data_glm, params);

figure
plot(cat(1,trial_data_glm.vel))
yyaxis right
plot(cat(1,trial_data_glm.nn_movementNN))
xlim([1000,2000])

%% getPCA: perform PCA on a signal in your trial_data struct
paramPCA = struct('signals', 'S1_spikes');
trial_data_pca = getPCA(trial_data_movement, paramPCA); 
trial_data_pca = smoothSignals(trial_data_pca, struct('signals', 'S1_pca'));
trial_data = trialAverage(trial_data_pca, 'target_direction');

dirs = unique([trial_data.target_direction]);

[~,rightMove] = getTDidx(trial_data_pca, 'target_direction', 0);
[~,upMove] = getTDidx(trial_data_pca, 'target_direction', dirs(5));
[~,leftMove] = getTDidx(trial_data_pca, 'target_direction',  dirs(9));
[~,downMove] = getTDidx(trial_data_pca, 'target_direction', dirs(13));

s1PCARight = cat(1, rightMove.S1_pca);
s1PCAUp = cat(1, upMove.S1_pca);
s1PCALeft =cat(1,leftMove.S1_pca);
s1PCADown = cat(1,downMove.S1_pca);
colors = linspecer(4);

figure
plot3(s1PCARight(:,1), s1PCARight(:,2), s1PCARight(:,3), 'Color', colors(1,:))
hold on
plot3(s1PCAUp(:,1), s1PCAUp(:,2), s1PCAUp(:,3),'Color', colors(2,:))
plot3(s1PCALeft(:,1), s1PCALeft(:,2), s1PCALeft(:,3),'Color', colors(3,:))
plot3(s1PCADown(:,1), s1PCADown(:,2), s1PCADown(:,3),'Color', colors(4,:))

axis equal

%% getTDPDs: Get bootstrapped PDs for all of your neurons:
trial_data_pca = removeBadNeurons(trial_data_pca);
paramsPDs.out_signals      =  ['S1_spikes'];
paramsPDs.out_signal_names = trial_data_pca(1).S1_unit_guide;
pdTable= getTDPDs(trial_data_pca, paramsPDs);


%% How to get your hands on a TD
load('Butter_RW_20180405_1_CDS_sorted.mat')
trial_data_butter= parseFileByTrial(cds);

%% Other features
%% Timestamp and real time alignment
paramsCDS.include_ts = true;
paramsCDS.include_start = true;
trial_data_ts = parseFileByTrial(cds, paramsCDS);
% Notice this takes MUCH MUCH LONGER and much more space
% only use these features if you need them
%% Feature requests?