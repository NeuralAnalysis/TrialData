clc;
clear;
close all;

load('/Users/mattperich/Data/TrialDataFiles/Chewie_CO_FF_2016-09-15.mat')
% load('/Users/mattperich/Data/TrialDataFiles/old/Chewie_CO_VR_2016-09-14.mat')
% load('/Users/mattperich/Data/han_td.mat','trial_data');

[~,trial_data] = getTDidx(trial_data,'result','R');
for i = 1:length(trial_data),trial_data(i).bin_size = 0.01; end

signals = {'PMd_spikes', 'all'; 'M1_spikes','all'};
% signals = {'S1_spikes'};
td = removeBadNeurons(trial_data,struct('min_fr',1,'do_shunt_check',false));
td = removeBadTrials(td,struct('ranges',{{'idx_go_cue','idx_movement_on',[5 50];'idx_movement_on','idx_trial_end',[50,120]}}));
td = smoothSignals(td,struct('signals',{signals},'sqrt_transform',true,'do_smoothing',true,'kernel_SD',0.1));

% td = softNormalize(td);

% td = trimTD(td,{'idx_bumpTime',0},{'idx_bumpTime',50});
td = trimTD(td,{'idx_movement_on',-20},{'idx_movement_on',50});

[~,td] = getTDidx(td,'epoch',{'BL','AD'});

%%
clear blocks;
blocks{1} = getTDidx(td,'epoch','AD','range',[0 0.33]);
blocks{2} = getTDidx(td,'epoch','AD','range',[0.33 0.66]);
blocks{3} = getTDidx(td,'epoch','AD','range',[0.66 1]);
getDPCA(td,'target_direction',blocks,struct('signals',{signals}));