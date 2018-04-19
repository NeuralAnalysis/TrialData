%% Trial Data Generation
load('Butter_RW_20180405_1_CDS_sorted.mat')
td = parseFileByTrial(cds);
%% TD manipulation
td = binTD(td, 5);
params.go_cue_name ='idx_goCueTime';
params.end_name = 'idx_endTime';
td = getRWMovements(td, params);
td = removeBadTrials(td);
td = trimTD(td, {'idx_movement_on',0}, {'idx_endTime', -1});
%% 