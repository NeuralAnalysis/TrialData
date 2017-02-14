%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function trial_data = reorderTDfields(trial_data)
%
% This is my total OCD-appeasing function. It will:
%   1) Put metadata at the top
%   2) Put idx_ fields second, in order of value for each trial
%   3) put continuous/time signals third (EMG last)
%   4) put spikes and neural fields last (and group by array name)
%
% Written by Matt Perich. Updated Feb 2017.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function trial_data = reorderTDfields(trial_data)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get names for the different types
fn        = fieldnames(trial_data);
fn_idx    = getTDfields(trial_data,'idx');
fn_neural = getTDfields(trial_data,'neural');
fn_cont   = getTDfields(trial_data,'cont');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get indices for everything

% the easy ones
meta_idx = find(~(ismember(fn,fn_idx) | ismember(fn,fn_neural) | ismember(fn,fn_cont)));
neural_idx = find(ismember(fn,fn_neural));

% sort idx indices by increasing bin count
idx_idx = find(ismember(fn,fn_idx));
% DO SORT

% put EMG at end of continuous
%   Because, aesthetically, I like having emg_names after the rest
cont_idx = find(ismember(fn,fn_cont));

% do reordering
master_idx = [meta_idx; idx_idx; cont_idx; neural_idx];

% order them
trial_data = orderfields(trial_data,master_idx);
