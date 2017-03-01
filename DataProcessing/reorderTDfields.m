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
fn_neural = [getTDfields(trial_data,'neural'); getTDfields(trial_data,'unit_guides')];
fn_cont   = getTDfields(trial_data,'cont');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get indices for everything

% the easy ones
meta_idx = find(~(ismember(fn,fn_idx) | ismember(fn,fn_neural) | ismember(fn,fn_cont)));
neural_idx = find(ismember(fn,fn_neural));

% sort idx indices by increasing bin count
idx_idx = find(ismember(fn,fn_idx));
% DO SORT
v = zeros(1,length(idx_idx));
for i = 1:length(idx_idx)
    temp = trial_data(1).(fn{idx_idx(i)});
    if isempty(temp) || any(isnan(temp))
        v(i) = Inf;
    else
        v(i) = temp(1); % in case there are multiple
    end
end
[~,I] = sort(v);
idx_idx = idx_idx(I);


% put EMG at end of continuous
%   Because, aesthetically, I like having emg_names after the rest
cont_idx = find(ismember(fn,fn_cont));
if isfield(trial_data,'emg')
    emg_idx = cellfun(@(x) ~isempty(x),strfind(fn_cont,'emg'));
    cont_idx = [cont_idx; cont_idx(emg_idx)];
    cont_idx(emg_idx) = [];
end

% do reordering
master_idx = [meta_idx; idx_idx; cont_idx; neural_idx];

% order them
trial_data = orderfields(trial_data,master_idx);
