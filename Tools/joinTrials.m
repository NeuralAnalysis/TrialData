%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function td_j = joinTrials(trial_data,params)
%
%   Joins trials into one mega trial. May include future functionality to join
% by key meta fields, but this just joins meta fields into a cell array if not
% all the same
%
% INPUTS:
%   trial_data : the struct
%   params     : parameter struct
%       .key_fields    : (string) name of meta fields to join on
%                               (default: '' - join all trials)
%
% OUTPUTS:
%   td_j : the struct joined together
%
% Written by Raeed Chowdhury. Updated November 2018.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function td_j = joinTrials(trial_data,params)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PARAMETER DEFAULTS
% key_fields = '';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Some undocumented extra parameters
verbose = false;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin > 1
    assignParams(who,params);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
trial_data = check_td_quality(trial_data);
% if ~iscell(key_fields), key_fields = {key_fields}; end

td_j = struct();

fn_time = getTDfields(trial_data,'time');
fn_labels = getTDfields(trial_data,'labels');
fn_meta = getTDfields(trial_data,'meta');
fn_idx = getTDfields(trial_data,'idx');

% add signals
for fn_ctr = 1:length(fn_time)
    td_j.(fn_time{fn_ctr}) = getSig(trial_data,fn_time{fn_ctr});
end

% add signal labels
for fn_ctr = 1:length(fn_labels)
    % figure out if labels are all the same
    if isequal(trial_data.(fn_labels{fn_ctr}))
        labels = trial_data(1).(fn_labels{fn_ctr});
    else
        labels = {trial_data.(fn_labels{fn_ctr})};
    end
    td_j.(fn_labels{fn_ctr}) = labels;
end

% add meta stuff
for fn_ctr = 1:length(fn_meta)
    % figure out if meta are all the same
    if isequal(trial_data.(fn_meta{fn_ctr}))
        meta = trial_data(1).(fn_meta{fn_ctr});
    else
        meta = {trial_data.(fn_meta{fn_ctr})};
    end
    td_j.(fn_meta{fn_ctr}) = meta;
end

% add idx stuff
for fn_ctr = 1:length(fn_idx)
    idx = [];
    last_end = 0;
    for trialnum = 1:length(trial_data)
        new_idx = trial_data(trialnum).(fn_idx{fn_ctr});
        idx = [idx new_idx+last_end];

        timevec = trial_data(trialnum).(fn_time{1});
        last_end = last_end + size(timevec,1);
    end

    td_j.(fn_idx{fn_ctr}) = idx;
end

td_j = reorderTDfields(td_j);
