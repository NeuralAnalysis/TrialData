function [trial_data,bad_trials] = pruneBadTrials(trial_data,params)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This will identify and prune bad trials
%
% INPUTS:
%   trial_data : trial data struct
%   params     : (struct) has parameter values
%       .min_trial_time : minimum time from start to end (in # bins)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin > 1
    if isfield(params,'min_trial_time'), min_trial_time = params.min_trial_time; else, min_trial_time = 10; end
end

bad_idx = true(1,length(trial_data));
for iTrial = 1:length(trial_data)
    err = false;
    
    % loop along all indices and make sure they aren't NaN
    fn = fn(cellfun(@(x) strcmpi(x(1:3),'idx'),fieldnames(td)));
    if any(cellfun(@(x) isnan(td.(x)),fn))
        err = true;
    end
    
    td = trial_data(iTrial);
    % there are a few problems I've run into that make things break
    if td.idx_trial_end - td.idx_trial_start < min_trial_time
        err = true;
    end
    
    if err, bad_idx(iTrial) = true; end
end
disp(['Pruning ' num2str(sum(bad_idx)) ' trials with bad trial info...']);
bad_trials = trial_data(bad_idx);
trial_data(bad_idx) = [];