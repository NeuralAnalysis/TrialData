%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [trial_data, bad_trials] = pruneBadTrials(trial_data, params);
%
% This function will identify and remove bad trials. Not much supported at
% the moment, but functionality can be expanded. Will remove any trials
% with any idx_ fields that are NaN.
%
% INPUTS:
%   trial_data : trial data struct
%   params     : (struct) has parameter values
%       .trial_time : [min,max] time from start to end (in # bins)
%
% OUTPUTS:
%   trial_data : struct with bad trials removed
%   bad_trials : struct containing said bad trials
%
% Written by Matt Perich. Updated Feb 2017.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [trial_data,bad_trials] = pruneBadTrials(trial_data,params)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
trial_time = [0,Inf];
if nargin > 1
    eval(structvars(length(fieldnames(params)),params)); %get parameters
end

bad_idx = false(1,length(trial_data));
for iTrial = 1:length(trial_data)
    err = false;
    
    td = trial_data(iTrial);
    
    % loop along all indices and make sure they aren't NaN
    fn = getTDfields(td,'idx');
    if any(cellfun(@(x) isnan(td.(x)),fn))
        err = true;
    end
    
    % there are a few problems I've run into that make things break
    if td.idx_trial_end - td.idx_target_on < trial_time(1) || ...
            td.idx_trial_end - td.idx_target_on > trial_time(2)
        err = true;
    end
    
    if err, bad_idx(iTrial) = true; end
end
disp(['Pruning ' num2str(sum(bad_idx)) ' trials with bad trial info...']);
bad_trials = trial_data(bad_idx);
trial_data(bad_idx) = [];