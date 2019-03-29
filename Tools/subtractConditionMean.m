%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% trial_data = subtractConditionMean(trial_data, params)
%
%   Finds the mean across all conditions for each time point and subtracts
% this value. Works on all time-varying signals. This is what Churchland
% and Co. do for their jPCA work etc.
%
% INPUTS:
%   trial_data : the struct
%   params     : parameter struct
%     .cond_idx : which trials to use for mean
%                   Note: will subtract this mean from all entries
%
% OUTPUTS:
%   trial_data : the struct, with mean subtracted for each signal
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function trial_data = subtractConditionMean(trial_data,params)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cond_idx  =  1:length(trial_data);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Some undocumented extra parameters
verbose = false;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin > 1, assignParams(who,params); end % overwrite parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
trial_data = check_td_quality(trial_data);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fn = getTDfields(trial_data,'time');

td = trial_data(cond_idx);

for i = 1:length(fn)
    m = mean(cat(3,td.(fn{i})),3);
    for trial = 1:length(trial_data)
        trial_data(trial).(fn{i}) = trial_data(trial).(fn{i}) - m;
    end
end


