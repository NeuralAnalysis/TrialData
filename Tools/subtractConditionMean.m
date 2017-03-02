%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% trial_data = subtractConditionMean(trial_data, params)
%
%   Finds the mean across all conditions for each time point and subtracts
% this value. Works on all time-varying signals. This is what Churchland
% and Co. do for their jPCA work etc. Note that this is designed to be done
% AFTER a call to trialAverage, and thus assumes that each entry represents
% a condition that you want to average over.
%
% INPUTS:
%   trial_data : the struct
%   params     : parameter struct
%     .cond_idx : which condition entries to use for mean
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
record_flag = true; % will add a flag field
if nargin > 1, assignParams(who,params); end % overwrite parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fn = getTDfields(trial_data,'time');

td = trial_data(cond_idx);

for i = 1:length(fn)
    m = mean(cat(3,td.(fn{i})),3);
    for trial = 1:length(trial_data)
        trial_data(trial).(fn{i}) = trial_data(trial).(fn{i}) - m;
    end
end

% let it be known that on this day a mean was subtracted!
if record_flag
    for trial = 1:length(trial_data)
        trial_data(trial).cond_mean_subtracted = true;
    end
end
