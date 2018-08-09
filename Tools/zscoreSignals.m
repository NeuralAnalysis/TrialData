%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function trial_data = softNormalize(trial_data,params)
% 
% This will zscore signals. Simple. Why can't life always be this easy?
%
% HINT: Can just pass SIGNALS input instead of params struct
%
% INPUTS:
%   trial_data : the struct
%   params     : parameter struct
%      .signals    : which fields (string or cell array of strings)
%                       Defaults to all _spikes fields
%
% OUTPUTS:
%   trial_data : the struct with all signals fields zscored
% 
% Written by Matt Perich. Updated July 2018.
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function trial_data = zscoreSignals(trial_data,params)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PARAMETER DEFAULTS:
signals  =  getTDfields(trial_data,'spikes');
if nargin > 1
    if ~isstruct(params) % must be signals input
        signals = params;
    else % overwrite parameters
        assignParams(who,params);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
signals = check_signals(trial_data,signals);
signals = signals(:,1); % you don't need the idx if they exist, just do it for them all

for i = 1:size(signals,1)
    % compute normalization factors
    for trial = 1:length(trial_data)
        trial_data(trial).(signals{i}) = zscore(trial_data(trial).(signals{i}));
    end
end