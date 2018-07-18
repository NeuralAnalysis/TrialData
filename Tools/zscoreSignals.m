%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function trial_data = softNormalize(trial_data,params)
% 
% This will zscore signals. Simple. Why can't life always be this easy?
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
if nargin > 1, assignParams(who,params); end % overwrite parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isstruct(trial_data), error('First input must be trial_data struct!'); end

for i = 1:length(signals)
    % compute normalization factors
    for trial = 1:length(trial_data)
        trial_data(trial).(signals{i}) = zscore(trial_data(trial).(signals{i}));
    end
end