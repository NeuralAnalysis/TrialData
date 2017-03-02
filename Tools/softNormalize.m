%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function trial_data = softNormalize(trial_data,params)
% 
% This will soft normalize a la Churchland papers
%   normalization factor = firing rate range + alpha (Default = 5)
%
% INPUTS:
%   trial_data : the struct
%   params     : parameter struct
%      .signals    : which fields (string or cell array of strings)
%                       Defaults to all _spikes fields
%      .alpha      : normalization factor (default 5)
%
% OUTPUTS:
%   trial_data : the struct with all signals fields normalized
% 
% Written by Matt Perich. Updated Feb 2017.
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function trial_data = softNormalize(trial_data,params)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PARAMETER DEFAULTS:
signals  =  getTDfields(trial_data,'spikes');
alpha    =  5;
if nargin > 1, assignParams(who,params); end % overwrite parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:length(signals)
    % compute normalization factors
    normfac = range(cat(1,trial_data.(signals{i}))) + alpha;
    for trial = 1:length(trial_data)
        trial_data(trial).(signals{i}) = trial_data(trial).(signals{i})./normfac;
    end
end