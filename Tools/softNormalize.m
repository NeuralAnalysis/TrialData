%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function trial_data = softNormalize(trial_data,alpha)
% 
% This will soft normalize a la Churchland papers
%   normalization factor = firing rate range + alpha (Default = 5)
%
% INPUTS:
%   trial_data : the struct
%   alpha      : for normalization factor
%
% OUTPUTS:
%   trial_data : the struct with all _spikes fields normalized
% 
% Written by Matt Perich. Updated Feb 2017.
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function trial_data = softNormalize(trial_data,alpha)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PARAMETER DEFAULTS:
if nargin == 1, alpha = 5; end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fn_spikes = getTDfields(trial_data,'spikes');

for i = 1:length(fn_spikes)
    
    % compute normalization factors
    normfac = range(cat(1,trial_data.(fn_spikes{i}))) + alpha;
    for trial = 1:length(trial_data)
        trial_data(trial).(fn_spikes{i}) = trial_data(trial).(fn_spikes{i})./normfac;
    end
end