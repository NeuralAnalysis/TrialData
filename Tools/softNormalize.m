%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function trial_data = softNormalize(trial_data,params)
% 
% This will soft normalize a la Churchland papers
%   normalization factor = firing rate range + alpha (Default = 5)
%
% HINT: instead of params struct, can just pass SIGNALS input if you are
% okay just using the defaults
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
if nargin > 1
    if ~isstruct(params) % must be signals as input
        signals = params;
    else % overwrite parameters
        assignParams(who,params);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isstruct(trial_data), error('First input must be trial_data struct!'); end

signals = check_signals(trial_data,signals);

for i = 1:size(signals,1)
    % compute normalization factors
    normfac = range(get_vars(trial_data,signals)) + alpha;
    for trial = 1:length(trial_data)
        trial_data(trial).(signals{i,1}) = get_vars(trial_data(trial),signals)./normfac;
    end
end
