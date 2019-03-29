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
% Written by Matt Perich. Updated March 2019.
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function trial_data = softNormalize(trial_data,params)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PARAMETER DEFAULTS:
signals  =  getTDfields(trial_data,'spikes');
alpha    =  5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Some undocumented extra parameters
field_extra  =  '';   % if empty, defaults to input field name(s)
verbose      =  false;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin > 1
    if ~isstruct(params) % must be signals as input
        signals = params;
    else % overwrite parameters
        assignParams(who,params);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
trial_data = check_td_quality(trial_data);
signals = check_signals(trial_data,signals);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% check output field addition
field_extra  = check_field_extra(field_extra,signals);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
for iSig = 1:size(signals,1)
    % compute normalization factors
    normfac = range(get_vars(trial_data,signals(iSig,:))) + alpha;
    for trial = 1:length(trial_data)
        trial_data(trial).([signals{iSig,1} field_extra{iSig}]) = get_vars(trial_data(trial),signals(iSig,:))./normfac;
    end
end
