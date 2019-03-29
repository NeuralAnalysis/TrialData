%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function trial_data = sqrtTransform(trial_data,signals)
%
% This function will smooth and/or square root transform spikes
%
% INPUTS:
%   trial_data : the struct
%   params     : params struct
%       .signals  :  which signals to work on. Defaults to all _spikes
%
% OUTPUTS:
%   trial_data : the struct with all (signals{}) fields transformed
%
% Written by Matt Perich. Updated March 2019.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function trial_data = sqrtTransform(trial_data,params)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
signals  =  '';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Some undocumented extra parameters
field_extra  =  '';   % if empty, defaults to input field name(s)
verbose      =  false;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin > 1
    if isstruct(params)
        assignParams(who,params); % overwrite parameters
    else % default to signals as input
        signals = params;
    end
else
    disp('No signals provided, square root transforming spikes');
    signals = getTDfields(trial_data,'spikes');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
trial_data = check_td_quality(trial_data);
% make sure the signal input formatting is good
signals = check_signals(trial_data(1),signals);
signals = signals(:,1); % don't need the idx if they exist
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% check output field addition
field_extra  = check_field_extra(field_extra,signals);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
for trial = 1:length(trial_data)
    for iSig = 1:size(signals,1)
        data = getSig(trial_data(trial),signals{iSig,:});
        if any(data < 0), warning('Negative values in signals. Sqrt_transform will be imaginary.'); end
        trial_data(trial).([signals{iSig,1} field_extra{iSig}]) = sqrt(data);
    end
end
