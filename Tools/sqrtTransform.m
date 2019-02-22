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
% Written by Matt Perich. Updated Feb 2017.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function trial_data = sqrtTransform(trial_data,params)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Some undocumented extra parameters
signals = '';
verbose = false;
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
if ~isstruct(trial_data), error('First input must be trial_data struct!'); end
if ~iscell(signals), signals = {signals}; end
% make sure the signal input formatting is good
signals = check_signals(trial_data(1),signals);
signals = signals(:,1); % don't need the idx if they exist
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for trial = 1:length(trial_data)
    for i = 1:length(signals)
        data = trial_data(trial).(signals{i});
        if any(data < 0), warning('Negative values in signals. Sqrt_transform will be imaginary.'); end
        trial_data(trial).(signals{i}) = sqrt(data);
    end
end
