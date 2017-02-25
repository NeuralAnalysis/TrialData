%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function trial_data = sqrtTransform(trial_data,signals)
%
% This function will smooth and/or square root transform spikes
%
% INPUTS:
%   trial_data : the struct
%   signals    : which signals to work on
%
% OUTPUTS:
%   trial_data : the struct with all (signals{}) fields transformed
%
% Written by Matt Perich. Updated Feb 2017.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function trial_data = sqrtTransform(trial_data,signals)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 2, error('Must provide one or more signals to transform.'); end
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
