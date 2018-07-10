%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function trial_data = getSpeed(trial_data)
%
%   Adds speed as an entry to each trial.
%
% INPUTS:
%   trial_data : the struct
%
% OUTPUTS:
%   trial_data : the struct with .speed fields (reordered logically)
%
% Written by Matt Perich. Updated Feb 2017.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function trial_data = getSpeed(trial_data,params)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
velocity_variable = 'vel';
if nargin > 1, assignParams(who,params); end % overwrite defaults
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isstruct(trial_data), error('First input must be trial_data struct!'); end

for trial = 1:length(trial_data)
    temp = trial_data(trial).(velocity_variable);
    trial_data(trial).speed = sqrt(temp(:,1).^2 + ...
        temp(:,2).^2);
end

% restore logical order
trial_data = reorderTDfields(trial_data);