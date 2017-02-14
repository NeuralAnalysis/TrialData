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
function trial_data = getSpeed(trial_data)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for trial = 1:length(trial_data)
    trial_data(trial).speed = sqrt(trial_data(trial).vel(:,1).^2 + ...
        trial_data(trial).vel(:,2).^2);
end

% restore logical order
trial_data = reorderTDfields(trial_data);