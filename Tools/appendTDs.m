%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function combined = appendTDs(varargin)
%
%   Stitches together an arbitrary number of trial data structs in time.
% This is useful if you want to plot discontinuous data, for instance 500
% ms after target presentation and 500 ms after movement onset in a
% variable delay center out task.
%
% To do:
%   1) check that all have same metadata
%   2) deal with idx_?
%
% INPUTS:
%   varargin : any number of trial_data structs
%       Note: must have same trials in the same order
%
% OUTPUTS:
%   combined : a new trial_data struct with each trial combined
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function combined = appendTDs(varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

combined = varargin{1};
if length(varargin) > 1
    % check that all have same number of trials
    if length(unique(cellfun(@(x) length(x),varargin))) > 1
        error('All inputs must have same number of trials');
    end
    % get the list of time variables
    fn_time = getTDfields(combined,'time');
    
    for i = 2:length(varargin)
        td = varargin{i};
        for trial = 1:length(combined)
            for j = 1:length(fn_time)
                combined(trial).(fn_time{j}) = [combined(trial).(fn_time{j}); td(trial).(fn_time{j})];
            end
        end
    end
else
    warning('Only one trial_data struct provided. Nothing to append.');
end
