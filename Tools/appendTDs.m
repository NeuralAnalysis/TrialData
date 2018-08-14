%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function combined = appendTDs(varargin)
%
%   Stitches together an arbitrary number of trial data structs in time.
% This is useful if you want to plot discontinuous data, for instance 500
% ms after target presentation and 500 ms after movement onset in a
% variable delay center out task.
%
% NOTE: the way this works, there are some best practices. You should make
% every struct you append have no overlapping time values (unless you want
% them to repeat). You should also pass them in chronologically, so the
% idx_ fields continue to have something resembling meaning (unless you
% don't want them to).
%
% To do:
%   1) check that all have same metadata
%
% INPUTS:
%   varargin : any number of trial_data structs
%       Note: must have same trials in the same order
%
% OUTPUTS:
%   combined : a new trial_data struct with each trial combined
%                adds field 'stitch_marks', which has the bin numbers
%                corresponding to the beginning of each discontinuity
%
% Written by Matt Perich. Updated Feb 2017.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function combined = appendTDs(varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check the inputs
if isempty(varargin)
    error('No inputs provided!');
end
% check for any empty inputs
idx_empty = false(size(varargin));
for i = 1:length(varargin)
    if isempty(varargin{i})
        disp(['appendTDs: Input ' num2str(i) ' is empty. Skipping.']);
        idx_empty(i) = true;
    end
end
varargin = varargin(~idx_empty);
if isempty(varargin)
    error('All inputs were empty!')
end
% make sure they're all trial_data structs
if ~all(cellfun(@isstruct,varargin))
    error('All inputs must be trial_data struct!');
end

% now combine the  structs
combined = varargin{1};
for trial = 1:length(combined)
    combined(trial).stitch_marks = 1;
end

if length(varargin) > 1
    % check that all have same number of trials
    if length(unique(cellfun(@(x) length(x),varargin))) > 1
        error('All inputs must have same number of trials');
    end
    % get the list of time variables
    fn_time = getTDfields(combined,'time');
    fn_idx  = getTDfields(combined,'idx');
    
    for i = 2:length(varargin)
        td = varargin{i};
        for trial = 1:length(combined)
            combined(trial).stitch_marks = [combined(trial).stitch_marks, size(combined(trial).(fn_time{1}),1)];
            
            % append the time variables
            for var = 1:length(fn_time)
                combined(trial).(fn_time{var}) = [combined(trial).(fn_time{var}); td(trial).(fn_time{var})];
            end
            
            % adjust the idx_values
            for var = 1:length(fn_idx)
                if ~isempty(td(trial).(fn_idx{var}))
                    combined(trial).(fn_idx{var}) = td(trial).(fn_idx{var})+size(combined(trial).(fn_time{1}),1);
                end
            end
        end
    end
else
    disp('appendTDs: only one trial_data struct provided. Nothing to append.');
end

% ensure logical order
combined = reorderTDfields(combined);
