%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [return_idx, trial_data] = getTDidx(trial_data, varargin)
%
%   This function will return indices for trials in trial_data that meet
% the criteria defined by varargin. varargin should be in pairs, where the
% first entry is the field name and the second is the desired value. It's
% pretty much agnostic to what fields trial_data contains, at least for all
% of the use cases that I've come across so far.
%
% Basically, first entry is the field to use. Second is an appropriate
% piece of information. There are a couple more useful things:
%   1) Use ...'range',[start,end]... to get a range of consecutive trials
%       that meet all other criteria. If start and end are integers, will
%       return everything in between (ends included). However, if start/end
%       are fractions, it will treat them like percentages
%   2) Use ...'rand',N... to get N random trials that meet all other criteria
%       (including the range specified above)
%
% OUTPUTS:
%   return_idx : list of indices that meet the desired criteria
%   trial_data : the trials of the struct that correspond to return_idx
%
% EXAMPLES:
%  e.g. to get all baseline trials:
%       idx = getTDidx(trial_data,'epoch','BL');
%  e.g. to get all force field trials to the pi/2 direction
%       idx = getTDidx(trial_data,'epoch','AD','perturbation','FF','target_direction',pi/2)
%  e.g. to get 10 random trials from the first half of baseline
%       idx = getTDidx(trial_data,'epoch','BL','range',[0 0.5],'rand',10);
%
% Written by Matt Perich. Updated Feb 2017.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [return_idx, trial_data] = getTDidx(trial_data,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if rem(length(varargin),2) ~= 0
    error('Inputs must be provided in pairs stating ...'' VARIABLE '',''VALUE'',...');
end

fn = varargin(1:2:length(varargin));
fv = varargin(2:2:length(varargin));

idx = ones(size(trial_data));

% Check for the trials that match each criterion
for i = 1:length(fn)
    if ~strcmpi(fn{i},'rand') && ~strcmpi(fn{i},'range')
        if ischar(fv{i})
            idx = idx & strcmpi({trial_data.(fn{i})},fv{i});
        elseif iscell(fv{i})
            idx = idx & ismember({trial_data.(fn{i})},fv{i});
        else
            idx = idx & [trial_data.(fn{i})] == fv{i};
        end
    end
end

return_idx = find(idx);

% see if a range was requested and return it
idx = strcmpi(fn,'range');
if any(idx)
    bounds = fv{idx};
    if bounds(1) >= bounds(2)
        warning('Make range monotonic increasing. Flipping order.');
        bounds = fliplr(bounds);
    end
    if bounds(2) <= 1 % it's a percent
        if bounds(1) < 0
            warning('No negative percent. Set to 0.');
            bounds(1) = 0;
        end
        return_idx = return_idx( 1+floor(bounds(1)*length(return_idx)) : floor(bounds(2)*length(return_idx)) );
    else
        if bounds(1) > length(return_idx)
            warning('Requested lower bound too large. No trials will match');
        end
        if bounds(2) > length(return_idx)
            warning('Requested upper bound too large. Returning max instead.');
            bounds(2) = length(return_idx);
        end
        return_idx = return_idx(bounds(1):bounds(2));
    end
end

% see if a random return is requested
idx = strcmpi(fn,'rand');
if any(idx)
    n = fv{idx};
    if n > length(return_idx)
        warning('Requested too many random trials. Receiving all');
        n = length(return_idx);
    end
    if n < 1
        warning('Requested too few random trials. Receiving one.');
        n = 1;
    end
    % get random trials
    temp = randperm(length(return_idx));
    return_idx = return_idx(temp(1:n));
end

% prepare outputs, in case this is desired
trial_data = trial_data(return_idx);

