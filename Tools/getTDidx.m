function idx = getTDidx(trial_data,varargin)
% this function will return indices for trials in trial_data that meet the
% criteria defined by varargin. varargin should be in pairs, where the
% first entry is the field name and the second is the desired value.
%
%  e.g. to get all baseline trials:
%       idx = get_trial_data_indices(trial_data,'epoch','BL');
%  e.g. to get all force field trials to the pi/2 direction
%       idx = get_trial_data_indices(trial_data,'epoch','AD','perturbation','FF','target_direction',pi/2)

if rem(length(varargin),2) ~= 0
    error('Inputs must be provided in pairs stating ...'' VARIABLE '',''VALUE'',...');
end

idx = ones(size(trial_data));

for i = 1:2:length(varargin)
    if ischar(varargin{i+1})
        idx = idx & strcmpi({trial_data.(varargin{i})},varargin{i+1});
    elseif iscell(varargin{i+1})
        idx = idx & ismember({trial_data.(varargin{i})},varargin{i+1});
    else
        idx = idx & [trial_data.(varargin{i})] == varargin{i+1};
    end
end

