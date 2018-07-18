%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function td_s = splitTD(trial_data,params)
%
%   Separates trials based on arbitrary events. Can take some time
% before/after each go cue, but note that this will result in overlapping
% time points across trials.
%
% Will also add idx_trial_start and idx_trial_end fields
%
% Note that you can just a string in as params which refers to the split
% idx, saving the need to define a struct, and it will just use the default
% params for everything else.
%
% INPUTS:
%   trial_data : the struct
%   params     : parameter struct
%       .split_idx_name    : (string) name of idx to split on
%                               (default: 'idx_go_cue')
%       .linked_fields     : (cell of strings) which fields are linked
%                               to the split idx (e.g. a trial RESULT info
%                               tied to each TRIAL_START event)
%       .extra_bins        : [TIME_BEFORE, TIME_AFTER] (in # bins, must be positive)
%
% Note: if the length(trial_data) == 1, it will add a "is_continuous" field
% to each trial to signify that this was a single piece of data that was
% split
%
% OUTPUTS:
%   td_s : the struct separated by movements
%
% Written by Matt Perich. Updated February 2018.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function td_s = splitTD(trial_data,params)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PARAMETER DEFAULTS
split_idx_name     =  'idx_trial_start';
linked_fields      =  {}; % list of meta etc fields linked to split_idx
extra_bins         =  [0 0];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Some parameters to overwrite that aren't documented
start_name         =  'idx_trial_start';
end_name           =  'idx_trial_end';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin > 1
    if ischar(params) % someone just passed in a string name. Very convenient
        split_idx_name = params;
    elseif isstruct(params)
        assignParams(who,params);
    else
        error('Not sure what to do with this params input.');
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isstruct(trial_data), error('First input must be trial_data struct!'); end
if ~iscell(linked_fields), linked_fields = {linked_fields}; end
if any(extra_bins < 0)
    disp('Extra bins must be positive: [BINS_BEFORE, BINS_AFTER]. Taking absolute value...');
    extra_bins = abs(extra_bins);
end

num_moves = sum(~isnan([trial_data.(split_idx_name)]));
td_s = repmat(struct(),1,num_moves);

fn_time = getTDfields(trial_data,'time');
fn_array = getTDfields(trial_data, 'arrays');
fn_idx = getTDfields(trial_data,'idx');
% remove the requested idx field
fn_idx = setdiff(fn_idx,split_idx_name);
fn_meta = getTDfields(trial_data,'meta');
% remove the linked fields from meta
fn_meta = setdiff(fn_meta,linked_fields);

if length(trial_data) == 1
    is_continuous = true;
else
    is_continuous = false;
end
count = 0;
for trial = 1:length(trial_data)
    td = trial_data(trial);
    split_idx = td.(split_idx_name);
    
    t_max = size(td.(fn_time{1}),1);
    
    
    for idx = 1:length(split_idx)
        if ~isnan(split_idx(idx)) && split_idx(idx) <= t_max
            count = count + 1;
            
            % first split up attached fields
            if ~isempty(linked_fields)
                for i = 1:length(linked_fields)
                    if length(td.(linked_fields{i})) == length(split_idx)
                        temp = td.(linked_fields{i});
                        temp = temp(idx);
                        if iscell(temp), temp = temp{1}; end
                        td_s(count).(linked_fields{i}) = temp;
                    end
                end
            end
            
            % copy over the meta data
            for i = 1:length(fn_meta)
                td_s(count).(fn_meta{i}) = td.(fn_meta{i});
            end
            
            % get the start and end indices
            idx_start = split_idx(idx)-extra_bins(1);
            if idx < length(split_idx)
                idx_end = split_idx(idx+1)+extra_bins(2)-1;
            else
                idx_end = Inf;
            end
            
            % check that they aren't outside the bounds of time
            if idx_start < 0, idx_start = 0; end
            if idx_end > t_max, idx_end = t_max; end
            
            % now add the new idx
            td_s(count).(split_idx_name) = extra_bins(1)+1;
            for i = 1:length(fn_idx)
                % find idx that are within idx_start and idx_end
                temp = td.(fn_idx{i});
                td_s(count).(fn_idx{i}) = temp(temp >= idx_start & temp < idx_end) - idx_start;
            end
            
            td_s(count).(start_name) = extra_bins(1)+1;
            td_s(count).(end_name)   = idx_end - extra_bins(2) - split_idx(idx)+1;
            
            % now add array information
            for i = 1:length(fn_array)
                td_s(count).([fn_array{i}, '_unit_guide']) = trial_data(1).([fn_array{i}, '_unit_guide']);
            end
            % check that the index won't crash
            if idx_end > size(td.(fn_time{1}),1)
                disp('Requested time extended beyond available trial data. Defaulting to last bin.');
                idx_end = size(td.(fn_time{1}),1);
            end
            % add time signals
            for i = 1:length(fn_time)
                temp = td.(fn_time{i});
                td_s(count).(fn_time{i}) = temp(idx_start:idx_end,:);
            end
            if is_continuous
                td_s(count).is_continuous = true;
            end
        end
    end
end

% find bad trials
td_s(cellfun(@(x) size(x,1)==0,{td_s.(fn_time{1})})) = [];

td_s = reorderTDfields(td_s);

