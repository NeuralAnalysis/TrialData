%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function trial_data = dupeAndShift(trial_data,varargin)
%
%   Will duplicate any specified variables and shift them some number of
% bins. This is mainly to let you preserve history while cutting and
% pasting chunks and data.
%
% Creates new fields in the trial_data struct of format:
%       [ORIGINAL_VAR_NAME '_shift']
% And the size of this new array will be:
%       ORIGINAL_SIZE * NUMBER OF SHIFTS
%
% Will truncate to remove all non-overlapping segments, and will adjust the
% idx fields to reflect that. Note you can only shift back in time for now.
%
% INPUTS:
%   trial_data : (struct) obvious
%   varargin   : pairs of ...'field',SHIFT,...
%       SHIFT is given in number of bins
%       List of possible variables are hard coded for now for stability
%
% OUTPUTS:
%   trial_data : original struct with added _shift fields
%
% Written by Matt Perich. Updated Feb 2017.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function trial_data = dupeAndShift(trial_data,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
idx_end_name = 'idx_trial_end';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if length(varargin) == 1 && iscell(varargin)
    varargin = varargin{1};
else
    if rem(length(varargin),2) ~= 0
        error('Must give pairs of inputs as ...field,shift,...');
    end
end
fn = fieldnames(trial_data);
the_shifts = [varargin{2:2:end}];
which_fields = {varargin{1:2:end}};

if any(~ismember(which_fields,fn)), error('Field not recognized'); end

[the_shifts,sort_idx] = sort(the_shifts,2,'Ascend');
which_fields = which_fields(sort_idx);

% get the max, as all signals will be truncated to this
max_shift = max(the_shifts);

% get names of the kinematic and spiking fields
fn_time = getTDfields(trial_data,'time');
fn_idx = getTDfields(trial_data,'idx');

if ~isfield(trial_data,'is_continuous') || ~any([trial_data.is_continuous]) % normal trial-based operation
    for trial = 1:length(trial_data)
        for j = 1:length(which_fields)
            if the_shifts(j) > 0
                temp = trial_data(trial).(which_fields{j});
                
                temp_shift = NaN(size(temp,1)+max_shift,size(temp,2)*(1+the_shifts(j)));
                temp_shift(1:end-max_shift,1:size(temp,2)) = temp;
                for k = 1:the_shifts(j)
                    temp_shift(k+1:k+size(temp,1),1+size(temp,2)*k:size(temp,2)*(k+1)) = temp;
                end
                
                % remove padding
                temp_shift = temp_shift(max_shift+1:end-max_shift,:);
                trial_data(trial).([which_fields{j} '_shift']) = temp_shift(:,size(temp,2)+1:end);
            else
                warning('You gave me a shift <= 0...');
            end
        end
        
        % remove extra time from other time varying signals
        for j = 1:length(fn_time)
            temp = trial_data(trial).(fn_time{j});
            trial_data(trial).(fn_time{j}) = temp(max_shift+1:end,:);
        end
        
        % remove extra time from index to keep all signals on same timeframe
        for j = 1:length(fn_idx)
            trial_data(trial).(fn_idx{j}) = trial_data(trial).(fn_idx{j}) - max_shift;
            if trial_data(trial).(fn_idx{j}) <= 0, trial_data(trial).(fn_idx{j}) = NaN; end
        end
    end
elseif isfield(trial_data,'is_continuous') && all([trial_data.is_continuous]) % preserve continuous data
    
    for j = 1:length(which_fields)
        if the_shifts(j) > 0
            temp = cat(1,trial_data.(which_fields{j}));
            trial_starts = cumsum([1,cellfun(@(x) size(x,1),{trial_data.(fn_time{1})})]);
            
            temp_shift = NaN(size(temp,1)+max_shift,size(temp,2)*(1+the_shifts(j)));
            temp_shift(1:end-max_shift,1:size(temp,2)) = temp;
            for k = 1:the_shifts(j)
                temp_shift(k+1:k+size(temp,1),1+size(temp,2)*k:size(temp,2)*(k+1)) = temp;
            end
            
            % now break it back into trials
            for trial = 1:length(trial_data)
                idx = trial_starts(trial):trial_starts(trial+1)-1;
                trial_data(trial).([which_fields{j} '_shift']) = temp_shift(idx,size(temp,2)+1:end);
            end
            
        else
            warning('You gave me a shift <= 0...');
        end
    end
    
    fn_time = getTDfields(trial_data,'time');
    % only need to remove time from the beginning of the first trial
    for j = 1:length(fn_time)
        temp = trial_data(1).(fn_time{j});
        trial_data(1).(fn_time{j}) = temp(max_shift+1:end,:);
    end
    
    % remove extra time from index to keep all signals on same timeframe
    for j = 1:length(fn_idx)
        trial_data(1).(fn_idx{j}) = trial_data(1).(fn_idx{j}) - max_shift;
        if trial_data(1).(fn_idx{j}) <= 0, trial_data(1).(fn_idx{j}) = NaN; end
    end
    
else
    error('Some trials are from continuous data, some are not. dupeAndShift is confused.');
end

% restore logical order
trial_data = reorderTDfields(trial_data);

