%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function trial_data = dupeAndShift(trial_data,varargin)
%
%   Will duplicate any specified variables and shift them some number of
% bins. This is mainly to let you preserve history while cutting and
% pasting chunks and data.
%
%   KEEP THIS IN  MIND:
%       -> Negative shift values slide the signal backwards, effectively
%           adding HISTORY to the current time bin
%       -> Positive shift values slide the signal forwards, effectively
%           adding THE FUTURE to the current time bin
%
% Creates new fields in the trial_data struct of format:
%       [ORIGINAL_VAR_NAME '_shift']
% And the size of this new array will be:
%       ORIGINAL_SIZE * NUMBER OF SHIFTS
%
% Will truncate to remove all non-overlapping segments, and will adjust the
% idx fields to reflect that.
%
% INPUTS:
%   trial_data : (struct) obvious
%   varargin   : pairs of ...'field',SHIFT,...
%       SHIFT is given in number of bins
%   params     : (struct) must be the last input!
%
% OUTPUTS:
%   trial_data : original struct with added _shift fields
%
% Written by Matt Perich. Updated March 2019.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function trial_data = dupeAndShift(trial_data,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% some extra parameters that can be overwritten
field_extra =  '_shift';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
params = [];
if length(varargin) == 1 && iscell(varargin)
    varargin = varargin{1};
else
    if rem(length(varargin),2) ~= 0
        % check if the last one is a params struct
        if isstruct(varargin{end})
            params = varargin{end};
            varargin = varargin(1:end-1);
        else
            error('Must give pairs of inputs as ...field,shift,...');
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isempty(params)
    assignParams(who,params);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
trial_data  =  check_td_quality(trial_data);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fn = fieldnames(trial_data);
all_shifts = {varargin{2:2:end}};
signals = {varargin{1:2:end}};
signals = check_signals(trial_data,signals);
signals = signals(:,1); % we don't need the idx
field_extra =  check_field_extra(field_extra,signals');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if all(~ismember(signals,fn)), error('Field not recognized'); end

% check the  shifts. They must all shift forward or all shift back in time
if ~( all(cellfun(@(x) all(x > 0), all_shifts)) || all(cellfun(@(x) all(x < 0), all_shifts)) )
    error('Currently shifts must be all  positive or all negative...');
end

% build a flag to say if it's positive or negative
direction_flag = all(cellfun(@(x) all(x > 0), all_shifts));

% get the max, as all signals will be truncated to this
max_shift = max(cellfun(@(x) max(abs(x)), all_shifts));

% get names of the kinematic and spiking fields
fn_time = getTDfields(trial_data,'time');
fn_idx = getTDfields(trial_data,'idx');


% loop  along trials

for trial = 1:length(trial_data)
    
    for j = 1:length(signals)
        the_shifts = all_shifts{j};
        the_field  = signals{j};
        [the_shifts,~] = sort(the_shifts,2,'Ascend');
        
        n_shifts =  length(the_shifts);
        
        temp = trial_data(trial).(the_field);
        
        if direction_flag % shift backward in time (i.e. add future)
            temp_shift = NaN(size(temp,1),size(temp,2)*(1+n_shifts));

            for k = 1:length(the_shifts)
                temp_shift(1:size(temp,1)- the_shifts(k), 1+size(temp,2)*(k):size(temp,2)*(k+1)) = temp(the_shifts(k)+1:end,:);
            end
            
            trial_data(trial).([the_field field_extra{j}]) = temp_shift(:,size(temp,2)+1:end);
            
        else % shift forward in time (i.e.add history)
            temp_shift = NaN(size(temp,1)+max_shift,size(temp,2)*(1+n_shifts));
            the_shifts_new  = sort(-the_shifts);
            for k = 1:length(the_shifts_new)
                temp_shift(the_shifts_new(k)+1:size(temp,1),1+size(temp,2)*k:size(temp,2)*(k+1)) = temp(1:end-the_shifts_new(k),:);
            end
            
            % remove padding
            temp_shift = temp_shift(1:end-max_shift,:);
            trial_data(trial).([signals{j} field_extra{j}]) = temp_shift(:,size(temp,2)+1:end);
            
        end
        
        
        trial_data(trial).([signals{j} field_extra{j} '_vals']) = the_shifts;
    end
    
    
    % update list of the kinematic and spiking fields
    fn_time = getTDfields(trial_data(trial),'time');
    
    % remove extra time from other time varying signals
    if direction_flag
        for k = 1:length(fn_time)
            temp = trial_data(trial).(fn_time{k});
            trial_data(trial).(fn_time{k}) = temp(1:size(temp,1)-max_shift,:);
        end
        
        % remove extra time from index to keep all signals on same timeframe
        for k = 1:length(fn_idx)
            temp = trial_data(trial).(fn_idx{k});
            temp( temp > size(trial_data(trial).(fn_time{1}),1)) = NaN;
            trial_data(trial).(fn_idx{k}) = temp;
        end
        
    else
        for k = 1:length(fn_time)
            temp = trial_data(trial).(fn_time{k});
            trial_data(trial).(fn_time{k}) = temp(max_shift+1:end,:);
        end
        
        % make sure idx are not larger
        for k = 1:length(fn_idx)
            temp = trial_data(trial).(fn_idx{k});
            temp = temp  - max_shift;
            temp( temp < 1 ) = NaN;
            trial_data(trial).(fn_idx{k}) = temp;
        end
        
    end
    
end


% restore logical order
trial_data = reorderTDfields(trial_data);

