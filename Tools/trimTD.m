%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function trial_data = trimTD(trial_data, idx_start, idx_end)
%
%   Will truncate all of the time-signals of each trial_data trial.
%
% INPUTS:
%   trial_data : (struct) trial_data struct
%   idx_start  : (cell) {'idx_to_align_start',num_bins_after}
%   idx_end    : (cell) {'idx_to_align_end',num_bins_after}
%                   For the above, can put just 'start' or 'end' to do
%                   first or last bin available. Also, second entry isn't
%                   necessary. Will default to 0 for num_bins_after.
%
% Note bin number for alignment can be negative to go before idx
%
% Written by Matt Perich. Updated March 2017.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function trial_data = trimTD(trial_data,idx_start,idx_end)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 3, error('Must provide start and end points for trimming.'); end
if ~iscell(idx_start), idx_start = {idx_start}; end
if ~iscell(idx_end), idx_end = {idx_end}; end

fn_time = getTDfields(trial_data,'time');
fn_idx = getTDfields(trial_data,'idx');

for trial = 1:length(trial_data)
    % use any time signal to get the amount of time
    t = 1:size(trial_data(trial).(fn_time{1}),1);
    
    % parse the input to get the start idx
    if length(idx_start) == 2
        if strcmpi(idx_start{1},'start')
            t_start = 1+idx_start{2};
        elseif strcmpi(idx_start{1},'end')
            t_start = size(trial_data(trial).(fn_time{1}),1) + idx_start{2};
        else
            t_start = floor(trial_data(trial).(idx_start{1}) + idx_start{2});
        end
    elseif length(idx_start) == 1 && strcmpi(idx_start,'start')
        t_start = 1;
    elseif length(idx_start) == 1
        if ischar(idx_start{1}) % if it's a idx field name
            t_start = floor(trial_data(trial).(idx_start{1}));
        else % it's a numerical index
            t_start = int32(idx_start{1});
        end
        
    else
        error('End input not formatted properly.');
    end
    
    % parse the input to get the end idx
    if length(idx_end) == 2
        if strcmpi(idx_end{1},'start')
            t_end = 1+idx_end{2};
        elseif strcmpi(idx_end{1},'end')
            t_end = size(trial_data(trial).(fn_time{1}),1) + idx_end{2};
        else
            t_end = ceil(trial_data(trial).(idx_end{1}) + idx_end{2});
        end
    elseif length(idx_end) == 1 && strcmpi(idx_end,'end')
        t_end = size(trial_data(trial).(fn_time{1}),1);
    elseif length(idx_end) == 1
        if ischar(idx_end{1}) % if it's a idx field name
            t_end = ceil(trial_data(trial).(idx_end{1}));
        else % it's a numerical index
            t_end = int32(idx_end{1});
        end
    else
        error('End input not formatted properly.');
    end
    
    if t_end > t(end)
        warning('Requested end time went beyond trial time...')
        t_end = length(t);
    end
    if isnan(t_start) || isnan(t_end)
        error('Cannot trim, because some necessary indices are NaN.');
    end
    t_new = t_start:t_end;
    
    % process time fields
    for iSig = 1:length(fn_time)
        temp = trial_data(trial).(fn_time{iSig});
        trial_data(trial).(fn_time{iSig}) = temp(t_new,:);
    end
    
    % process idx fields
    for iIdx = 1:length(fn_idx)
        temp = trial_data(trial).(fn_idx{iIdx});
        if temp > t_end, temp = NaN; end
        if temp < 0, temp = NaN; end
        if ~isnan(temp)
            % in cases like the RW go cues, there can be multiple idx_, so
            % loop along them
            temp_adjust = zeros(size(temp));
            for i = 1:size(temp,2)
                if isempty(temp(i)) || (temp(i) < t_new(1) || temp(i) > t_new(end))
                    temp_adjust(i) = NaN;
                else
                    temp_adjust(i) = find(t_new <= temp(i),1,'last');
                end
            end
            trial_data(trial).(fn_idx{iIdx}) = temp_adjust;
        else
            trial_data(trial).(fn_idx{iIdx}) = NaN; % Should this be NaN or []?
        end
        
    end
    
    
end

% if you're trimming, it's no longer continuous
if isfield(trial_data,'is_continuous')
    trial_data = rmfield(trial_data,'is_continuous');
end

