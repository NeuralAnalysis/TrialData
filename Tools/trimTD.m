%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function trial_data = trimTD(trial_data, varargin)
%
%   Will truncate all of the time-signals of each trial_data trial.
%
%   trial_data = trimTD(trial_data, params)
%       OR
%   trial_data = trimTD(trial_data, idx_start, idx_end)
%
%   With struct input, must define both idx_start and idx_end parameters.
%
% INPUTS:
%   trial_data : (struct) trial_data struct
%   idx_start  : (cell) {'idx_to_align_start',num_bins_after}
%   idx_end    : (cell) {'idx_to_align_end',num_bins_after}
%                   For the above, can put just 'start' or 'end' to do
%                   first or last bin available. Also, second entry isn't
%                   necessary. Will default to 0 for num_bins_after.
%
%   params:
%       zero_pad : if true, will zero pad signals to make all trials the
%           same length (default: false)
%       remove_short : if true, will remove the trials where requested
%           end time is after trial end (default: false)
% Note bin number for alignment can be negative to go before idx
%
% Written by Matt Perich. Updated July 2020.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function trial_data = trimTD(trial_data,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
idx_start = {};
idx_end   = {};
zero_pad = false;
remove_short = false;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Some undocumented extra parameters
verbose = true;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
trial_data = check_td_quality(trial_data);
if nargin == 2 && ~isstruct(varargin{1}), error('Must provide start and end points for trimming.'); end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin == 2 % it's a struct
    assignParams(who,varargin{1});
elseif nargin == 3 % start and end were defined
    idx_start = varargin{1};
    idx_end   = varargin{2};
else
    error('Too many input arguments.');
end
if ~iscell(idx_start), idx_start = {idx_start}; end
if ~iscell(idx_end), idx_end = {idx_end}; end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fn_time = getTDfields(trial_data,'time');
fn_idx = getTDfields(trial_data,'idx');

bad_idx = false(1,length(trial_data));
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
        error('Start input not formatted properly.');
    end
    
    if isempty(t_start)
        warning('No start time found. Excluding trial.');
        bad_idx(trial) = true;
    else
        
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
        
        
        % if there are multiple, assume the start needs the first idx and the
        % end needs the last idx
        if length(t_start) > 1
            disp('Multiple starting indices found. Using first...');
            t_start = t_start(1);
        end
        if length(t_end) > 1
            disp('Multiple ending indices found. Using last...');
            t_end = t_end(end);
        end
        
        if t_start < 1
            if verbose
                warning('Requested start time is < 1. Defaulting to first index...');
            end
            t_start = 1;
        end
        
        if t_end > t(end)
            if verbose
                warning('Trial %d: Requested end time went beyond trial time...',trial)
            end
            if ~zero_pad
                t_end = length(t);
            end
            if remove_short
                warning('Removing short trial (%d)',trial)
                bad_idx(trial) = true;
            end
        end
        
        if isnan(t_start) || isnan(t_end)
            error('Cannot trim, because some necessary indices are NaN.');
        end
        t_new = t_start:t_end;
        
        if t_start > t_end
            error('Start time is after end time. Cannot trim!');
        end
        
        % process time fields
        for iSig = 1:length(fn_time)
            temp = trial_data(trial).(fn_time{iSig});
            if isempty(temp) || ( length(temp) == 1 && isnan(temp) )
                if verbose
                    warning(['Time field ' fn_time{iSig} ' is empty or NaN. Propagating empty or NaN forward.']);
                end
                trial_data(trial).(fn_time{iSig}) = temp;
            elseif length(temp)<t_end
                if ~zero_pad
                    error('Something went wrong with the zero-pad thing Raeed added...Talk to him.')
                else
                    % zero pad temp to full length
                    temp_padded = zeros(t_end,size(temp,2));
                    temp_padded(1:length(temp),:) = temp;
                    temp = temp_padded;
                end
                
                trial_data(trial).(fn_time{iSig}) = temp(t_new,:);
            end
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
end

trial_data(bad_idx) = [];


