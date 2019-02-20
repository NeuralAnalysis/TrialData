function [trial_data,td_params,error_flag] = convertDataToTD(varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function provides a more flexible and experiment-agnostic framework
% to bring data into the trial_data world. Basic functionality for
% processing common signals (spiking neurons, EMG, soon LFP) is
% built-in, but it can be expanded to any arbitrary data type, format, or
% setup using the 'routine' and 'operation' inputs (see initSignalStruct
% for some usage details).
%
% INPUTS:
%   1) signal_info : see initSignalStruct
%   2) params      : (optional struct) allows you to overwrite defaults
%       .meta            : (struct) adds all fields as meta info to trials
%       .bin_size        : (default: 0.01) the standard time vector
%       .trigger_thresh  : (default: 1) threshold for trigger signals
%
% OUTPUTS:
%   trial_data : the struct
%   td_params  : contains all parameters, including inputs
%
% Things to do:
%   - implement LFP
%   - Implement read_waveforms and include_spike_times
%
% Written by Matt Perich. Updated December 2018.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEFAULT PARAMETERS
bin_size       =  0.01;    % start at 10ms bins by default
trigger_thresh  =  1;     % to identify time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Some parameters that CAN be overwritten, but are mostly intended to be
% hard coded and thus aren't documented in the function header
add_spike_times = false; % (not implemented) add unbinned spike times
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin == 1
    signal_info = varargin{1};
    params      = struct();
elseif nargin == 2
    signal_info = varargin{1};
    params      = varargin{2};
    if ~isempty(params)
        if ~isfield(params,'meta'), disp('WARNING: no meta information provided.'); end
        assignParams(who,params); % overwrite parameters
    end
elseif nargin > 2
    error('Too many input arguments.');
else
    error('No inputs provided.');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set things up
if ~iscell(signal_info), signal_info = {signal_info}; end
if add_spike_times, error('Spike times not implemented yet'); end
error_flag = false;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% loop along the files

% initalize by looking ahead at how many names
num_signals = sum(cellfun(@(x) length(x.name),signal_info));
sig_data = repmat(struct('type',[],'meta',[],'name',[],'duration',[],'samprate',[],'labels',[],'data',[],'t_bin',[]),1,num_signals);

% loop along the signal inputs
count = 1;
for iFile = 1:length(signal_info)
    which_file    = signal_info{iFile}.filename;
    which_routine = signal_info{iFile}.routine;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % process using specified routine
    if ~isempty(which_routine)
        % Inputs to function must be 1) filename, 2) params struct
        % Output needs format of struct with fields:
        %   t            : time vector. All signals will be aligned at 0! To
        %                  account for sync signals your code should make the
        %                  final sync/start time 0, and any time before should
        %                  be negative. Signals will be truncated to shortest
        %                  time duration
        %   cont_data    : rows are time, cols are signals)
        %   event_data   : cell array of event times
        %   cont_labels  : cell array, one entry for each col of cont_data
        %   event_labels : cell array, one entry for each cell of event_data
        %   meta         : optional field, a struct holding any info to add
        %
        % Note that data for spikes is different: it's a cell array of
        %  spike times of the same length as labels, rather than a matrix.
        % Also, events functions in this same way.
        %
        % If there is not simultaneous continuous and event, you can just
        % have the struct contain 'data' and 'labels' and it will determine
        % the correct way to use them based on the data type.
        %
        % If there's an error, add an error_flag field to the output and
        % this will all fail gracefully
        file_data = which_routine(which_file,signal_info{iFile});
    else
        error_flag = true;
        disp(['ERROR: ' mfilename ': routine not specified!']);
    end
    
    if isfield(file_data,'error_flag'), error_flag = file_data.error_flag; end
    
    if error_flag
        trial_data = [];
        td_params = [];
        if nargout == 3 % the user asked for the error flag so they probably will use it appropriately
            return;
        else
            error('There was an error! See the description above.');
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % get the signals which apply to this file
    for iSig = 1:length(signal_info{iFile}.name)
        
        % get the info we need from the struct
        which_name = signal_info{iFile}.name{iSig};
        which_type = signal_info{iFile}.type{iSig};
        which_label = signal_info{iFile}.label{iSig};
        which_operation = signal_info{iFile}.operation{iSig};
        
        % copy to make sure everything we do won't destroy the original
        file_data_temp = file_data;
        
        
        %%%%%%%%%%%%%%%%%%%%%%%
        % parse out the data contents if needed
        if isfield(file_data_temp,'cont_data') || isfield(file_data_temp,'event_data')
            if isfield(file_data_temp,'data')
                warning('The file data struct contains possibly redundant data fields. Ignoring ''data'' and using ''cont_data'' or ''event_data''.');
            end
            % determine if we are looking at continuous or event data and
            % reassigning things appropriately
            if strcmpi(which_type,'spikes') || strcmpi(which_type,'event')
                % it's an event type, so use event_data
                if isfield(file_data_temp,'event_data') && isfield(file_data_temp,'event_labels')
                    file_data_temp.data = file_data_temp.event_data;
                    file_data_temp.labels = file_data_temp.event_labels;
                    file_data_temp = rmfield(file_data_temp,{'event_data','event_labels'});
                else
                    error_flag = true;
                    disp('ERROR: File data struct was missing event_data or event_labels');
                end
            else % it's continuous
                % check to make sure the right fields are present
                if isfield(file_data_temp,'cont_data') && isfield(file_data_temp,'cont_labels')
                    file_data_temp.data = file_data_temp.cont_data;
                    file_data_temp.labels = file_data_temp.cont_labels;
                    file_data_temp = rmfield(file_data_temp,{'cont_data','cont_labels'});
                else
                    error_flag = true;
                    disp('ERROR: File data struct was missing cont_data or cont_labels');
                end
                
            end
            % if there was an error, end it all
            if error_flag
                trial_data = [];
                td_params = [];
                if nargout == 3 % the user asked for the error flag so they probably will use it appropriately
                    return;
                else
                    error('There was an error! See the description above.');
                end
            end
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if isfield(file_data_temp,'meta')
            sig_data(count).meta = file_data_temp.meta;
        end
        sig_data(count).name = which_name;
        sig_data(count).type = which_type;
        sig_data(count).duration = file_data_temp.t(end);
        sig_data(count).samprate = 1/mode(diff(file_data_temp.t)); % assumption: time vectors have uniform sampling
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % parse out the labels associated with this signal
        % if multiple signals are specified, might want to do some
        % operation with them (e.g. compute joint angles, differential EMG)
        if ~isempty(which_label)
            if ~iscell(which_label), which_label = {which_label}; end
        end
        
        if ~isempty(which_label)
            if ischar(which_label{1})
                if strcmpi(which_type,'spikes')
                    idx = 1:size(file_data_temp.labels,1);
                else
                    [~,idx] = ismember(which_label,file_data_temp.labels);
                end
                if any(idx==0)
                    error_flag = true;
                    disp(['ERROR: ' mfilename ': No label found! Here is the list: ' file_data_temp.labels]);
                end
            elseif length(which_label) == 1
                temp_label = which_label{1};
                if iscell(temp_label) && ischar(temp_label{1}) % it's multiple text entries
                    [~,idx] = ismember(temp_label,file_data_temp.labels);
                    if any(idx==0)
                        error_flag = true;
                        disp(['ERROR: ' mfilename ': No label found! Here is the list: ' file_data_temp.labels]);
                    end
                elseif isnumeric(temp_label) % it's not a char, so it must be an array of numbers
                    idx = temp_label;
                else
                    error_flag = true;
                    disp(['ERROR: ' mfilename ': no idea what to do with this label!']);
                end
            else
                error_flag = true;
                disp(['ERROR: ' mfilename ': Cannot parse label. It appears to be a cell array of numbers. It should either be an array of numbers or a cell array of strings.']);
            end
        else
            idx = 1:size(file_data_temp.data,2);
        end
        
        if error_flag
            trial_data = [];
            td_params = [];
            if nargout == 3 % the user asked for the error flag so they probably will use it appropriately
                return;
            else
                error('There was an error! See the description above.');
            end
        end
        
        
        data = file_data_temp.data(:,idx);
        if strcmpi(which_type,'spikes')
            sig_data(count).labels = file_data_temp.labels(idx,:);
        else
            sig_data(count).labels = file_data_temp.labels(idx);
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % do operation if requested
        if ~isempty(which_operation)
            % pass in a function handle that only takes data as input
            % and assumes that you work along the first row. THIS NEEDS
            % WORK OKAY IT'S A FIRST PASS DEAL WITH IT
            data = which_operation(data')';
        end
        
        % remove negative time values, thus aligning all signals
        switch sig_data(count).type
            case 'spikes' % it's a spike, so time
                for i = 1:length(data)
                    % shift them back in time by the amount post-sync
                    idx_keep = data{i} >= 0;
                    data{i} = data{i}(idx_keep);
                end
            case 'event'
                if iscell(data) % events are a time
                    for i = 1:length(data)
                        % shift them back in time by the amount post-sync
                        idx_keep = data{i} >= 0;
                        data{i} = data{i}(idx_keep);
                    end
                else % events are binned
                    idx_keep = file_data_temp.t >= 0;
                    data = data(idx_keep,:);
                end
            otherwise
                idx_keep = file_data_temp.t >= 0;
                data = data(idx_keep,:);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Make all of these generic based on routine, and put this code in a
        % subfunc maybe?
        
        % first get new unified time vector
        t_bin = (0:bin_size:sig_data(count).duration)';
        
        % now resample/bin the actual data
        switch lower(which_type)
            case 'spikes' % these are timestamps of spikes in seconds
                data_bin = bin_spikes(data,1:length(data),t_bin);
                
            case 'emg' % filter appropriately and resample/interpolate
                temp_params = params;
                temp_params.samprate = sig_data(count).samprate;
                data = process_emg(data,temp_params);
                
                % rebin the signals at the new sampling rate (given by bin_size)
                data_resampled = resample_signals(data,file_data_temp.t,struct('bin_size',bin_size,'samprate',sig_data(count).samprate));
                % then interpolate to unified time vector
                data_bin = interp1(t_resamp,data_resampled,t_bin);
                
            case 'lfp' % filter into requested bands and downsample
                error('LFP not supported yet.');
                
            case 'trigger' % look for threshold crossing and call it an event
                % turn trigger into a timestamp
                trig_ts = find([0; diff(data > trigger_thresh) > 0])/sig_data(count).samprate;
                data_bin = histcounts(trig_ts,t_bin)'; % bin
                
            case 'event' % these are timestamps (in seconds)
                % will be cell if it's multiple labels, but otherwise
                if ~iscell(data), data = {data}; end
                data_bin = bin_events(data,t_bin);

            case 'generic' % resample/interpolate
                % resample signals at new sampling rate
                [data_resampled,t_resamp] = resample_signals(data,file_data_temp.t,struct('bin_size',bin_size,'samprate',sig_data(count).samprate));
                % then interpolate to unified time vector
                data_bin = interp1(t_resamp,data_resampled,t_bin);
                
            otherwise
                error_flag = true;
                disp(['ERROR: ' mfilename ': type not recognized']);
        end
        
        if error_flag
            trial_data = [];
            td_params = [];
            if nargout == 3 % the user asked for the error flag so they probably will use it appropriately
                return;
            else
                error('There was an error! See the description above.');
            end
        end
        
        % package up data and move on
        sig_data(count).data = data_bin;
        sig_data(count).t_bin = t_bin;
        
        count = count + 1;
    end
end

% check the durations for all of them and make any missing ones the max
for iSig = 1:length(sig_data)
    if isempty(sig_data(iSig).duration)
        sig_data(iSig).duration = max(cell2mat({sig_data.duration}));
    end
end

% Also, check all of the names. Replace ' ' with '_' so that it doesn't
% throw errors on the struct field naming
for iSig = 1:length(sig_data)
    sig_data(iSig).name = strrep(sig_data(iSig).name,' ','_');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% we should now have a bunch of signals. Group appropriately and put in a
% master binned data file
file_types = unique({sig_data.type});

% check for spikes
if ismember('spikes',file_types)
    % there could be multiple arrays, but by design each spiking thing
    % should be only a single entry. Breaking this assumption will require
    % some additional coding, put it on the to do list!
    idx = find(strcmpi({sig_data.type},'spikes'));
    for i = 1:length(idx)
        trial_data.([sig_data(idx(i)).name '_spikes'])     = sig_data(idx(i)).data;
        trial_data.([sig_data(idx(i)).name '_unit_guide']) = sig_data(idx(i)).labels;
        trial_data.([sig_data(idx(i)).name '_spikes_t'])   = sig_data(idx(i)).t_bin;
    end
end

% group all EMG signals together
if ismember('emg',file_types)
    idx = find(strcmpi({sig_data.type},'emg'));
    trial_data.emg       = cat(2,sig_data(idx).data);
    trial_data.emg_names = {sig_data(idx).name};
    trial_data.emg_t = sig_data(idx(1)).t_bin;
end

% group all LFP signals together
if ismember('lfp',file_types)
    error('lfp is not implemented yet');
end

% generic signals are easy
if ismember('generic',file_types)
    idx = find(strcmpi({sig_data.type},'generic'));
    for i = 1:length(idx)
        trial_data.(sig_data(idx(i)).name) = sig_data(idx(i)).data;
        trial_data.([sig_data(idx(i)).name '_t']) = sig_data(idx(i)).t_bin;
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% trim the loaded data to be the shortest time vector
%   Note: right now assumes that all signals are on a unified bin_size
%   timeframe and that they all start in sync at 1
%   Also note: we do this before turning the events into bin counts. Here
%   they are still on a binned time vector
trial_data = trim_time(trial_data);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% find indices for trigger times to turn them into events
if ismember('trigger',file_types) || ismember('event',file_types) % triggers will need an extra "find"
    idx = find(strcmpi({sig_data.type},'trigger') | strcmpi({sig_data.type},'event'));
    for i = 1:length(idx)
        temp = find(sig_data(idx(i)).data)';
        % it's easier to see the numbers in the GUI if they are in a row
        if size(temp,1) > 1
            temp = temp';
        end
        trial_data.(['idx_' sig_data(idx(i)).name]) = temp;
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% now add in meta data from the sig_data entries
for iSig = 1:length(sig_data)
    if ~isempty(sig_data(iSig).meta)
        fn = fieldnames(sig_data(iSig).meta);
        for i = 1:length(fn)
            trial_data.(fn{i}) = sig_data(iSig).meta.(fn{i});
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Add in meta data to each trial
for iTrial = 1:length(trial_data)
    % loop along all meta fields
    if isfield(params,'meta')
        fields = fieldnames(params.meta);
        for fn = 1:numel(fields)
            trial_data(iTrial).(fields{fn}) = params.meta.(fields{fn});
        end
    end
    
    % and some defaults
    trial_data(iTrial).bin_size = bin_size;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ensure standard field order
trial_data = reorderTDfields(trial_data);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Package up parameter output
td_params = struct( ...
    'datetime', datetime, ...
    'bin_size', bin_size, ...
    'signal_info', signal_info);
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUB FUNCS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function binned_spikes = bin_spikes(units,unit_idx,t_bin)

% note that histcounts outputs rows
binned_spikes = zeros(size(unit_idx,1),length(t_bin));

% histcounts works weirdly and needs an extra bin
t_bin = [t_bin; t_bin(end)+mode(diff(t_bin))];

for unit = 1:length(unit_idx)
    % get the spikes for that cell in the current time window
    ts = units{unit_idx(unit)};
    binned_spikes(unit,:) = histcounts(ts,t_bin);
end
% must transform to have same dimensions as kinematics etc
binned_spikes = binned_spikes';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function binned_events = bin_events(data,t_bin)
% here, the events start already as indices but we want to make sure they
% align with the timeframe provided, so we put them in bins again. This
% also helps streamline the processing code by reducing conditional
% processing on the data types, without a great hit to processing speed
% note that histcounts outputs rows
binned_events = zeros(length(data),length(t_bin));

if ~isempty(data)
    if length(data{1}) == length(t_bin) % these are already bin idx
        for e = 1:length(data)
            if ~isempty(data{e})
                binned_events(e,:) = data{e};
            end
        end
        
    else % it's a timestamp and we need to bin it
        % histcounts works weirdly and needs an extra bin
        t_bin_temp = [t_bin; t_bin(end)+mode(diff(t_bin))];
        
        for e = 1:length(data)
            if ~isempty(data{e})
                %data_ts = t_bin(find(data{e}));
                % get the event times for that cell in the current time window
                binned_events(e,:) = histcounts(data{e},t_bin_temp);
            end
        end
        % must transform to have same dimensions as kinematics etc
        binned_events = binned_events';
    end
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function trial_data = trim_time(trial_data)

% find other possible time vectors
%   assumes format '_t'
fn = fieldnames(trial_data);
t_fn = fn(cellfun(@(x) ~isempty(x),strfind(fn,'_t')));

% get the shortest length
t_lengths = cellfun(@(x) size(trial_data.(x),1),t_fn);

t_min = min(t_lengths);

for i = 1:length(fn)
    % see if this time length is greater than min
    if size(trial_data.(fn{i}),1) > t_min
        temp = trial_data.(fn{i});
        trial_data.(fn{i}) = temp(1:t_min,:);
    end
end

% now remove all of the time fields
trial_data = rmfield(trial_data,t_fn);

end




