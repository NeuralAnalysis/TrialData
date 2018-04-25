function [trial_data,td_params] = convertDataToTD(varargin)
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
%       .bin_size        : (default: 0.01) the standard time vector
%       .trigger_thresh  : (default: 1) threshold for trigger signals
%
% OUTPUTS:
%   trial_data : the struct
%   td_params  : contains all parameters, including inputs
%
% Things to do:
%   - implement LFP. It's important to ensure indices match up with spikes!
%   - read digital data from cerebus NEV files
%   - deal with case where signals have lower sampling rate than bin_size.
%       Do we upsample those signals or increase the bin_size?
%   - Implement read_waveforms and include_spike_times
%
% Written by Matt Perich. Updated April 2018.
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
        %   duration : in s, how long is the file
        %   t_bin    : time vector
        %   data     : (rows are time, cols are signals)
        %   labels   : (cell array, one entry for each col of data)
        %   samprate : (sampling rate in Hz)
        %   meta     : (optional field, a struct holding any info to add)
        %
        % Note that data for spikes is different: it's a cell array of
        %  spike times of the same length as labels, rather than a matrix.
        % Also, events functions in this same way.
        %
        % Can be any time length, but first index for all signals (e.g. NEV,
        % NSx, and these, etc) has to all be aligned at the start so
        % this processing function has to handle syncing somehow
        % !!!!!BUT I AM FIXING THIS!!!!!
        file_data = which_routine(which_file,signal_info{iFile}.routine_params);
    else
        error('routine not specified!');
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
        
        if isfield(file_data_temp,'meta')
            sig_data(count).meta = file_data_temp.meta;
        end
        sig_data(count).name = which_name;
        sig_data(count).type = which_type;
        sig_data(count).duration = file_data_temp.duration;
        sig_data(count).samprate = file_data_temp.samprate;
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % parse out the labels associated with this signal
        % if multiple signals are specified, might want to do some
        % operation with them (e.g. compute joint angles, differential EMG)
        if ~isempty(which_label)
            if ~iscell(which_label), which_label = {which_label}; end
        end
        
        if ~isempty(which_label)
            if ischar(which_label{1})
                idx = ismember(file_data_temp.labels,which_label);
                if isempty(idx), error(['No label found! Here is the list: ' file_data_temp.labels]); end
            elseif length(which_label) == 1 % it's not a char, so it must be an array of numbers
                idx = which_label{1};
            else
                error('Cannot parse label. It appears to be a cell array of numbers. It should either be an array of numbers or a cell array of strings.');
            end
        else
            idx = 1:size(file_data_temp.data,2);
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
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Make all of these generic based on routine, and put this code in a
        % subfunc maybe?
        switch lower(which_type)
            case 'spikes'
                % for binning we need an extra on the end
                t_bin = (0:bin_size:file_data_temp.duration+bin_size)';
                data_bin = bin_spikes(data,1:length(data),t_bin);
                
            case 'emg' % filter appropriately and downsample
                t_bin = (0:bin_size:file_data_temp.duration)';
                data_bin = process_emg(data,samprate,params);
                
            case 'lfp' % filter into requested bands and downsample
                t_bin = (0:bin_size:file_data_temp.duration)';
                error('LFP not supported yet.');
                
            case 'trigger' % look for threshold crossing and call it an event
                t_bin = (0:bin_size:file_data_temp.duration)';
                % turn trigger into a timestamp
                trig_ts = find([0; diff(data > trigger_thresh) > 0])/file_data_temp.samprate;
                data_bin = histcounts(trig_ts,t_bin)'; % bin
                
            case 'event' % these are timestamps (in seconds)
                t_bin = (0:bin_size:file_data_temp.duration)';
                % will be cell if it's multiple labels, but otherwise
                if ~iscell(data), data = {data}; end
                data_bin = bin_events(data,t_bin);
                
            case 'meta' % it's a bit of a hack but meta goes in as if it's a signal
                t_bin = [];
                data_bin = data;
                
            case 'generic' % just downsample and call it a day
                t_bin = (0:bin_size:file_data_temp.duration)';
                data_bin = zeros(ceil(size(data,1)/round(bin_size*file_data_temp.samprate)),size(data,2));
                for i = 1:size(data,2)
                    data_bin(:,i) = decimate(double(data(:,i)),round(bin_size*file_data_temp.samprate));
                end
                
            otherwise
                error('type not recognized');
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

% group all EMG signals together
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
trial_data = trim_time(trial_data);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% find indices for trigger times to turn them into events
if ismember('trigger',file_types) || ismember('event',file_types) % triggers will need an extra "find"
    idx = find(strcmpi({sig_data.type},'trigger') | strcmpi({sig_data.type},'event'));
    for i = 1:length(idx)
        trial_data.(['idx_' sig_data(idx(i)).name]) = find(sig_data(idx(i)).data)';
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
function binned_emg = process_emg(data,samprate,params)
% Process EMG (high pass, rectify, low pass)
%   default: high pass at 10 Hz, rectify, low pass at 20 Hz
% filter
emg_LPF_cutoff  =  20;    % for EMG butterworth filter
emg_HPF_cutoff  =  10;    % for EMG butterworth filter
emg_n_poles     =  4;     % for EMG butterworth filter
if ~isempty(params), assignParams(who,params); end

[blow,alow] = butter(emg_n_poles,emg_LPF_cutoff/samprate);
[bhigh,ahigh] = butter(emg_n_poles,emg_HPF_cutoff/samprate,'high');
% !!! note the rectification step in the following command:
data = filtfilt(blow,alow,abs(filtfilt(bhigh,ahigh,double(data))));

binned_emg = zeros(ceil(size(data,1)/round(bin_size*samprate)),size(data,2));
for i = 1:size(data,2)
    binned_emg(:,i) = decimate(data(:,i),round(bin_size*samprate))';
end
clear temp_data;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function binned_events = bin_events(data,t_bin)
% here, the events start already as indices but we want to make sure they
% align with the timeframe provided, so we put them in bins again. This
% also helps streamline the processing code by reducing conditional
% processing on the data types, without a great hit to processing speed
% note that histcounts outputs rows
binned_events = zeros(length(data),length(t_bin));

% histcounts works weirdly and needs an extra bin
t_bin = [t_bin; t_bin(end)+mode(diff(t_bin))];

for e = 1:length(data)
    % get the event times for that cell in the current time window
    binned_events(e,:) = histcounts(data{e},t_bin);
end
% must transform to have same dimensions as kinematics etc
binned_events = binned_events';

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




