%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [trial_data,td_params] = parseFileByTrial(data, params)
%
%   Wrapper function to create trial_data structs. Currently only supports
% CDS.
% 
% NOTE: Discretizes everything to the dt of the CDS kinematics (by default
% should be 10 ms bins, I believe). Can rebin later if desired using binTD.
%
% INPUTS:
%   data     : a CDS file or BDF file (not yet supported)
%   params   : a struct containing parameters
%     .meta          : a struct with a field for each meta parameter you want attached
%                       to this file. This can handle any arbitrary information!
%     .event_list   : Which cds.trials events to add to struct
%                       Format: {'CDS_TRIAL_TABLE_NAME','ALIAS'; ... etc ...}
%                       Can ignore ALIAS and give Nx1 cell vector if you want
%                       By default, assumes startTime and endTime exist,
%                       and attempts to add tgtOnTime and goCueTime if possible
%     .array_alias   : Aliases for renaming arrays from CDS names
%                       Format: {'CDS_NAME','NEW_NAME'; ...etc...}
%     .exclude_units : ID for which units to exclude (Default: [0,255])
%                       NOTE: this default gets rid of unsorted!
%     .trial_results : which reward codes to use ('R','A','F','I') (Default: {'R'})
%     .extra_time    : [time before, time after] beginning and end of trial (default [0.2 0.2] sec)
%     .all_points    : flag to include all data points. Thus, disregards extra_time
%                       and trial_results. Each trial ends at trial_start of the one after
%     .pos_offset    : offset (in units of cds.pos) to zero position (default [0,0])
%     .include_ts    : Flag to denote whether to include timestamps aligned
%                       start of the trial or not. If so, includes a cell
%                       array of timestamps for each trial, one for each
%                       source of neural data.
%     .include_start :  Flag to denot whether to include an extra column
%                       which contains the 'real time' start for each trial
%                       Useful for comparing in non-standard tasks
%       
%     .include_naming: Flag to denote whether to include a conversion
%                       matrix between the labels that appear on the screen when you are doing sensory mapping
%                       versus the labels that are given to the recorded
%                       file.
% OUTPUTS:
%   trial_data : the struct! Huzzah!
%
% TIPS:
%   I recommend the following extra function calls afterwards
%       trial_data = getMoveOnsetAndPeak(trial_data);
%       trial_data = removeBadTrials(trial_data);
%       trial_data = getCommonUnits(trial_data);
%
% Written by Matt Perich. Updated Feb 2017.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [trial_data,td_params] = parseFileByTrial(data,params)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin == 1, params = []; end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isa(data,'commonDataStructure')
    [trial_data,td_params] = parseFileByTrial_cds(data,params);
else
    error('BDF not currently supported.');
    % when someone is ready to implement this, can write the following func
    %   trial_data = parseFileByTrial_bdf(data,params);
end

% ensure standard field order
trial_data = reorderTDfields(trial_data);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [trial_data,td_params] = parseFileByTrial_cds(cds,inputArgs)
%
% INPUTS:
%   cds    : CDS object
%   params : a struct containing parameters (see above)
%
% OUTPUTS:
%   trial_data : the trial_data struct
%
% Written by Matt Perich. Updated Feb 2017.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [trial_data,td_params] = parseFileByTrial_cds(cds,params)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEFAULT PARAMETERS
event_list     =  {};
array_alias    =  {};
trial_results  =  {'R'};
exclude_units  =  [0,255];
extra_time     =  [0.2, 0.2];
all_points     =  false;
pos_offset     =  [0,0];
include_ts     =  false;
include_naming =  false;
marker_data    =  [];
if ~isfield(params,'meta'), disp('WARNING: no meta information provided.'); end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Some parameters that CAN be overwritten, but are mostly intended to be
% hard coded and thus aren't documented in the function header
LPF_cutoff  =  20; % for EMG butterworth filter
HPF_cutoff  =  10; % for EMG butterworth filter
n_poles     =  4;  % for EMG butterworth filter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin > 1 && ~isempty(params)
    assignParams(who,params); % overwrite parameters
else
    params = struct();
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do some input processing
if ~iscell(trial_results), trial_results = {trial_results}; end
switch size(event_list,2)
    case 0, event_alias = {};
    case 1, event_alias = {};
    case 2, event_alias = event_list; event_list = event_list(:,1);
    otherwise, error('event_list should have two columns {name, alias (if desired)}');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% see what array data is present
if ~isempty(cds.units)
    arrays = strsplit(cds.meta.array,', ');
    % get info on neurons
    unit_idx = cell(1,length(arrays));
    for array = 1:length(arrays)
        unit_idx{array} = find(~ismember([cds.units.ID],exclude_units) & strcmpi({cds.units.array},arrays{array}));
        if include_naming
            chanNames = cds.units(~cellfun(@isempty,([strfind({cds.units.array},arrays{array})])));
            sortedUnits = chanNames([chanNames.ID]>0 & [chanNames.ID]<255);
            elecNames = unique([sortedUnits.chan]);
            screenNames = {sortedUnits.label};
            labelNames = zeros(length(sortedUnits),1);
            for i= 1:length(sortedUnits)
               labelNames(i) = str2num(screenNames{i}(5:end)); 
            end
            labels = unique(labelNames);
            conversion{array} = [elecNames', labels];
        end
    end
end


% process some of the trial information
fn = cds.trials.Properties.VariableNames;
if ~all(ismember({'startTime','endTime'},fn)), error('Must have start and end times in CDS.'); end
event_list = union({'startTime';'endTime'},event_list);
if ismember({'goCueTime'},fn), event_list = union({'goCueTime'},event_list); end
if ismember({'tgtOnTime'},fn), event_list = union({'tgtOnTime'},event_list); end
% determine which signals are time-varying and which are parameter values
%   There was a CDS bug where start/end times didn't have units, but I know
%   they are supposed to be here so it's hard coded for now
%   Some parameter values have 's' as units but are not time events.
%   Include them in the time_event_exceptions below.
time_event_exceptions = {'ctrHold'};
time_event_exc_idx = false(size(fn));
for exc = 1:length(time_event_exceptions)
    time_event_exc_idx = time_event_exc_idx | contains(fn,time_event_exceptions{exc});
end
extra_time_events = fn( strcmpi(cds.trials.Properties.VariableUnits,'s') & ~time_event_exc_idx );
time_events = union({'startTime','endTime'},extra_time_events);

% get trial list and initialize
if all_points % we want everything
    idx_trials = 1:length(cds.trials.result);
else
    idx_trials = find(ismember(cds.trials.result,trial_results));
end
trial_data = repmat(struct(),1,length(idx_trials));

% find the bin size of the CDS kinematics
bin_size = round(1000*mode(diff(cds.kin.t)))/1000;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Add basic data and metadata
for i = 1:length(idx_trials)
    iTrial = idx_trials(i);
    % add some meta data about the trial
    trial_data(i).monkey = cds.meta.monkey;
    trial_data(i).date = datestr(cds.meta.dateTime,'mm-dd-yyyy');
    trial_data(i).task = cds.meta.task;
    
   
    switch lower(cds.meta.task)
        case 'co' % center out apparently has broken target dir, so use target id
            % this assumes there are exactly eight targets!!!
            % some CDS trials have wrong vals which is a bug likely
            targ_angs = [pi/2, pi/4, 0, -pi/4, -pi/2, -3*pi/4, pi, 3*pi/4];
            if ~isnan(cds.trials.tgtID(iTrial)) && cds.trials.tgtID(iTrial) < length(targ_angs)
                trial_data(i).target_direction = targ_angs(cds.trials.tgtID(iTrial)+1);
            else
                disp('This CO trial had a bad target ID...');
                trial_data(i).target_direction = NaN;
            end
        case 'rw' % In random walk, target_direction doesn't make sense
            trial_data(i).target_center = reshape(cds.trials.tgtCtr(iTrial,:),size(cds.trials.tgtCtr(iTrial,:),2)/2,2);
        case 'trt' % TRT denotes targets in the same way as random walk
            trial_data(i).target_center = reshape(cds.trials.tgtCtr(iTrial,:),2,size(cds.trials.tgtCtr(iTrial,:),2)/2)';
        otherwise
            if any(abs(cds.trials.tgtDir) > 2*pi) % good assumption that it's deg
                trial_data(i).target_direction = pi/180*cds.trials.tgtDir(iTrial);
            else % should be rad
                trial_data(i).target_direction = minusPi2Pi(cds.trials.tgtDir(iTrial));
            end
    end
    
    trial_data(i).trial_id = iTrial;
    trial_data(i).result = cds.trials.result(iTrial);
    trial_data(i).bin_size = bin_size;
    
    % if it's all_points, add a flag for later use
    if all_points
        trial_data(i).is_continuous = true;
    end
    
    % loop along all meta fields
    if isfield(params,'meta')
        fields = fieldnames(params.meta);
        for fn = 1:numel(fields)
            trial_data(i).(fields{fn}) = params.meta.(fields{fn});
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bin CDS <- Can add support for any signal in this section
%   1) Check for kinematics and bin those
%   2) Check for force and bin that
%   3) Check for EMG and process/bin those
%   4) Check for OpenSim data and bin those
%   5) Check for marker data and bin those
%   6) Turn the trial table into a binned version
kin_list = {'t','x','y','vx','vy','ax','ay'};
% force_list = {'t','fx','fy'};
force_list = cds.force.Properties.VariableNames;
cds_bin = struct();
cds_bin.kin = decimate_signals(cds.kin,kin_list,bin_size);
cds_bin.force = decimate_signals(cds.force,force_list,bin_size);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Process EMG (high pass, rectify, low pass)
%   default: high pass at 10 Hz, rectify, low pass at 20 Hz
if ~isempty(cds.emg)
    emg_list = cds.emg.Properties.VariableNames;
    
    emg=cds.emg;
    % find sampling rate
    samprate = 1/mode(diff(emg.t));
    % filter
    [blow,alow] = butter(n_poles,LPF_cutoff/samprate);
    [bhigh,ahigh] = butter(n_poles,HPF_cutoff/samprate,'high');
    idx_emg = contains(emg.Properties.VariableNames,'EMG');
    emg{:,idx_emg} = filtfilt(blow,alow,abs(filtfilt(bhigh,ahigh,emg{:,idx_emg})));
    cds_bin.emg = decimate_signals(emg,emg_list,bin_size);
    
    clear emg;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Process OpenSim data
% Figure out if we have opensim data
opensim_analog_idx = 0;
for i=1:length(cds.analog)
    header = cds.analog{i}.Properties.VariableNames;
    if any(contains(header,'_ang')) || any(contains(header,'_vel')) || any(contains(header,'_len')) || any(contains(header,'_muscVel'))
        opensim_analog_idx = i;
        break
    end
end
if opensim_analog_idx > 0
    opensim=cds.analog{opensim_analog_idx};
    opensimList = opensim.Properties.VariableNames;
    
    % TODO?: replace this with something that separates out the different types of
    % opensim data, e.g. separate joints and muscles, kinematics and
    % dynamics
    
    % Assign to a new 'opensim' variable
    cds_bin.opensim = decimate_signals(opensim,opensimList,bin_size);
    clear opensim;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Process Markers data
% Figure out if we have marker data
marker_analog_idx = 0;
for i=1:length(cds.analog)
    header = cds.analog{i}.Properties.VariableNames;
    if any(contains(header,'Frame')) || any(contains(header,'Marker'))
        marker_analog_idx = i;
        break
    end
end
if marker_analog_idx > 0
    markers=cds.analog{marker_analog_idx};
    markers = markers(:,2:end); % get rid of frames column
    
    % recondition table to have only single column variables
    new_markers = table(markers.t,'VariableNames',{'t'});
    marker_postfix = {'_x','_y','_z'};
    marker_names = markers.Properties.VariableNames;
    for marker_ctr = 2:width(markers)
        % for each column in each marker
        for col_ctr = 1:3
            % add new column to new_markers
            col_name = [marker_names{marker_ctr} marker_postfix{col_ctr}];
            new_col = table(markers.(marker_names{marker_ctr})(:,col_ctr),'VariableNames',{col_name});
            new_markers = [new_markers new_col];
        end
    end
    markersList = new_markers.Properties.VariableNames;
    clear marker_analog_idx header markers marker_postfix marker_names marker_ctr col_ctr col_name new_col
    
    % Assign to a new 'markers' variable
    cds_bin.markers = resample_signals(new_markers,markersList,bin_size,cds.kin.t(1),cds.kin.t(end));
    clear new_markers
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check the time vectors for these signals to make a master one
fn = fieldnames(cds_bin);
cds_bin.t = roundTime(cds_bin.(fn{1}).t);
% check the time vector lengths just to be sure
for i = 2:length(fn)
    if ~isempty(cds_bin.(fn{i})) && (size(cds_bin.(fn{i}).t,1) ~= size(cds_bin.t,1))
        error('Time is different!');
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Neural data now
for array = 1:length(arrays)
    [bs,sg] = bin_spikes(cds.units,unit_idx{array},cds_bin.t);
    cds_bin.([arrays{array} '_spikes']) = bs; clear bs;
    cds_bin.([arrays{array} '_unit_guide']) = sg; clear sg;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% convert events in trial table into bin indices
cds_bin.trials = bin_events(cds.trials,intersect(event_list,time_events),cds_bin.t);
param_events = setdiff(event_list,time_events);
% add in the non-time events
for var = 1:length(param_events)
    if ismember(param_events{var},cds.trials.Properties.VariableNames)
        cds_bin.trials.(param_events{var}) = cds.trials.(param_events{var});
    else
        warning(['Requested event ' param_events{var} ' not found in CDS.']);
    end
end
% get a new master event list in case any weren't found
event_list = fieldnames(cds_bin.trials);

extra_time_temp = extra_time;
% This is a little "hack" in case all data is desired
if all_points % here we want to include all data
    cds_bin.trials.endTime = [cds_bin.trials.startTime(2:end)-1; length(cds_bin.t)];
    extra_time = [0 0];
else
    extra_time = round(extra_time/bin_size); % convert to number of bins
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute the timestamps for raster plots (if wanted)
if include_ts
    for i = 1:length(idx_trials)
        trialStart(i) = cds_bin.t(cds_bin.trials.startTime(idx_trials(i)));
        trialEnd(i) = cds_bin.t(cds_bin.trials.endTime(idx_trials(i)));
        for j = 1:length(arrays)
            for k = 1:length(unit_idx{j})
                timestamps{i,j,k} = cds.units(unit_idx{j}(k)).spikes.ts(cds.units(unit_idx{j}(k)).spikes.ts> trialStart(i)-extra_time_temp(1) & cds.units(unit_idx{j}(k)).spikes.ts < trialEnd(i)+extra_time_temp(2)) - trialStart(i);
            end
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Loop along trials and break into struct entries
for i = 1:length(idx_trials)
    iTrial = idx_trials(i);
    
    % find trial start/end times
        t_start = cds_bin.trials.startTime(iTrial) - extra_time(1);
        t_end = cds_bin.trials.endTime(iTrial) + extra_time(2);

    idx = t_start:t_end-1;
    % check if any trials have idx<0, i.e. first trial starts within
    % extra_time samples of beginning of file
    if any(idx<0) || any(idx>length(cds_bin.kin.t))
        % skip trial
        warning(['Trial ' num2str(iTrial) ' extends outside of file, skipping'])
        continue
    end
    
    % adjust start time for use next
    cds_bin.trials.startTime(iTrial) = cds_bin.trials.startTime(iTrial)+1;
    
    if ~isempty(idx)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Add trial index markers
        for e = 1:length(event_list)
            temp = cds_bin.trials.(event_list{e});
            temp = temp(iTrial,:);
            
            % check to see if there's an alias and use it
            temp_name = event_list{e};
            if ~isempty(event_alias)
                alias_idx = find(strcmpi(event_alias(:,1),event_list{e}));
                if ~isempty(alias_idx)
                    temp_name = event_alias{alias_idx,2};
                end
            end
            
            if ismember(event_list{e},time_events) % adjust to be relative to first bin
                trial_data(i).(['idx_' temp_name]) = temp - idx(1);
            else % take parameter value
                trial_data(i).(temp_name) = temp;
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % ADD kinematics
        if isfield(cds_bin,'kin') && ~isempty(cds_bin.kin)
            trial_data(i).pos = [cds_bin.kin.x(idx)-pos_offset(1),cds_bin.kin.y(idx)-pos_offset(2)];
            trial_data(i).vel = [cds_bin.kin.vx(idx),cds_bin.kin.vy(idx)];
            trial_data(i).acc = [cds_bin.kin.ax(idx),cds_bin.kin.ay(idx)];
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Add force
        if isfield(cds_bin,'force') && ~isempty(cds_bin.force)
            trial_data(i).force = [cds_bin.force.fx(idx),cds_bin.force.fy(idx)];
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Add emg
        if isfield(cds_bin,'emg') && ~isempty(cds_bin.emg)
            fn = cds.emg.Properties.VariableNames;
            fn = fn(~strcmpi(fn,'t'));
            
            trial_data(i).emg = zeros(length(idx),length(fn));
            % loop along the muscles to decimate
            for muscle = 1:length(fn)
                temp = cds_bin.emg.(fn{muscle});
                trial_data(i).emg(:,muscle) = temp(idx);
            end
            % add emg names
            trial_data(i).emg_names = fn;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Add opensim
        if isfield(cds_bin,'opensim') && ~isempty(cds_bin.opensim)
            fn = opensimList;
            fn = fn(~strcmpi(fn,'t'));
            
            trial_data(i).opensim = zeros(length(idx),length(fn));
            % loop along the opensim to decimate
            for entry = 1:length(fn)
                temp = cds_bin.opensim.(fn{entry});
                trial_data(i).opensim(:,entry) = temp(idx);
            end
            % add opensim names
            trial_data(i).opensim_names = fn;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Add markers
        if isfield(cds_bin,'markers') && ~isempty(cds_bin.markers)
            fn = markersList;
            fn = fn(~strcmpi(fn,'t'));
            
            trial_data(i).markers = zeros(length(idx),length(fn));
            % loop along the markers to decimate
            for entry = 1:length(fn)
                temp = cds_bin.markers.(fn{entry});
                trial_data(i).markers(:,entry) = temp(idx);
            end
            % add opensim names
            trial_data(i).marker_names = fn;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Add array data
        for array = 1:length(arrays)
            use_array_name = arrays{array};
            if ~isempty(array_alias)
                temp_idx = ismember(array_alias(:,1),arrays{array});
                if sum(temp_idx) == 0
                    warning([arrays{array} ': not found in alias list. Using original name instead.']);
                else
                    use_array_name = array_alias{temp_idx,2};
                end
            end
            binned_spikes = cds_bin.([arrays{array} '_spikes']);
            trial_data(i).([use_array_name '_spikes']) = binned_spikes(idx,:);
            trial_data(i).([use_array_name '_unit_guide']) = cds_bin.([arrays{array} '_unit_guide']);
            if include_ts
                trial_data(i).([use_array_name '_ts']) = squeeze(timestamps(i, array, 1:length(unit_idx{array})));
                trial_data(i).trial_start_time = trialStart(i);
            end
            if include_naming
                trial_data(i).([use_array_name '_naming']) = conversion{array};
            end
        end
    end
end

%%%%% Check for bad trials that got skipped
fn = getTDfields(trial_data,'idx');
trial_data = trial_data(~cellfun(@isempty,{trial_data.(fn{1})}));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Package up parameter output
td_params = struct( ...
    'event_list',{event_list}, ...
    'array_alias',{array_alias}, ...
    'event_alias',{event_alias}, ...
    'trial_results',{trial_results}, ...
    'exclude_units',exclude_units, ...
    'bin_size',bin_size, ...
    'extra_time',extra_time, ...
    'all_points',all_points);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BUNCHA SUB FUNCS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = decimate_signals(data,var_list,bin_size)
% decimates continuous data from CDS
% assumes data is cds field, e.g. cds.kin
if ~isempty(data)
    out = struct();
    
    dt = round(mode(diff(data.t)),7);
    for var = 1:length(var_list)
        out.(var_list{var}) = decimate(data.(var_list{var}),round(bin_size/dt));
    end
else
    out = [];
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = resample_signals(data,var_list,bin_size,start_time,end_time)
% resamples continuous data from CDS that may be irregularly sampled
% assumes data is cds field, e.g. cds.markers
% Uses spline to reconstruct missing points
% alters time vector to fit given start_time, end_time, and bin_size
% also assumes each variable in data table is single column
if ~isempty(data)
    out = struct();
    
    for var = 1:length(var_list)
        x = data.(var_list{var});
        tx = data.t;
        
        if size(x,2)>1
            error('Data in table must have single column variables to resample correctly')
        end
        % set edges of signal to 0 to remove edge effects from spline fit
        % (see https://www.mathworks.com/help/signal/examples/resampling-nonuniformly-sampled-signals.html)
        x(1) = x(find(~isnan(x(:,1)),1,'first'));
        x(end) = x(find(~isnan(x(:,1)),1,'last'));
        a(1) = (x(end)-x(1)) / (tx(end)-tx(1));
        a(2) = x(1);
        xdetrend = x - polyval(a,tx);

        % Make sure time vector endpoints are start_time and end_time
        if start_time>tx(1)
            xdetrend = xdetrend(tx>start_time);
            tx = tx(tx>start_time);
        end
        tx = [start_time; tx];
        xdetrend = [xdetrend(1,:); xdetrend];
        if end_time<tx(end)
            xdetrend = xdetrend(tx<end_time);
            tx = tx(tx<end_time);
        end
        tx = [tx;end_time];
        xdetrend = [xdetrend; xdetrend(end,:)];
        
        % resample signal
        [ydetrend,ty] = resample(xdetrend,tx,1/bin_size,'spline');
        
        % check time vector (sometimes ty is one sample too long?)
        tol = 1e-6;
        if abs(ty(end)-end_time)>tol
            % check what's wrong
            if ty(end)>end_time
                while ty(end)-end_time>tol
                    % probably an extra sample, remove it
                    ty=ty(1:end-1);
                    ydetrend = ydetrend(1:end-1);
                end
            else
                warning('Something screwy going on with the end of ty in resample_signals...')
            end
        end
        if abs(ty(1)-start_time)>tol
            % check what's wrong
            if ty(1)<start_time
                while start_time-ty(1)>tol
                    % probably an extra sample, remove it
                    ty=ty(2:end);
                    ydetrend = ydetrend(2:end);
                end
            else
                warning('Something screwy going on with the start of ty in resample_signals...')
            end
        end
        
        % add back trend line
        out.(var_list{var}) = ydetrend + polyval(a,ty);
    end
    out.t = ty;
else
    out = [];
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = bin_events(trials,event_list,t_bin)
% bins events from CDS
out = struct();
for e = 1:length(event_list)
    all_events = trials.(event_list{e});
    nan_bins = NaN(size(all_events));
    for i = 1:size(all_events,2) % e.g. in RW and TRT there can be multiple go cues
        temp = find(histcounts(all_events(:,i),t_bin));
        nan_bins(~isnan(all_events(:,i)),i) = temp;
    end
    
    out.(event_list{e}) = nan_bins;
end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [binned_spikes,sg] = bin_spikes(units,unit_idx,t_bin)

% note that histcounts outputs rows
binned_spikes = zeros(size(unit_idx,1),length(t_bin));

% histcounts works weirdly and needs an extra bin
t_bin = [t_bin; t_bin(end)+mode(diff(t_bin))];

sg = zeros(length(unit_idx),2);
for unit = 1:length(unit_idx)
    % get the spikes for that cell in the current time window
    ts = units(unit_idx(unit)).spikes.ts;
    binned_spikes(unit,:) = histcounts(ts,t_bin);
    sg(unit,:) = [units(unit_idx(unit)).chan, units(unit_idx(unit)).ID];
end
% must transform to have same dimensions as kinematics etc
binned_spikes = binned_spikes';
end

