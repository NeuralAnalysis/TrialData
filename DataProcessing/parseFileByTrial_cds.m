%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function trial_data = parseFileByTrial_cds(cds,inputArgs)
%
% INPUTS:
%   cds    : CDS object
%   params : a struct containing parameters
%     .meta         : a struct with a field for each meta parameter you want attached 
%                       to this file. This can handle any arbitrary information!
%     .excludeUnits : ID for which units to exclude (Default: [0,255])
%     .trialResults : which reward codes to use ('R','A','F','I')
%     .binSize      : default 0.01 sec
%     .extraTime    : [time before, time after] beginning and end of trial (default [0.5 0.3] sec)
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
trialResults  =  {'R'};
excludeUnits  =  [0,255];
binSize       =  0.01;
extraTime     =  [0.2, 0.2];
if nargin > 1
    eval(structvars(length(fieldnames(params)),params)); %overwrite parameters
else
    params = struct();
end
if ~isfield(params,'meta'), warning('WARNING: no meta information provided.'); end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% do some input processing
if ~iscell(trialResults), trialResults = {trialResults}; end

% see what array data is present
arrays = strsplit(cds.meta.array,', ');

% get info on neurons
unit_idx = cell(1,length(arrays));
for iArray = 1:length(arrays)
    unit_idx{iArray} = find(~ismember([cds.units.ID],excludeUnits) & strcmpi({cds.units.array},arrays{iArray}));
end

idx_trials = find(ismember(cds.trials.result,trialResults));

% loop along trials
trial_data = repmat(struct(),1,length(idx_trials));
for i = 1:length(idx_trials)
    iTrial = idx_trials(i);
    % add some meta data about the trial
    trial_data(i).monkey = cds.meta.monkey;
    trial_data(i).date = datestr(cds.meta.dateTime,'mm-dd-yyyy');
    trial_data(i).task = cds.meta.task;
    if any(abs(cds.trials.tgtDir) > 2*pi)
        trial_data(i).target_direction = minusPi2Pi(pi/180*cds.trials.tgtDir(iTrial));
    else
        trial_data(i).target_direction = minusPi2Pi(cds.trials.tgtDir(iTrial));
    end
    trial_data(i).trial_id = iTrial;
    trial_data(i).result = cds.trials.result(iTrial);
    trial_data(i).bin_size = binSize;
    
    % loop along all meta fields
    if isfield(params,'meta')
        for fn = fieldnames(params.meta)
            trial_data(i).(fn) = params.meta.(fn);
        end
    end
    
    % find trial start/end times
    t_start = cds.trials.startTime(iTrial) - extraTime(1);
    t_end = cds.trials.endTime(iTrial) + extraTime(2);
    
    % get kinematics
    if ~isempty(cds.kin)
        dt_kin = cds.kin.t(2)-cds.kin.t(1);
        t_dec = decimate(cds.kin.t(idx),round(binSize/dt_kin));
        trial_data(i).pos = [decimate(cds.kin.x(idx),round(binSize/dt_kin)) decimate(cds.kin.y(idx),round(binSize/dt_kin))];
        trial_data(i).vel = [decimate(cds.kin.vx(idx),round(binSize/dt_kin)) decimate(cds.kin.vy(idx),round(binSize/dt_kin))];
        trial_data(i).acc = [decimate(cds.kin.ax(idx),round(binSize/dt_kin)) decimate(cds.kin.ay(idx),round(binSize/dt_kin))];
    end
    
    % get force
    if ~isempty(cds.force)
        dt_force = cds.force.t(2)-cds.force.t(1);
        if ~exist('t_dec','var') % we only need one of these
            t_dec = decimate(cds.force.t(idx),round(binSize/dt_force));
        end
        idx = cds.force.t >= t_start & cds.force.t <= t_end;
        trial_data(i).force = [decimate(cds.force.fx(idx),round(binSize/dt_force)) decimate(cds.force.fy(idx),round(binSize/dt_force))];
        if size(trial_data(i).force,1) ~= size(trial_data(i).pos,1)
            error('Decimated force size does not match kinematics.');
        end
    end
    
    % get EMG
    if ~isempty(cds.emg)
        dt_emg = cds.emg.t(2)-cds.emg.t(1);
        if ~exist('t_dec','var') % we only need one of these
            t_dec = decimate(cds.emg.t(idx),round(binSize/dt_emg));
        end
        idx = cds.emg.t >= t_start & cds.emg.t <= t_end;
        
        fn = cds.emg.Properties.VariableNames;
        trial_data(i).emg = decimate(table2array(cds.emg(idx,~strcmpi(fn,'t'))),round(binSize/dt_emg));
        trial_data(i).emg_names = fn(~strcmpi(fn,'t'));
        
         if size(trial_data(i).emg,1) ~= size(trial_data(i).pos,1)
            error('Decimated force size does not match kinematics.');
        end
    end
    
    % get time vector for binned spikes and events
    t_bins = [t_dec', t_dec(end)+binSize];
    
    % put trial markers (target on etc) in bins for each spikes
    trial_data(i).idx_trial_start = find(histcounts(cds.trials.startTime(iTrial),t_bins));
    trial_data(i).idx_target_on = find(histcounts(cds.trials.tgtOnTime(iTrial),t_bins));
    trial_data(i).idx_go_cue = find(histcounts(cds.trials.goCueTime(iTrial),t_bins));
    trial_data(i).idx_trial_end = find(histcounts(cds.trials.endTime(iTrial),t_bins));
    
    for iArray = 1:length(arrays)
        binned_spikes = zeros(size(unit_idx{iArray},1),length(t_bins)-1);
        sg = zeros(length(unit_idx{iArray}),2);
        for unit = 1:length(unit_idx{iArray})
            % get the spikes for that cell in the current time window
            ts = cds.units(unit_idx{iArray}(unit)).spikes.ts;
            ts = ts(ts >= t_start & ts <= t_end);
            
            binned_spikes(unit,:) = histcounts(ts,t_bins);
            sg(unit,:) = [cds.units(unit_idx{iArray}(unit)).chan, cds.units(unit_idx{iArray}(unit)).ID];
        end
        
        % check to make sure all is well
        if size(binned_spikes,2) ~= size(trial_data(i).pos,1)
            error('Binned spike size does not match decimated kinematics size')
        end
        %   transpose binned_spikes to be consistent with kin
        trial_data(i).([arrays{iArray} '_spikes']) = binned_spikes';
        trial_data(i).([arrays{iArray} '_unit_guide']) = sg;
    end
    clear binned_spikes;
end

end
