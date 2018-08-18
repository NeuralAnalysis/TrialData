%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function trial_data = binTD(trial_data, num_bins)
%
%   Will adjust bin size of trial_data struct.
%
%   Can provide num_bins as a char input for some special functionality.
% Currently it's only implemented to do 'average', which  will return a
% single bin representing the average
%
% INPUTS:
%   trial_data : (struct) trial_data struct
%   num_bins   : (int) how many bins to group together
%                       OR
%                (char) 'average': reduces to a single bin by averaging for
%                                   continuous signals and counting for
%                                   discrete signals (e.g. spikes). In this
%                                   case bin_size now contains the total
%                                   time for each trial. Adds a flag too!
%
% Note bin number for alignment can be negative to go before idx
% Also note that start is assumed to always come before end
%
% EXAMPLE:
%   trial_data = binTD(trial_data, 5);
%
% Written by Matt Perich. Updated August 2018.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function trial_data = binTD(trial_data,num_bins)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isstruct(trial_data), error('First input must be trial_data struct!'); end
if nargin < 2, num_bins = 1; end
if ischar(num_bins)
    switch lower(num_bins)
        case 'average'
            disp('Returning single bin representing average (or sum for spikes)');
            do_avg =  true;
    end
else
    do_avg = false;
end

fn_spikes = getTDfields(trial_data,'spikes');
fn_time = getTDfields(trial_data,'time');
fn_idx = getTDfields(trial_data,'idx');

if ~isfield(trial_data,'is_continuous') || ~any([trial_data.is_continuous])
    for trial = 1:length(trial_data)
        % get the time vectors for this trial
        t = 1:size(trial_data(trial).(fn_time{1}),1);
        if do_avg % make it a single big bin
            t_bin = [1 t(end)+1];
            num_bins = t(end); % this will tell you how  many
        else
            t_bin = 1:num_bins:t(end)+1;
        end
        
        % update entry to new bin size
        trial_data(trial).bin_size = num_bins * trial_data(trial).bin_size;
        
        % process spike fields
        for iArray = 1:length(fn_spikes)
            temp = trial_data(trial).(fn_spikes{iArray});
            % fr is size bins x neurons
            fr = zeros(length(t_bin)-1,size(temp,2));
            for iBin = 1:length(t_bin)-1
                fr(iBin,:) = sum(temp(t_bin(iBin):t_bin(iBin+1)-1,:),1);
            end
            trial_data(trial).(fn_spikes{iArray}) = fr;
        end
        
        % process other time fields
        for iSig = 1:length(fn_time)
            % make sure it's not spikes because I already did that
            if isempty(strfind(fn_time{iSig},'_spikes'))
                temp = trial_data(trial).(fn_time{iSig});
                % fr is size bins x neurons
                kin = zeros(length(t_bin)-1,size(temp,2));
                for iBin = 1:length(t_bin)-1
                    kin(iBin,:) = mean(temp(t_bin(iBin):t_bin(iBin+1)-1,:),1);
                end
                trial_data(trial).(fn_time{iSig}) = kin;
            end
        end
        
        % process idx fields
        for iIdx = 1:length(fn_idx)
            temp = trial_data(trial).(fn_idx{iIdx});
            if temp > length(t), temp = length(t); end
            if temp <= 0, temp = NaN; end
            if ~isnan(temp)
                temp(temp > max(t)) = [];
                temp = t(temp);
                % in cases like the RW go cues, there can be multiple idx_, so
                % loop along them
                temp_adjust = zeros(size(temp));
                for i = 1:size(temp,2)
                    if isempty(temp(i)) || (temp(i) < t_bin(1) || temp(i) > t_bin(end))
                        temp_adjust(i) = NaN;
                    else
                        temp_adjust(i) = find(t_bin <= temp(i),1,'last');
                    end
                end
                trial_data(trial).(fn_idx{iIdx}) = temp_adjust;
            else
                trial_data(trial).(fn_idx{iIdx}) = NaN; % Should this be NaN or []?
            end
            
        end
    end

    % put in a flag if it's an average
    if do_avg
        for trial = 1:length(trial_data)
            trial_data(trial).is_time_averaged = true;
        end
    end
    
elseif isfield(trial_data,'is_continuous') && all([trial_data.is_continuous]) && ~do_avg % preserve continuous data
    % code here is very similar to above, but I use a global time vector
    % instead of per trial
    t = 1:size(cat(1,trial_data.(fn_time{1})),1);
    t_bin = 1:num_bins:t(end)+1;
    trial_starts = cumsum([1,cellfun(@(x) size(x,1),{trial_data.(fn_time{1})})]);
    trial_starts_bin = zeros(size(trial_starts));
    for trial = 1:length(trial_starts)
        trial_starts_bin(trial) = find(t_bin <= trial_starts(trial),1,'last');
    end
    
    % process spike fields
    for iArray = 1:length(fn_spikes)
        temp = cat(1,trial_data.(fn_spikes{iArray}));
        % fr is size bins x neurons
        fr = zeros(length(t_bin)-1,size(temp,2));
        for iBin = 1:length(t_bin)-1
            fr(iBin,:) = sum(temp(t_bin(iBin):t_bin(iBin+1)-1,:),1);
        end
        % break fr back out by trial
        for trial = 1:length(trial_data)
            idx = trial_starts_bin(trial):trial_starts_bin(trial+1)-1;
            trial_data(trial).(fn_spikes{iArray}) = fr(idx,:);
        end
    end
    
    % process other time fields
    for iSig = 1:length(fn_time)
        % make sure it's not spikes because I already did that
        if isempty(strfind(fn_time{iSig},'_spikes'))
            temp = cat(1,trial_data.(fn_time{iSig}));
            % kin is size bins x signals
            kin = zeros(length(t_bin)-1,size(temp,2));
            for iBin = 1:length(t_bin)-1
                kin(iBin,:) = mean(temp(t_bin(iBin):t_bin(iBin+1)-1,:),1);
            end
            for trial = 1:length(trial_data)
                idx = trial_starts_bin(trial):trial_starts_bin(trial+1)-1;
                trial_data(trial).(fn_time{iSig}) = kin(idx,:);
            end
        end
    end
    
    % process idx fields
    for iIdx = 1:length(fn_idx)
        
        % does not like empty fields
        for trial = 1:length(trial_data)
            if isempty(trial_data(trial).(fn_idx{iIdx}))
                trial_data(trial).(fn_idx{iIdx}) = NaN;
            end
            if length(trial_data(trial).(fn_idx{iIdx})) > 1
                error('Multiple events are not supported yet. Sorry!');
            end
        end
        events = [trial_data.(fn_idx{iIdx})];
        events = events + trial_starts(1:end-1)-1;
        events(events > length(t)) = length(t);
        
        for trial = 1:length(events)
            temp = events(trial);
            if ~isnan(temp)
                % in cases like the RW go cues, there can be multiple idx_, so
                % loop along them
                temp_adjust = zeros(size(temp));
                for i = 1:size(temp,2)
                    if isempty(temp(i)) || (temp(i) < t_bin(1) || temp(i) > t_bin(end))
                        temp_adjust(i) = NaN;
                    else
                        temp_adjust(i) = find(t_bin <= temp(i),1,'last');
                    end
                end
                trial_data(trial).(fn_idx{iIdx}) = temp_adjust-trial_starts(trial)+1;
            else
                trial_data(trial).(fn_idx{iIdx}) = NaN; % Should this be NaN or []?
            end
        end
        
    end
    
    for trial = 1:length(trial_data)
        % update meta entry to new bin size
        trial_data(trial).bin_size = num_bins * trial_data(trial).bin_size;
    end
    
    % put in a flag if it's an average
    if do_avg
        for trial = 1:length(trial_data)
            trial_data(trial).is_time_averaged = true;
        end
    end
else
    error('Some trials are from continuous data, some are not. binTD is confused.');
end






