%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function trial_data = binTD(trial_data, num_bins)
%
%   Will adjust bin size of trial_data struct.
%
% INPUTS:
%   trial_data : (struct) trial_data struct
%   num_bins   : (int) how many bins to group together
%   idx_start  : (cell) {'idx_to_align_start',num_bins_after}
%   idx_end    : (cell) {'idx_to_align_end',num_bins_after}
%
% Note bin number for alignment can be negative to go before idx
% Also note that start is assumed to always come before end
%
% EXAMPLE:
%   trial_data = binTD(trial_data, 5);
%
% Written by Matt Perich. Updated March 2017.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function trial_data = binTD(trial_data,num_bins)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 2, num_bins = 1; end

fn_spikes = getTDfields(trial_data,'spikes');
fn_time = getTDfields(trial_data,'time');
fn_idx = getTDfields(trial_data,'idx');

if ~isfield(trial_data,'is_continuous') || ~any([trial_data.is_continuous])
    for trial = 1:length(trial_data)
        % get the time vectors for this trial
        t = 1:size(trial_data(trial).(fn_time{1}),1);
        t_bin = 1:num_bins:t(end);
        
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
            if temp < 0, temp = NaN; end
            if ~isnan(temp)
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
    
elseif isfield(trial_data,'is_continuous') && all([trial_data.is_continuous]) % preserve continuous data
    % code here is very similar to above, but I use a global time vector
    % instead of per trial
    t = 1:size(cat(1,trial_data.(fn_time{1})),1);
    t_bin = 1:num_bins:t(end);
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
    
    
else
    error('Some trials are from continuous data, some are not. binTD is confused.');
end


