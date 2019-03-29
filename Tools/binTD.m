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
%                                   time for each trial.
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
trial_data  =  check_td_quality(trial_data);
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

% see if spikes are spikes or if they have been messed with
spikes_are_spikes = true(1,length(fn_spikes));
for iArray = 1:length(fn_spikes)
    spikes_are_spikes(iArray) = all(all(mod(getSig(trial_data,fn_spikes{iArray}),1) == 0));
end

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
            % if we think it's a  spike, do a sum
            if spikes_are_spikes
                fr(iBin,:) = sum(temp(t_bin(iBin):t_bin(iBin+1)-1,:),1);
            else % if we think it's not, do a mean
                fr(iBin,:) = mean(temp(t_bin(iBin):t_bin(iBin+1)-1,:),1);
            end
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


