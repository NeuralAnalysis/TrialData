%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function trial_data = truncateAndBin(trial_data, varargin)
%
%   Will truncate all of the time-signals of each trial_data trial and also
% compute new bins if desired. Can do either one of those alone by simply
% passing in the inputs required for each. Note that the truncation happens
% BEFORE the re-binning.
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
%   trial_data = truncateAndBin(trial_data, 5, {'idx_target_on',0}, {'idx_trial_end',-3});
%
% Written by Matt Perich. Updated Feb 2017.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function trial_data = truncateAndBin(trial_data,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin == 2 % only a bin size was passed. Just rebin.
    do_trunc = false;
    num_bins = varargin{1};
elseif nargin == 3 % only the start/end indices were passed. Just truncate.
    do_trunc = true;
    idx_start = varargin{1};
    idx_end = varargin{2};
    num_bins = 1;
elseif nargin == 4
    do_trunc = true;
    temp = find(cellfun(@(x) iscell(x),varargin));
    if length(temp) ~= 2, error('Need to give two and only two alignment inputs'); end
    idx_start = varargin{temp(1)};
    idx_end = varargin{temp(2)};
    
    num_bins = varargin{cellfun(@(x) ~iscell(x),varargin)};
else
    error('Inputs not recognized.');
end

fn_spikes = getTDfields(trial_data,'spikes');
fn_kin = getTDfields(trial_data,'cont');
fn_idx = getTDfields(trial_data,'idx');

%%%%
% THIS IS A HACK FOR NOW BECAUSE MATT IS TOO LAZY TO REMAKE HIS SAVED
% STRUCTS FROM SCRATCH. CAN BE DELETED ONCE IT GOES PUBLIC.
if ~isfield(trial_data,'bin_size')
    disp('Bin size not found in struct. Defaulting to 10 ms.');
    [trial_data.bin_size] = deal(0.01);
end
%%%%

for trial = 1:length(trial_data)
    % assumes there will always be a pos
    t = 1:size(trial_data(trial).pos,1);
    
    if do_trunc
        t_start = floor(trial_data(trial).(idx_start{1}) + idx_start{2});
        t_end = ceil(trial_data(trial).(idx_end{1}) + idx_end{2});
        if t_end > t(end)
            warning('Requested end time went beyond trial time...')
            t_end = length(t);
        end
    else
        t_start = 1;
        t_end = length(t);
    end
    
    t_bin = t(t_start):num_bins:t(t_end);
    
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
    
    % process kinematic fields
    for iKin = 1:length(fn_kin)
        temp = trial_data(trial).(fn_kin{iKin});
        % fr is size bins x neurons
        kin = zeros(length(t_bin)-1,size(temp,2));
        for iBin = 1:length(t_bin)-1
            kin(iBin,:) = mean(temp(t_bin(iBin):t_bin(iBin+1)-1,:),1);
        end
        trial_data(trial).(fn_kin{iKin}) = kin;
    end
    
    % process idx fields
    for iIdx = 1:length(fn_idx)
        temp = t(trial_data(trial).(fn_idx{iIdx}));
        if temp < t_bin(1) || temp > t_bin(end)
            trial_data(trial).(fn_idx{iIdx}) = [];
        else
            trial_data(trial).(fn_idx{iIdx}) = find(t_bin <= temp,1,'last');
        end
    end
end



