%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function emg_data = process_emg(data,params)
% Process EMG - computes envelope of emg signal and then bins it. Written
% to work with convertDataToTD.
%
% (high pass, rectify, low pass)
%   default: high pass at 10 Hz, rectify, low pass at 50 Hz
% filter
emg_LPF_cutoff  =  50;    % for EMG butterworth filter
emg_HPF_cutoff  =  [10 900];    % for EMG butterworth filter
emg_n_poles     =  4;     % for EMG butterworth filter
samprate        =  [];
if ~isempty(params), assignParams(who,params); end


[blow,alow] = butter(emg_n_poles,emg_LPF_cutoff/(samprate/2));
if length(emg_HPF_cutoff) > 1
    [bhigh,ahigh] = butter(emg_n_poles,emg_HPF_cutoff/(samprate/2));
else
    [bhigh,ahigh] = butter(emg_n_poles,emg_HPF_cutoff/(samprate/2),'high');
end

% !!! note the rectification step in the following command:
% only run operations on non-nan samples (since there might be nans at
% beginning and end if data doesn't exist at those points)
data_idx = ~any(isnan(data),2);

% check data to see if there are any random NaNs in the middle
block_starts = find(diff([0;data_idx;0])>0);
block_ends = find(diff([0;data_idx;0])<0);
assert(numel(block_starts)==1 && numel(block_ends)==1,'EMG data block is not continuous')

temp = filtfilt(bhigh,ahigh,double(data(data_idx,:)));
temp = 2*temp.*temp;
temp = filtfilt(blow,alow,temp);
data(data_idx,:) = abs(sqrt(temp));
data(~data_idx,:) = NaN;

emg_data = data;

end
