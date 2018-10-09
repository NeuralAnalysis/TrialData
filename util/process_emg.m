%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function binned_emg = process_emg(data,params)
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
bin_size        =  0.01;
if ~isempty(params), assignParams(who,params); end


[blow,alow] = butter(emg_n_poles,emg_LPF_cutoff/(samprate/2));
if length(emg_HPF_cutoff) > 1
    [bhigh,ahigh] = butter(emg_n_poles,emg_HPF_cutoff/(samprate/2));
else
    [bhigh,ahigh] = butter(emg_n_poles,emg_HPF_cutoff/(samprate/2),'high');
end

% !!! note the rectification step in the following command:
data = filtfilt(bhigh,ahigh,double(data));
data = 2*data.*data;
data = filtfilt(blow,alow,data);
data = abs(sqrt(data));

binned_emg = zeros(ceil(size(data,1)/round(bin_size*samprate)),size(data,2));
for i = 1:size(data,2)
    binned_emg(:,i) = decimate(data(:,i),round(bin_size*samprate))';
end
clear temp_data;
end
