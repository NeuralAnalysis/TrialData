%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = extract_ns5_spikes(file_data,params)
%%%%%%%%%%%
% some parameters
spiking_chans   = 1:96;
spike_method    = 'threshold'; % 'threshold', later: 'template'?
spike_threshold = 3;   % x RMS
n_poles         = 0;   % for high pass filter on threshold (0 means no filter)
HPF_cutoff      = 750; % in Hz
%%%%%%%%%%%
if nargin == 1
    params = struct();
end
assignParams(who,params); % overwrite parameters
%%%%%%%%%%%

% define HPF parameters
if n_poles > 0
    [bhigh,ahigh] = butter(n_poles,HPF_cutoff/file_data.samprate,'high');
end

%%%%%%%%%%%
% find spiking channels

data = file_data.data;

out = struct( ...
    'duration',file_data.duration, ...
    'samprate',file_data.samprate, ...
    'labels',zeros(length(spiking_chans),2), ...
    'data',{cell(1,length(spiking_chans))}, ...
    'wf',{{}});
% process each channel
for i = 1:size(data,1)
    switch lower(spike_method)
        case 'threshold'
            d = data(i,:);
            % HPF
            
            if n_poles > 0
                d = filtfilt(bhigh,ahigh,double(d));
            end
            
            sd = -rms(d);
            idx = find([0 diff(d < spike_threshold*sd)] > 0);
            
            % enforce refactory period
            if ~isempty(idx)
                idx([false diff(idx)] < 0.0011*file_data.samprate) = [];
            end
    end
    
    out.data{i} = double(idx)/file_data.samprate;
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
