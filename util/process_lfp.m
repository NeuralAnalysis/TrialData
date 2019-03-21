%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [lfp_data,t_fft] = process_lfp(data,t,params)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Process LFP - compute power of LFP channels. Works with convertDataToTD.
% Ideas for default freq_bands = [ ...
%                       0 4; ...
%                       4 8; ...
%                       8 12; ...
%                       12 18; ...
%                       18 25; ...
%                       25 50; ...
%                       50 80; ...
%                       80 150; ...
%                       150 300; ...
%                       300 500; ...
%                       500 1200; ...
%                       1200 2000];% in Hz
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
freq_bands        =  [];       % first column is low cutoff, second is high
samprate          =  [];       % sampling rate of data
lfp_method        =  'fft';    % 'trfft' or 'taper' for now
downsample_fac    =  0;        % factor for decimate
remove_common_avg =  true;     % whether to do common average referencing
remove_time_avg   =  true;     % whether to center  data in time
filt_line_noise   =  true;     % whether to filter line noise
line_noise_freq   =  50;       % 60 for US, 50 for Europe
filt_poles        =  2;        % for butterworth filter
filt_high_freq    =  true;     % apply LPF at highest freq_band
fft_window_size   =  1024;     % in samples
fft_step          =  0.01;     % in seconds (should be bin size)
fft_win_fun       =  @hamming; % windowing function handle
bandwidth         =  50; % max frequency. Add  case to default to  fs/2?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isempty(params), assignParams(who,params); end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if size(freq_bands,2) ~= 2
    error('LFP frequency range input should be two columns: [LOW, HIGH]');
end

data = double(data);

% downsample
if downsample_fac > 0
    data_ds = zeros(ceil(size(data,1)/downsample_fac),size(data,2));
    for iSig = 1:size(data,2)
        data_ds(:,iSig) = decimate(data(:,iSig),downsample_fac);
    end
    samprate = samprate/downsample_fac;
    data = data_ds;
    t_ds = decimate(t,downsample_fac);
    if t(1) == 0 && t_ds(1) ~= 0
        t_ds = t_ds-t_ds(1);
    end
    t = t_ds;
end


% only run operations on non-nan samples (since there might be nans at
% beginning and end if data doesn't exist at those points)
data_idx = ~any(isnan(data),2);

% check data to see if there are any random NaNs in the middle
block_starts = find(diff([0;data_idx;0])>0);
block_ends   = find(diff([0;data_idx;0])<0);
assert(numel(block_starts)==1 && numel(block_ends)==1,'EMG data block is not continuous')


% filter the line noise
if filt_line_noise
    [B,A]             = butter(filt_poles,( line_noise_freq + [-2 2] )/(samprate/2),'stop');
    temp              = filtfilt(B,A,data(data_idx,:));
    data(data_idx,:)  = temp;
    data(~data_idx,:) = NaN;
end

% low pass filter with cutoff at max freq
if filt_high_freq
    [B,A]             = butter(filt_poles,max(max(freq_bands))/(samprate/2));
    temp              = filtfilt(B,A,data(data_idx,:));
    data(data_idx,:)  = temp;
    data(~data_idx,:) = NaN;
end

% remove common average
if remove_common_avg
    data = data - repmat(mean(data,2),1,size(data,2));
end
if remove_time_avg
    data = data -  repmat(mean(data,1),size(data,1),1);
end

% set up the window for the  FFT
win_func = fft_win_fun(fft_window_size);
win_func =  win_func./norm(win_func);


% anticipate the size of the fft output
N = floor(size(data,1)/round(fft_step*samprate));
lfp_data = zeros(N,size(freq_bands,1)*size(data,2));
for iSig = 1:size(data,2)
    
tic;
    switch lower(lfp_method)
        case 'taper'
            [~, data_fft, freq, ~]= multitaperSpectrum_univariate(data(:,iSig), ...
                samprate, ...
                bandwidth, ...
                size(data,1), ...
                remove_time_avg, ...
                remove_common_avg, ...
                []);
        case 'fft'
            % get FFT over time
            [data_fft, freq, t_fft] = trFFT(data(:,iSig), ...
                fft_window_size, ...
                round(fft_step*samprate), ...
                samprate, ...
                win_func);
        otherwise
            error('LFP method not  recognized');
    end
    
    data_fft = abs(data_fft);
    
    % find average power in each  band for this channel
    for iBand = 1:size(freq_bands,1)
        idx_freq = freq >= freq_bands(iBand,1)  &  freq  <= freq_bands(iBand,2);
        temp = mean(data_fft(idx_freq,:),1)';
        lfp_data(:,size(data,2)*(iBand-1)+iSig) = temp;
    end
    
toc;
end

t_fft = t_fft/samprate;

% new time vector assumes it starts at 0, so subtract (or add)
%   any offset identified in the original  processing
t_fft = t_fft + t(1);


end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [data_fft,freq,winCenter]=trFFT( ...
    data, ...
    window_size, ...
    step, ...
    samplerate, ...
    win_func)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function trFFT - time-resolved fourier-transformation
%
% parameter:
% data        - signal; first dimension: time; second dimension: trials
% windowSize  - length of time window
% step        - time steps in which time window is moved
% sampleRate  - sampling rate
% winFunction - window function
%
% return values:
% FFTofData   - fourie transformed signal, first dimension: frequency; second dimension: time; third dimension: trials
% frequencies - vector containing frequencies
% winCenter   - vector containing times of centers of windows
%

% Tomislav Milekovic, 06/19/2008

freq = samplerate/2*linspace(0,1,ceil(window_size/2+1));
winCenter=window_size:step:size(data,1);

data_fft=nan([length(freq) length(winCenter) size(data,2)]);


for iSig=1:size(data,2)
    for iWin=1:length(winCenter)
        
        temp_data=data(winCenter(iWin)-window_size+1:winCenter(iWin),iSig);
        temp_data=win_func.*temp_data;
        FFTrez=fft(temp_data);
        data_fft(:,iWin,iSig) = FFTrez(1:ceil(window_size/2+1));
        
    end
end

% pad with NaN
data_fft = cat(2, ...
    NaN(size(data_fft,1),floor(window_size/2),size(data_fft,3)), ...
    data_fft, ...
    NaN(size(data_fft,1),floor(window_size/2)-1,size(data_fft,3)));


winCenter=winCenter-floor(window_size/2);

winCenter = cat(2, ...
    0:step:floor(window_size/2)-1, ...
    winCenter,  ...
    winCenter(end)+step:step:winCenter(end)+floor(window_size/2)-1);


end



