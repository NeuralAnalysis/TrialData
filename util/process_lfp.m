%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [lfp_data,t_fft,freq_bands] = process_lfp(data,t,params)
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
freq_bands = [ ...% first column is low cutoff, second is high
    0 4; ...
    4 8; ...
    8 12; ...
    12 18; ...
    18 25; ...
    25 50; ...
    50 80; ...
    80 150; ...
    150 300];% in Hz
samprate          =  [];       % sampling rate of data
lfp_method        =  'fft';    % 'fft' for now
downsample_fac    =  0;        % factor for decimate
remove_common_avg =  true;     % whether to do common average referencing
remove_time_avg   =  true;     % whether to center  data in time
filt_line_noise   =  true;     % whether to filter line noise
line_noise_freq   =  50;       % 60 for US, 50 for Europe
filt_poles        =  2;        % for butterworth filter
filt_high_freq    =  true;     % apply LPF at highest freq_band
fft_window_size   =  2048;     % in samples
fft_step          =  0.01;     % in seconds (should be bin size)
fft_win_fun       =  @hamming; % windowing function handle
bandwidth         =  50; % max frequency. Add  case to default to  fs/2?
do_LMP            =  false;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isempty(params), assignParams(who,params); end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% In case LMP is wanted, add a new row to the frequency bands to save the
% LMP data, this will have 0 at both frequency limits
if do_LMP
    freq_bands = cat(1,[0 0], freq_bands);
end

if size(freq_bands,2) ~= 2
    error('LFP frequency range input should be two columns: [LOW, HIGH]');
end

data = double(data);


% step for FFT sliding window
step = round(fft_step*samprate);


% downsample
if downsample_fac > 0
    data_ds = zeros(ceil(size(data,1)/downsample_fac),size(data,2));
    for iSig = 1:size(data,2)
        data_ds(:,iSig) = decimate(data(:,iSig),downsample_fac);
    end
    samprate = samprate/downsample_fac;
    data = data_ds;
    t = decimate(t,downsample_fac);
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
lfp_data = nan(N,size(freq_bands,1)*size(data,2));

% chop the data into blocks, to avoid memory overflows!
block_size = fft_step*samprate*100000; % seemed reasonable to take 100000 time bins blocks; could be turn into a param
nblocks = ceil(N/block_size);



for iSig = 1:size(data,2)
    disp(['Processing LFP from channel ' num2str(iSig)]);
    tic;
    
    for b = 1:nblocks
        
        % carefully select the start of each window, so when we stitch the
        % data, we don't have discontinuities or NaNs --this implies that
        % the FFT winfows overlap by a number of points that depends on the
        % FFT window size
        bstart = single( b + block_size*(b-1) - 2*floor(fft_window_size/2)*(b-1) );
        
        if b ~= nblocks
            bend = bstart + block_size - 1;
        else
            bend = N;
        end
        
        switch lower(lfp_method)
            case 'fft'
                
                % get FFT over time
                [data_fft, freq, ~] = trFFT(data(bstart:bend,iSig), ...
                    fft_window_size, ...
                    step, ...
                    samprate, ...
                    win_func);
                

            otherwise
                error('LFP method not  recognized');
        end
        
        if do_LMP
            % LMP is the moving average of the signal using 50ms window
            data_lmp = movmean(data(bstart:bend,iSig),samprate*0.05);
        end

        data_fft = abs(data_fft);

        % find average power in each  band for this channel
        for iBand = 1:size(freq_bands,1)
            if freq_bands(iBand,1) == 0 && freq_bands(iBand,2) == 0
                temp = data_lmp;
            else
                idx_freq = freq >= freq_bands(iBand,1)  &  freq  <= freq_bands(iBand,2);
                temp = mean(data_fft(idx_freq,:),1)';
                % lfp_data(:,size(data,2)*(iBand-1)+iSig) = temp;
            end
            
            % find idx for ~NaN bins (trFFT introduces NaNs at the
            % beginning and end of data_fft because the FFT is computed in windows ofc)
            temp_data = temp( floor(fft_window_size/2)+1 : end-floor(fft_window_size/2)+1 );
            
            % find idx for where the data belong in the lfp_data matrix
            % --remember, this is to compensate for the NaN padding at the
            % beginning and end of data_fft
            idx_start = single( bstart + floor(fft_window_size/2) );
            % idx_end = single( bstart + block_size - floor(fft_window_size/2) );
            idx_end = idx_start + length(temp_data) - 1;
            lfp_data(idx_start:idx_end,size(freq_bands,1)*(iSig-1)+iBand) = temp_data;
        end

    end

    toc;
end


t_fft = 0:step/samprate:(size(lfp_data,1)-1)/samprate;

% new time vector assumes it starts at 0, so subtract (or add)
%   any offset identified in the original  processing
t_fft = t_fft + t(1);


end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [data_fft,freq,window_centers]=trFFT( ...
    data, ...
    window_size, ...
    step, ...
    samplerate, ...
    win_func )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function trFFT - time-resolved fourier-transformation
%
% parameter:
% data        - signal; first dimension: time; second dimension: trials
% windowSize  - length of time window
% step        - time steps in which time window is moved
% sampleRate  - sampling rate
% win_func : window function (typically  a normalized hamming window)
%
% return values:
% FFTofData   - fourie transformed signal, first dimension: frequency; second dimension: time; third dimension: trials
% frequencies - vector containing frequencies
% winCenter   - vector containing times of centers of windows
%

% Tomislav Milekovic, 06/19/2008

freq = samplerate/2*linspace(0,1,ceil(window_size/2+1));
window_centers = window_size:step:size(data,1);


data_fft = nan([length(freq) length(window_centers) size(data,2)]);



for iSig=1:size(data,2)
    for iWin=1:length(window_centers)
        
        temp_data=data(window_centers(iWin)-window_size+1:window_centers(iWin),iSig);
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


window_centers=window_centers-floor(window_size/2);

window_centers = cat(2, ...
    0:step:floor(window_size/2)-1, ...
    window_centers,  ...
    window_centers(end)+step:step:window_centers(end)+floor(window_size/2)-1);

end
