%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function trial_data = smoothSignal(trial_data,params)
%
% This function will smooth and/or square root transform spikes
%
% INPUTS:
%   trial_data : the struct
%   params     :
%     .signals        : field names to smooth (single string or cell array)
%                          Note: defaults to all signals
%     .do_smoothing   : flag to convolve spikes with gaussian (default: true)
%     .kernel_SD      : kernel s.d. for smoothing (default: 0.05)
%     .sqrt_transform : flag to square root transform (default: false)
%                           Mostly meant for spiking data
%     .calc_rate      : flag to calculate rate (divide by bin size)
%                           Mostly meant for spiking data (default: false)
%
% OUTPUTS:
%   trial_data : the struct with all (signals{}) fields smoothed
%
% Written by Matt Perich. Updated Feb 2017.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function trial_data = smoothSignal(trial_data,params)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEFAULT PARAMETER VALUES
signals         =  getTDfields(trial_data,'time'); % default to all signals
sqrt_transform  =  false;
do_smoothing    =  true;
kernel_SD       =  0.05;
calc_rate       =  false;
if nargin > 1, assignParams(who,params); end % overwrite defaults
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~iscell(signals), signals = {signals}; end
bin_size = trial_data(1).bin_size;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if do_smoothing || sqrt_transform % if you don't want to do either just passes back trial_data
    for trial = 1:length(trial_data)
        for i = 1:length(signals)
            data = trial_data(trial).(signals{i});
            if sqrt_transform, data = sqrt(data); end
            if calc_rate, data = data./bin_size; end
            if do_smoothing
                data = smooth_data(data,bin_size,kernel_SD);
            end
            trial_data(trial).(signals{i}) = data;
        end
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function data_smooth = smoothSpikesForPCA(data,dt,kernel_SD)
%
%   Convolves spike counts or firing rates with a Guassian kernel to smooth
% the responses. Made for spikes but works well with other signals as well.
% Used by getPCA. Based on gaussian_smoothing code from Juan's proc folder.
%
% INPUTS:
%   data      : array of data (rows: time, columns: signals)
%   dt        : size of time steps in fr in s
%   kernel_SD : gaussian kernel standard deviation
%
% OUTPUTS:
%   fr_smooth : smoothed data
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function data_smooth = smooth_data(data,dt,kernel_SD)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% now apply smoothing
% get nbr of channels and nbr of samples
[nbr_samples, nbr_chs]  = size(data);
% preallocate return matrix
data_smooth = zeros(nbr_samples,nbr_chs);

% kernel half length is 3·SD out
kernel_hl = ceil( 3 * kernel_SD / (dt) );
% create the kernel --it will have length 2*kernel_hl+1
kernel = normpdf( -kernel_hl*(dt) : ...
    dt : kernel_hl*(dt), ...
    0, kernel_SD );
% compute normalization factor --this factor depends on the number of taps
% actually used
nm = conv(kernel,ones(1,nbr_samples))';

% do the smoothing
for i = 1:nbr_chs
    aux_smoothed_FR     = conv(kernel',data(:,i)) ./ nm;
    % cut off the edges so that the result of conv is same length as the
    % original data
    data_smooth(:,i)    = aux_smoothed_FR(kernel_hl+1:end-kernel_hl);
end
end