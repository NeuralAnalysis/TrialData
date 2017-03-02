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

% kernel half length is 3�SD out
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