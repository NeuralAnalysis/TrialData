function fr_smooth = smoothSpikesForPCA(fr,bin_size,kernel_SD)

% now apply smoothing
% get nbr of channels and nbr of samples
[nbr_samples, nbr_chs]  = size(fr);
% preallocate return matrix
fr_smooth = zeros(nbr_samples,nbr_chs);

% kernel half length is 3·SD out
kernel_hl = ceil( 3 * kernel_SD / (bin_size) );
% create the kernel --it will have length 2*kernel_hl+1
kernel = normpdf( -kernel_hl*(bin_size) : ...
                            bin_size : kernel_hl*(bin_size), ...
                            0, kernel_SD );
% compute normalization factor --this factor depends on the number of taps
% actually used 
nm = conv(kernel,ones(1,nbr_samples))';

% do the smoothing
for i = 1:nbr_chs
    aux_smoothed_FR     = conv(kernel,fr(:,i)) ./ nm;
    % cut off the edges so that the result of conv is same length as the
    % original data
	fr_smooth(:,i)    = aux_smoothed_FR(kernel_hl+1:end-kernel_hl);
end