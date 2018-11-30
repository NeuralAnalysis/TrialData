%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function binned_data = rebin_signals(data,params)
% Resample signals by upsampling and downsampling (with requisite filters)
%   some extra code in here to make sure we're not extrapolating, i.e. the
%   signals start and end at roughly the same time point.

samprate        =  [];
bin_size        =  0.01;
if ~isempty(params), assignParams(who,params); end

% get integers for resampling ratio
[P,Q] = rat((1/bin_size)/samprate,1e-7);

% figure out how much data we'll have at the end
resample_vector = ones(size(data,1),1);
resample_vector = downsample(upsample(resample_vector,P),Q);
resample_vector = interp1(find(resample_vector>0),resample_vector(resample_vector>0),(1:length(resample_vector))');

binned_data = zeros(length(resample_vector),size(data,2));
for i = 1:size(data,2)
    binned_data(:,i) = resample(data(:,i),P,Q);
end

% get rid of extrapolated points
% there may be nans left over at the end from upsampling (indicating
% extrapolation...we don't want that
extrap_idx = isnan(resample_vector);
binned_data(extrap_idx,:) = [];

end
