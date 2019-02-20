%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dataResampled,ty] = resample_signals(data,tx,params)
% Resample signals by upsampling and downsampling (with requisite filters)
%   some extra code in here to make sure we're not extrapolating, i.e. the
%   signals start and end at roughly the same time point.

samprate        =  [];
bin_size        =  [];
if ~isempty(params), assignParams(who,params); end

% get integers for resampling ratio
[P,Q] = rat((1/bin_size)/samprate,1e-7);

% need to detrend first...
% detrend first because resample assumes endpoints are 0
a = zeros(2,size(data,2));
dataDetrend = zeros(size(data,1),size(data,2));
for i = 1:size(data,2)
    % in case start or end are nans
    nanners = isnan(data(:,i));
    data_poly = data(~nanners,i);
    t_poly = tx(~nanners);
    
    a(1,i) = (data_poly(end)-data_poly(1))/(t_poly(end)-t_poly(1));
    a(2,i) = data_poly(1);

    dataDetrend(:,i) = data(:,i)-polyval(a(:,i),tx);
end

% resample
temp = resample(dataDetrend,P,Q);

% interpolate time vector
% using upsample -> downsample to save memory (it's the same thing
% as the reverse) but it adds extra points at the end that aren't
% in the resampled data
resamp_vec = ones(size(data,1),1);
resamp_vec = upsample(downsample(resamp_vec,Q),P);
ty=upsample(downsample(tx,Q),P);
ty=interp1(find(resamp_vec>0),ty(resamp_vec>0),(1:length(ty))');

% get rid of extrapolated points at the end
extrap_idx = isnan(ty);
ty(extrap_idx) = [];
temp(extrap_idx(1:size(temp,1)),:) = [];

% retrend...
dataResampled = zeros(size(temp,1),size(temp,2));
for i=1:size(dataDetrend,2)
    dataResampled(:,i) = temp(:,i)+polyval(a(:,i),ty(:,1));
end

end
