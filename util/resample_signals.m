%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dataResampled,ty] = resample_signals(data,tx,params)
% Resample signals by upsampling and downsampling (with requisite filters)
%   some extra code in here to make sure we're not extrapolating, i.e. the
%   signals start and end at roughly the same time point.


if size(tx,1) == 1 && size(tx,2) ~= 1
    tx = tx';
end

samprate        =  [];
bin_size        =  [];
if ~isempty(params), assignParams(who,params); end

% get integers for resampling ratio
[P,Q] = rat((1/bin_size)/samprate,1e-7);

% need to detrend first...
% detrend first because resample assumes endpoints are 0
a = zeros(2,size(data,2));
[s1, s2] = size(data);
dataDetrend = zeros(s1,s2);

for i = 1:size(data,2)
    % in case start or end are nans
    nanners = isnan(data(:,i));
    data_poly = data(~nanners,i);
    t_poly = tx(~nanners);
    
    if isempty(t_poly)
        error('Everything is NaN in this data column!')
    end
    
    a(1,i) = (data_poly(end)-data_poly(1))/(t_poly(end)-t_poly(1));
    a(2,i) = data_poly(1);

    dataDetrend(:,i) = double(data(:,i))-polyval(a(:,i),tx);
end
clear data;


% resample -- do in blocks to avoid memory overflow with huge matrices
% --e.g., LFPs (we need to make this code better!)
temp = res(dataDetrend,P,Q);
clear dataDetrend;


% interpolate time vector
% using upsample -> downsample to save memory (it's the same thing
% as the reverse) but it adds extra points at the end that aren't
% in the resampled data
resamp_vec = ones(s1,1);
resamp_vec = upsample(downsample(resamp_vec,Q),P);
ty=upsample(downsample(tx,Q),P);
ty=interp1(find(resamp_vec>0),ty(resamp_vec>0),(1:length(ty))');


% get rid of extrapolated points at the end
extrap_idx = isnan(ty);
ty(extrap_idx) = [];
% use matfile if data matrix is too large
if size(temp,2) > 0 && size(temp,2) < 100
    temp(extrap_idx(1:size(temp,1)),:) = [];
else
    out = matfile('temp_out.mat','Writable',true);
    [size1,size2] = size(out,'res_data');
    
    out.res_data(extrap_idx(1:size1),:) = [];
end


% retrend, using a matfile if data matrix is too large
if size(temp,2) > 0 && size(temp,2) < 100
    dataResampled = zeros(size(temp,1),size(temp,2));
    for i=1:s2
        dataResampled(:,i) = temp(:,i)+polyval(a(:,i),ty(:,1));
    end
else
    dataResampled = zeros(size1,size2);

    for i=1:s2
        
        dataResampled(:,i) = out.res_data(:,i)+polyval(a(:,i),ty(:,1)); 
    end
    
    % delete MAT files from hard drive
    delete('temp_in.mat');
    delete('temp_out.mat');
end

end




% -------------------------------------------------------------------------
% Fcn to resample the data --implemented as a different fcn to save memory
%
function res_data = res(data,p,q)

if size(data,2) < 100
    
    res_data = resample(data,p,q);
    
else
    
    blocksize = 20; % seemed reasonable
    nblocks = ceil(size(data,2)/blocksize);
    
    [s1, s2] = size(data);
    
    % Save input and output data into MAT files to avoid memory errors. Not
    % very elegant but oh well
    save('temp_in.mat','data','-v7.3');
    clear data;
    
    in = matfile('temp_in.mat');
    out = matfile('temp_out.mat','Writable',true);
    
    sr1 = ceil(s1/q*p); % time bins after downsampling
    
    
    for b = 1:nblocks
    
        disp(['Resampling block ' num2str(b) ' of ' num2str(nblocks)]);
        % what variables (columns) to resample
        bstart = 1+(b-1)*blocksize;
        bend = min( b*blocksize, s2 );
           
        out.res_data(1:sr1,bstart:bend) = resample(in.data(:,bstart:bend),p,q);
    end
    
    res_data = [];
end

end