%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function data_mvavg = moving_average(data, bins);
%
%   Computes moving average of input signal(s)
%
% THIS IS A WORK IN PROGRESS DON'T TRUST IT!
%                                   -Matt, March 2017
%
% INPUTS:
%   data : the data to smooth (rows are time, columns are signals)
%   bins : how many time bins to use for each side of the filter
%             so, bins = 2 gives a moving average filter of total width = 5
%
% OUTPUTS:
%   data_mvavg : smoothed data
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function data_mvavg = moving_average(data,bins)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if 0
    data_mvavg = zeros(size(data));
    for i = 1:size(data,2)
        L = filter(ones(2*bins+1,1)/(2*bins+1),1,[data(:,i); zeros(bins,1)]);
        data_mvavg(:,i) = L(bins+1:end);
    end
    
else
    data = [nan(bins/2,size(data,2)); data; nan(bins/2,size(data,2))];
    data_mvavg = zeros(size(data,1)-(bins),size(data,2));
    for i = 1:size(data,2)
        for j = 1:size(data,1)-(bins)
            data_mvavg(j,i) = nanmean(data(j:j+bins,i));
        end
    end
end