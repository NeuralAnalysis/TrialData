function VAF = compute_vaf(Y,Yhat,Yhat2)
% Calculate VAF for given data
% INPUTS: Y - Actual output variable(s) (each row corresponds to a
%               datapoint, each column a dimension)
%         Yhat - Model-predicted output variable(s)
%         Yhat2 - A second model predicted output variable if you want to
%         do a relative metric. Yhat should be the basic model and Yhat2
%         should be the full model
% Output: VAF - variance-accounted-for, between 0 and 1 (same as R^2)
%
% for multi-dimensional data, rows are observations, columns are variables


if nargin == 2 % normal VAF
    %calculate sum of square errors
    SSE = sum((Y-Yhat).^2,1);
    
    %calculate mean of actual output
    meanY = mean(Y,1);
    
    %find sum of squared deviations from mean
    SS = sum((Y-repmat(meanY,size(Y,1),1)).^2);
elseif nargin == 3 %relative VAF
    %calculate sum of square errors
    SSE = sum((Y-Yhat2).^2,1);
    SS = sum((Y-Yhat).^2,1);
end

%find VAF
VAF = 1-SSE./SS;
