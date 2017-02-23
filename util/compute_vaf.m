function VAF = compute_vaf(Y,Yhat)
    % Calculate VAF for given data
    % INPUTS: Y - Actual output variable(s) (each row corresponds to a
    %               datapoint, each column a dimension)
    %         Yhat - Model-predicted output variable(s)
    % Output: VAF - variance-accounted-for, between 0 and 1 (same as R^2)
    % 
    % for multi-dimensional data, rows are observations, columns are variables
    
    %calculate sum of square errors    
    SSE = sum((Y-Yhat).^2,1);
    
    %calculate mean of actual output
    meanY = mean(Y,1);
    
    %find sum of squared deviations from mean
    SS = sum((Y-repmat(meanY,size(Y,1),1)).^2);
    
    %find VAF
    VAF = 1-SSE./SS;
end