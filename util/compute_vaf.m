function VAF = compute_vaf(Y,Yhat)
    % Calculate VAF for given data
    % INPUTS: Y - Actual output variable(s) (each row corresponds to a
    %               datapoint, each column a dimension)
    %         Yhat - Model-predicted output variable(s)
    % Output: VAF - variance-accounted-for, between 0 and 1 (same as R^2)
    % WARNING: NOT TESTED FOR MULTI-DIMENSIONAL Y
    
    %calculate sum of square errors    
    SSE = (Y-Yhat)'*(Y-Yhat);
    
    %calculate mean of actual output
    meanY = mean(Y);
    
    %find sum of squared deviations from mean
    SS = (Y-repmat(meanY,length(Y),1))'*(Y-repmat(meanY,length(Y),1));
    
    %find VAF
    VAF = 1-SSE/SS;
end