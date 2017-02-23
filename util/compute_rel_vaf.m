function VAF = compute_rel_vaf(Y,Yhat1,Yhat2)
    % Calculate VAF for given data
    % INPUTS: Y - Actual output variable(s) (each row corresponds to a
    %               datapoint, each column a dimension)
    %         Yhat1 - Model-predicted reference variable
    %         Yhat2 - Model-predicted output variable
    % Output: VAF - variance-accounted-for, between 0 and 1 (same as R^2)
    % WARNING: NOT TESTED FOR MULTI-DIMENSIONAL Y
    
    %calculate sum of square errors    
    SSE1 = (Y-Yhat1)'*(Y-Yhat1);
    SSE2 = (Y-Yhat2)'*(Y-Yhat2);
    
    %find VAF
    VAF = 1-SSE2/SSE1;
end