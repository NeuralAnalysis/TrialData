function [orig,pred,c,b] = linPredSignal(td,in_var,out_var,num_dims,varargin)
% linear predictions of a signal
%
% WORK IN PROGRESS

fit_metric = 'vaf';

if nargin > 4
    b = varargin{1};
    new_model = false;
else
    new_model = true;
end

in_signal = cat(1,td.(in_var));
in_signal = [ones(size(in_signal,1),1), in_signal(:,1:num_dims)];
orig = cat(1,td.(out_var));

if new_model, zeros(size(in_signal,2),size(orig,2)); end
c = zeros(1,size(orig,2));
pred = zeros(size(orig));
for iDim = 1:size(orig,2)
    out_signal = orig(:,iDim);
    
    if new_model
        b(:,iDim) = in_signal\out_signal;
    end
    pred(:,iDim) = in_signal*b(:,iDim);
    
    switch lower(fit_metric)
        case 'corr'
            temp = corrcoef(out_signal,in_signal*b(:,iDim)); c(iDim) = temp(1,2);
        case 'r2'
            c(iDim) = CalculateR2(out_signal,in_signal*b(:,iDim));
        case 'vaf'
            c(iDim) = calcVAF(out_signal,in_signal*b(:,iDim));
    end
end

end