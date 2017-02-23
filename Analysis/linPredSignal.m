%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function 
% 	Compute linear predictions of a signal. For example. predict X/Y
% position using M1 firing rates. Doesn't incorporate history outright. If
% you want history you should use dupeAndShift or convBasisFunc.
%
% WORK IN PROGRESS. To do:
%   1) allow for multiple inputs at once (e.g. M1_spikes and M1_spikes_shift)
%   2) return predicted signals as field in trial_data, not as matrix
%
% INPUTS:
%   trial_data : the struct
%   in_var     : (string) name of input variable field
%   out_var    : (string) name out predicted variable field
%   params     : parameter struct
%     .fit_metric : which metric for evaluating fit ('corr','r2','vaf')
%     .num_dims   : how many dimensinos of in_var to use
%     .b          : weights from a previous call
%                       Used to use old model on new data
% OUTPUTS:
%   orig : time-varying original signal
%   pred : tim-varying predicted signal
%   c    : output of fit_metric
%   b    : weights computed for model
%
% Written by Matt Perich. Updated Feb 2017.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [orig,pred,c,b] = linPredSignal(trial_data,in_var,out_var,params)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fit_metric  =  'vaf';
num_dims    =  [];
b           =  [];
assignParams(who,params);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isempty(num_dims)
    num_dims = size(trial_data.(in_var),2);
end

if isempty(b)
    new_model = true;
end

in_signal = cat(1,trial_data.(in_var));
in_signal = [ones(size(in_signal,1),1), in_signal(:,1:num_dims)];
orig = cat(1,trial_data.(out_var));

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
            c(iDim) = compute_r2(out_signal,in_signal*b(:,iDim));
        case 'vaf'
            c(iDim) = compute_vaf(out_signal,in_signal*b(:,iDim));
    end
end

end