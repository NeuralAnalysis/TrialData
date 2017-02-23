function r2 = compute_r2(sig1,sig2,which_method)
%   sig1         : vector of real signal
%   sig2         : vector of model fit signal
%   which_method : (string) how to calculate
%     1) 'corr' : use corrcoef (e.g. 'r', but squared)
if nargin < 3
    which_method = 'corr';
end

numSigs = size(sig1,2);
r2 = zeros(numSigs,1);

switch lower(which_method)
    case 'r'
        for i = 1:numSigs
            %Calculate R2
            R=corrcoef(sig1(:,i),sig2(:,i));
            r2(i)=R(1,2).^2;
            clear R
        end
        
    otherwise
        error('which_method not recognized.');
end

end