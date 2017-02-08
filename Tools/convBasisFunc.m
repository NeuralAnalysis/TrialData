function trial_data = convBasisFunc(trial_data,which_vars,params)
% will convolve basis_func with which_vars fields of trial_data
%
% basis_funcs: {NUMBER,  % how many basis funcs
%               SPACING} % angle (typically pi/2)
dt = params.dt;
rcb_hpeaks = params.rcb_hpeaks;
rcb_b = params.rcb_b;
rcb_n = params.unit_lags;

if ~iscell(which_vars)
    which_vars = {which_vars};
end

% Inputs:
%     prs = param structure with fields:
%            ncols = # of basis vectors
%            hpeaks = 2-vector containg [1st_peak  last_peak], the peak
%                     location of first and last raised cosine basis vectors
%            b = offset for nonlinear stretching of x axis:  y = log(x+b)
%                     (larger b -> more nearly linear stretching)
%            absref = absolute refractory period (optional).  If specified,
%                     this param creates an additional "square" basis
%                     vector with support n [0,absref] (i.e., with a hard
%                     cutoff at absref)
%
%     dt = grid of time points for representing basis
%     iht (optional) = cut off time (or extend) basis so it matches this
%  Outputs:  iht = time lattice on which basis is defined
%            ihbas = orthogonalized basis
%            ihbasis = original (non-orthogonal) basis
%
%  ihbasprs.ncols = 5;
%  ihbasprs.hpeaks = [.1 2];
%  ihbasprs.b = .5;
%  ihbasprs.absref = .1;  %% (optional)
%  [iht,ihbas,ihbasis] = makeBasis_PostSpike(ihprs,dt);
[~, ~, b] = makeBasis_PostSpike(struct('ncols',rcb_n,'hpeaks',rcb_hpeaks,'b',rcb_b),dt);


% compute the raised cosine basis function
for iTrial = 1:length(trial_data) % loop along trials
    for iVar = 1:length(which_vars)
        temp = trial_data(iTrial).(which_vars{iVar});
        temp_conv = zeros(size(temp,1),size(temp,2)*size(b,2));
        for iFunc = 1:size(b,2)
            for i = 1:size(temp,2)
                temp_conv(:,i+(iFunc-1)*size(temp,2)) = conv(temp(:,i),b(:,iFunc),'same');
            end
        end
        trial_data(iTrial).([which_vars{iVar} '_shift']) = temp_conv;
    end
end
