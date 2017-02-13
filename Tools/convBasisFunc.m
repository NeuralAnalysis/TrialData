%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function trial_data = convBasisFunc(trial_data, which_vars, params)
% 
%   Convolves raised cosine basis functions of specified width and spacing
% with any signal. Adds field for each of which_vars called
% which_vars_shift, where each column is the original signal convolved with
% one of the bases. If which_var has N columns, then they will be
% next to each other in which_vars_shift (size of _shift is N*rcb_n)
%
% Note: based on Pillow's makeBasis_postSpike code (in TrialData/util).
%
% INPUTS:
%   trial_data : the struct
%   which_vars : name (or names, as cell array) of field of trial_data to convolve
%   params     : parameter struct (all are required)
%       .rcb_hpeaks : vector containing location of first and last RCB vectors
%       .rcb_b      : nonlinear stretching of basis axes (large be is more linear)
%       .rcb_n      : how many basis vectors
% 
% OUTPUTS:
%   trial_data : struct with fields added for convolved
% 
% Written by Matt Perich. Updated Feb 2017.
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function trial_data = convBasisFunc(trial_data,which_vars,params)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rcb_hpeaks = params.rcb_hpeaks;
rcb_b = params.rcb_b;
rcb_n = params.rcb_n;

bin_size = trial_data(1).bin_size;

if ~iscell(which_vars), which_vars = {which_vars}; end

% Here is the documentation from Pillow's code so you know what's going on
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
[~, ~, b] = makeBasis_PostSpike(struct('ncols',rcb_n,'hpeaks',rcb_hpeaks,'b',rcb_b),bin_size);


% do the convolution for each trial
for iTrial = 1:length(trial_data) % loop along trials
    for iVar = 1:length(which_vars) % loop along variables
        temp = trial_data(iTrial).(which_vars{iVar});
        temp_conv = zeros(size(temp,1),size(temp,2)*size(b,2));
        for iFunc = 1:size(b,2) % loop along basis funcs
            for i = 1:size(temp,2)
                temp_conv(:,i+(iFunc-1)*size(temp,2)) = conv(temp(:,i),b(:,iFunc),'same');
            end
        end
        trial_data(iTrial).([which_vars{iVar} '_shift']) = temp_conv;
    end
end
