%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function trial_data = convBasisFunc(trial_data, params)
%
%   Convolves raised cosine basis functions of specified width and spacing
% with any signal. Adds field for each of which_vars called
% which_vars_shift, where each column is the original signal convolved with
% one of the bases. If which_var has N columns, then they will be
% next to each other in which_vars_shift (size of _shift is N*rcb_n)
%
% Notes:
%   1) based on Pillow's makeBasis_postSpike code (in TrialData/util)
%   2) all basis funcs have amplitude of 1. Does NOT normalize by area
%   under curve or anything like that, so resultant convolution changes the
%   amplitude of the original signal depending on how wide the basis is
%
% INPUTS:
%   trial_data : the struct
%   params     : parameter struct (all are required)
%       .which_vars : name (or names, as cell array) of field of trial_data to convolve
%       .rcb_hpeaks : vector containing location of first and last RCB vectors
%       .rcb_b      : nonlinear stretching of basis axes (large be is more linear)
%       .rcb_n      : how many basis vectors
%       .flip_time  : flag to flip causality
%                       true: acausal; false: causal (standard)
%
% OUTPUTS:
%   trial_data : struct with fields added for convolved
%
% Written by Matt Perich. Updated March 2019.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function trial_data = convBasisFunc(trial_data,params)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEFAULTS
rcb_hpeaks    =  [];
rcb_b         =  [];
rcb_n         =  [];
flip_time     =  false;
signals       =  getTDfields(trial_data,'time'); % default to all time signals
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Some undocumented extra parameters
field_extra   =  {'_rcb'};   % if empty, defaults to input field name(s)
verbose       =  false;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin > 1
    assignParams(who,params); % overwrite parameters
    if isempty([rcb_hpeaks, rcb_b, rcb_n])
        error('Requires peaks/b/n parameters.');
    end
else
    error('Requires params input.');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
trial_data  =  check_td_quality(trial_data);
signals     =  check_signals(trial_data(1),signals);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% check output field addition
field_extra  = check_field_extra(field_extra,signals);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bin_size = trial_data(1).bin_size;

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

fn_time = getTDfields(trial_data,'time');

if rcb_n > 0
    % do the convolution for each trial
    for trial = 1:length(trial_data) % loop along trials
        for iSig = 1:size(signals,1) % loop along variables
            
            sig_name = signals{iSig,1};
            % check if it's an idx_ field
            if ~isempty(regexp(sig_name,'idx_','ONCE'))
                sig_name = sig_name(5:end);
                % make a time vector
                data = zeros(size(trial_data(trial).(fn_time{1}),1),1);
                data(trial_data(trial).(signals{iSig,1})) = 1;
            else % maybe later add check to ensure it's time varying
                data = trial_data(trial).(signals{iSig,1});
            end
            temp_conv = zeros(size(data,1),size(data,2)*size(b,2));
            for iFunc = 1:size(b,2) % loop along basis funcs
                for i = 1:size(data,2)
                    if flip_time % acausal
                        temp =  flipud(conv(flipud(data(:,i)),b(:,iFunc),'full'));
                        temp_conv(:,i+(iFunc-1)*size(data,2)) = temp(end-size(data,1)+1:end);
                    else % causal (default mode)
                        temp = conv(data(:,i),b(:,iFunc),'full');
                        temp_conv(:,i+(iFunc-1)*size(data,2)) = temp(1:size(data,1));
                    end
                end
            end
            trial_data(trial).([sig_name field_extra{1}]) = temp_conv;
        end
    end
else
    disp('Requested zero raised cosine basis vectors. Skipping...');
end

% restore logical order
trial_data = reorderTDfields(trial_data);