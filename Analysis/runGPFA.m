%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [trial_data, gpfa_out] = runGPFA(trial_data, params)
%
%   Run Gaussian Process Factor Analysis (GPFA) for a set of neural data.
% 
% INPUTS:
% Trial_data is a struct array where each element is a trial.
% Has the following required fields:
%   trialID      : a unique trial ID number
%   ARRAY_spikes : NxT array where N is number of neurons and T is number of time bins
%                    Each element is a count of binned spikes. ARRAY is currently 'M1' and/or 'PMd'
% 
% Fields for params input struct:
%   .arrays      : which arrays to use, e.g. 'M1' (put in cell for multiple)
%   .save_dir    : directory to save GPFA results in Byron's code. Pass in [] to skip (default)
%   .method      : gpfa, pca, etc. See Byron's code
%   .xdim        : assumed number of latent dimensions (default to 8, find optimal using CV)
%   .kernsd      : kernal width in s (default to 0.03, find optimal using CV)
%   .bin_w       : bin size desired for GPFA in s (default to 0.02)
%
% OUTPUTS:
% Field for gpfa_out output struct:
%   gpfa_out.trajectories          : trajectory information by trials (seqTrain in Byron's code)
%   gpfa_out.model                 : model parameter estimates (estParams in Byron's code)
%   gpfa_out.params.method         : 'gpfa','pca', etc (See Byron's code)
%   gpfa_out.params.arrays         : which arrays to use (copy of input)
%   gpfa_out.params.xdim           : how many assumed latent dimensions (copy of input)
%   gpfa_out.params.kernsd         : kernel width (copy of input)
%   gpfa_out.params.bin_width      : GPFA bin width (copy of input)
%
% Based on the GPFA Matlab codepack by Byron Yu. http://users.ece.cmu.edu/~byronyu/software.shtml
%   Note that I made some modifications to their code to make it work better with our data
%   These functions are copied in TrialData/util and have _td appended
% 
% Written by Matt Perich 12/2015. Updated 02/2017.
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [trial_data, gpfa_out] = runGPFA(trial_data, params)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isfield(params,'arrays'), arrays = params.arrays; else, error('Arrays not specified.'); end
save_dir    =  [];
method      =  'gpfa';
xDim        =  8;
kernSD      =  0.03;
bin_w       =  0.02;
numFolds    = 0;
runid       = '';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Some undocumented extra parameters
verbose = false;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
assignParams(who,params); % overwrite parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
trial_data  =  check_td_quality(trial_data);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~iscell(arrays), arrays = {arrays}; end
% get the desired spikes
dat = repmat(struct(),size(trial_data));
for iTrial = 1:length(dat)
    dat(iTrial).trialId = iTrial;
    temp = [];
    for arraynum = 1:length(arrays)
        temp = [temp, trial_data(iTrial).([arrays{arraynum} '_spikes'])];
    end
    dat(iTrial).spikes = temp';
end

% Basic extraction of neural trajectories
runIdx = ['-' [arrays{:}] runid];

% Extract neural trajectories
result = neuralTraj_td(runIdx, dat, save_dir, 'method', method, 'xDim', xDim,...
    'kernSDList', 1000*kernSD,'binWidth',1000*bin_w,'dataBinWidth',1000*trial_data(1).bin_size,'numFolds',numFolds);

% Orthonormalize neural trajectories
[estParams, seqTrain] = postprocess(result, 'kernSD', 1000*kernSD);

% Seq train will NOT be ordered by trial. This makes labeling by trial
% with indices quite difficult, so we reorder it to be sequential. We also
% keep the original (orig_SeqTrain) in case a shuffled form is desired
ord_seqTrain = seqTrain;
for tr = 1:length(seqTrain)
    tId = seqTrain(tr).trialId;
    ord_seqTrain(tId) = seqTrain(tr);
end
seqTrain = ord_seqTrain;

% assign useful things to output struct
gpfa_out.trajectories = seqTrain;
gpfa_out.model = estParams;
gpfa_out.params.method = method;
gpfa_out.params.arrays = arrays;
gpfa_out.params.xdim = xDim;
gpfa_out.params.kernsd = kernSD;
gpfa_out.params.bin_width = bin_w;

% add trajectory information to trial_data struct
for iTrial = 1:length(trial_data)
    trial_data(iTrial).([[arrays{:}] '_' method]) = seqTrain(iTrial).xorth';
end
