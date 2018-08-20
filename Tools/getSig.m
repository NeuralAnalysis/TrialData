%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function out = getSig(trial_data,signals)
%
%   Returns a matrix containing the requested signals for all time points.
% In effect, this is equivalent to cat(1,trial_data.SIGNAL), except it
% allows you to subselect columns of the signal, and also to concatenate
% multiple signals columnwise.
%
% INPUTS:
%   trial_data : the struct
%   signals    : which signals to use (either a string with a field name or
%                a Nx2 cell where N is the number of signals. The first
%                column is the field name and the second column is the
%                column indices that you want.
%
% OUTPUTS:
%   out        : a M x N matrix, where M is the total number of timepoints
%                across all trials and N is the total number of signals
%                requested.
%
% EXAMPLE:
%   Get the velocities for all trials
%       vel = getSig(trial_data,'vel');
%
%   Get the first 8 dimensions of M1 activity after PCA
%       m1 = getSig(trial_data,{'M1_pca',1:8});
%
%   Get all M1 and PMd activity
%       spikes = getSig(trial_data,{'M1_spikes','PMd_spikes'});
%
% Written by Matt Perich. Updated August 2018.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = getSig(trial_data,signals)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~iscell(signals)
    signals = {signals};
end
signals = check_signals(trial_data,signals);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% loop along signals and concatenate them
%   not the most efficient since I don't allocate memory but it should be
%   okay for reasonable file sizes
out = [];
for iSig = 1:size(signals,1)
    temp = cat(1,trial_data.(signals{iSig,1}));
    out = cat(2,out,temp(:,signals{iSig,2}));
end


