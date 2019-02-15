%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function out = getSigByTrial(trial_data,signals)
%
%   Returns a matrix containing the requested signals for all time points.
% In effect, this is equivalent to cat(3,trial_data.SIGNAL), except it
% allows you to subselect columns of the signal, and also to concatenate
% multiple signals columnwise.
%
%   Note that this requires every trial to be the same length! You can
% easily do this with trimTD or stretchSignals.
%
% INPUTS:
%   trial_data : the struct
%   signals    : which signals to use (either a string with a field name or
%                a Nx2 cell where N is the number of signals. The first
%                column is the field name and the second column is the
%                column indices that you want.
%
% OUTPUTS:
%   out        : a M x N x T matrix, where M is the total number of timepoints
%                across all trials and N is the total number of signals
%                requested and T is the total number of trials.
%
% EXAMPLE:
%   Get the velocities for all trials
%       vel = getSigByTrial(trial_data,'vel');
%
%   Get the first 8 dimensions of M1 activity after PCA
%       m1 = getSigByTrial(trial_data,{'M1_pca',1:8});
%
%   Get all M1 and PMd activity
%       spikes = getSigByTrial(trial_data,{'M1_spikes','PMd_spikes'});
%
% Written by Matt Perich. Updated August 2018.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = getSigByTrial(trial_data,signals)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~iscell(signals)
    signals = {signals};
end
signals = check_signals(trial_data,signals);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if length(unique(cellfun(@(x) size(x,1),{trial_data.(signals{1,1})}))) > 1
    error('Trials must all be the same length for getSigByTrial to work');
end

% initialize
out = zeros(size(trial_data(1).(signals{1,1}),1), ...
    length(cat(2,signals{:,2})), ...
    length(trial_data));

count = 0;
for iSig = 1:size(signals,1)
    temp = cat(3,trial_data.(signals{iSig,1}));
    out(:,count+1:count+length(signals{iSig,2}),:) = temp(:,signals{iSig,2},:);
    count = count + length(signals{iSig,2});
end







