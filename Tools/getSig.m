%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function out = getSig(trial_data,signals)
%
%   Returns a matrix containing the requested signals for all time points.
% In effect, this is equivalent to cat(1,trial_data.SIGNAL), except it
% allows you to subselect columns of the signal, and also to concatenate
% multiple signals columnwise.
%
%   If signals is idx_ fields, it will return a matrix with all timepoints
% as the rows and each column is an idx_ field where 1 indicates an event
% happened and 0 indicates no event. This is useful for plotting events
% when matched with getSig for normal signals. Note that at the moment you
% cannot mix idx_ and regular signals in the same function call.
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

% check if any are idx_
if any(cellfun(@(x) ~isempty(regexp(x,'idx_','ONCE')),signals(:,1)')) 
    if ~all(cellfun(@(x) ~isempty(regexp(x,'idx_','ONCE')),signals(:,1)'))
        error('At the moment, getSig cannot mix idx_ and signals');
    end
    
    % preallocate, because with massively large datasets its worth it
    fn = getTDfields(trial_data,'time');
    if isempty(fn)
        error('Could not find any time varying signals to initialize the matrix');
    end
    out = zeros(size(cat(1,trial_data.(fn{1})),1),size(signals,1));
    % piece everything together
    for iSig = 1:size(signals,1)
        idx=[];
        t_total = 0;
        for trial = 1:length(trial_data)
            idx =  [idx, ...
                trial_data(trial).(signals{iSig}) + t_total];
            t_total = t_total + size(trial_data(trial).(fn{1}),1);
        end
        % ignore NaNs
        idx(isnan(idx)) = [];
        out(idx,iSig) = 1;
    end
    
    out = ~~out;
    
else % they are normal signals
    
    out = get_vars(trial_data,signals);
    
end





