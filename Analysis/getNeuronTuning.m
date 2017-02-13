%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [tcs,fr,covar] = getNeuronTuning(trial_data,which_method,params)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% will compute tuning curves using desired method
% A WORK IN PROGRESS
%
% INPUTS:
%   trial_data   : the struct
%   which_method : how to get the tuning
%          1) 'regress'  : cosine regression
%          2) 'vonmises' : regression of von mises
%          3) 'glm'      : glm of velocity
%   params       : parameter struct
%     .covariate : which covariate to use (e.g. 'target')
%     .array     : which array to use. Will do all neurons independently
%     .win       : window {'idx_OF_START',BINS_AFTER; 'idx_OF_END', BINS_AFTER}
%     .dt        : bin_size for computing firing rates
%
% To do:
%   1) implement glm and von mises methods
%   2) support multiple arrays
%   3) suppor arbitrary covariates, somehow?
%
% OUTPUTS:
%   tcs   : tuning curves [MEAN, MODULATION DEPTH, PREFERRED DIRECTION]
%   fr    : firing rates for each trial (rows) and neuron (cols)
%   covar : covariate (e.g. target direction) for each trial
%               last two outpus are the data used for tcs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEFAULT PARAMETERS
covariate  =  'target';
win        =  params.window;
array      =  params.array;
dt         =  0.01;
eval(structvars(length(fieldnames(params)),params)); %overwrite parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch covariate
    case 'target'
        covar = [trial_data.target_direction]';
    case 'movement'
        % find the direction the hand went in win
        error('not implemented');
end

switch lower(which_method)
    case 'regress'
        % build firing rate matrix for the specified window
        fr = zeros(length(trial_data),size(trial_data(1).([array '_spikes']),2));
        for trial = 1:length(trial_data)
            temp = trial_data(trial).([array '_spikes']);
            idx = trial_data(trial).(win{1,1})+win{1,2}:trial_data(trial).(win{2,1})+win{2,2};
            fr(trial,:) = sum(temp(idx,:),1)./(dt*sum(idx));
            
        end
        %fr(trial,:) = cellfun(@(x) sum(x(idx,:),1)./(dt*sum(idx)),{trial_data.([array '_spikes'])});
        
        % this is my cosine tuning function
        tcs = regressTuningCurves(fr,covar,{'none'},'doparallel',false);
    case 'vonmises'
        error('not implemented');
    case 'glm'
        error('not implemented.');
    otherwise
        error('Method not recognized');
end