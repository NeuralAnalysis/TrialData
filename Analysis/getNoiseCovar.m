%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [noise_covar,td] = getNoiseCovar(trial_data, params)
%
% [noise_covar] = getNoiseCovar(trial_data, params);
%   Computes covariance matrix of noise, given actual signals and model
% predicted signals.
%
% INPUTS:
%   trial_data : the struct
%   params     : struct containing parameters
%     .actual_signals : (cell) actual signals. Two options:
%                           1) {'NAME1','NAME2',etc}
%                           2) {'NAME1',idx; 'NAME2',idx; etc}
%                                   Here idx is which columns to use
%                                   Note: can use 'all' as idx for all
%     .modeled_signals : (cell) modeled signals. Two options:
%                           1) {'NAME1','NAME2',etc}
%                           2) {'NAME1',idx; 'NAME2',idx; etc}
%                                   Here idx is which columns to use
%                                   Note: can use 'all' as idx for all
%     .model_name     : unique name under which residuals will be put
%                           (will be postfixed by '_residuals')
%     .use_trials     : which trials to use for noise_covar (default: all)
%     .do_plot        : flag to make pairwiseCorr plot (default: false)
%
% OUTPUTS:
%   noise_covar : noise covariance matrix for given signal
%   td : output trial_data structure with new field '*_residuals'
%
% EXAMPLES:
%
% Written by Raeed Chowdhury. Updated Oct 2017.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [noise_covar,td] = getNoiseCovar(trial_data,params)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEFAULT PARAMETER VALUES
actual_signals  =  {};
modeled_signals =  {};
model_name      =  'default';
use_trials      =  1:length(trial_data);
do_plot         =  false;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Some extra parameters you can change that aren't described in header
verbose = false;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin > 1, assignParams(who,params); end % overwrite parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Process inputs
if ~isstruct(trial_data), error('First input must be trial_data struct!'); end

if isempty(actual_signals) || isempty(modeled_signals)
    error('actual/modeled info must be provided');
end

actual_signals = check_signals(trial_data(1),actual_signals);
modeled_signals = check_signals(trial_data(1),modeled_signals);
trial_data = trial_data(use_trials);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

td = trial_data;
fn = [model_name,'_residuals'];
for i = 1:length(td)
    td(i).(fn) = get_vars(td(i),actual_signals)-get_vars(td(i),modeled_signals);
end

% build master matrix of spiking
noise_signals = check_signals(td,{fn});
data = get_vars(td,noise_signals);

% get pairwise covariance
noise_covar = cov(data);

if do_plot
    figure
    imagesc(noise_covar)
    colorbar
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
