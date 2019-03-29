%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [td] = addCorrelatedNoise(trial_data, params)
%
%   Add correlated noise to given signal
%
% INPUTS:
%   trial_data : the struct
%   params     : struct containing parameters
%     .signals : (cell) actual signals. Two options:
%                           1) {'NAME1','NAME2',etc}
%                           2) {'NAME1',idx; 'NAME2',idx; etc}
%                                   Here idx is which columns to use
%                                   Note: can use 'all' as idx for all
%     .noise_covar     : covariance matrix of noise to be added
%     .use_trials     : which trials to use
%
% OUTPUTS:
%   td : output trial_data structure with new field '*_noisy'
%
% EXAMPLES:
%
% Written by Raeed Chowdhury. Updated Oct 2017.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [td] = addCorrelatedNoise(trial_data,params)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEFAULT PARAMETER VALUES
signals         =  {};
noise_covar     = [];
use_trials      =  1:length(trial_data);
do_plot         =  false;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Some undocumented extra parameters
field_extra     =  {'_noisy'};
verbose         = false;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin > 1
    if isstruct(params) % overwrite parameters
        assignParams(who,params);
    else % assumes you  passed in signals
        signals = params;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
trial_data  =  check_td_quality(trial_data);
% Process inputs
td = trial_data(use_trials);
if isempty(signals)
    error('signals info must be provided');
end
signals = check_signals(td(1),signals);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isempty(noise_covar)
    warning('no noise covariance provided, using standard normal distribution')
    noise_covar = eye(sum(cellfun(@length,signals(:,2))));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% check output field addition
field_extra  = check_field_extra(field_extra,signals);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  

sig_names = cellfun(@(x) strrep(x,'_spikes',''),signals(:,1),'uni',0);
for trial = 1:length(td)
    data = get_vars(td(trial),signals);

    noisy_data = mvnrnd(data,noise_covar);

    td(trial).([[sig_names{:}] field_extra{1}]) = noisy_data;
end
