%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS:
%   params : parameter struct
%     .signals        : (cell) which signals. Two options:
%                           1) {'NAME1','NAME2',etc}
%                           2) {'NAME1',idx; 'NAME2',idx; etc}
%                                   Here idx is which columns to use
%                                   Note: can use 'all' as idx for all
%     .conditions     : cell array of idx for each condition
%     .decision_var   : string name of field for "decision"
%                           Default: 'target_direction'
%     .num_dims       : how many dPCA dimensions
%     .do_plot        : flag to make dPCA plot
%
% OUTPUTS:
%
%
% TO DO:
%   - add dPCA projections to trial table
%       Would be getting geometry from trial averaging and projecting
%       single trials (add a guide to say what each dPC is)
%   - support different numbers of decisions for different conditions?
%   - support an arbitrary number of condition structs, rather than one
%   plus the decision_var?
%   - figure out naming parameters and plotting parameters, make them
%   writeable
%
% Written by Matt Perich. Adapted from Juan's code. Updated Feb 2017.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [trial_data, dPCA_info] = getDPCA(trial_data,params)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEFAULT PARAMETERS
signals = [];
conditions = [];
decision_var = 'target_direction';
num_dims   = 15;
do_plot = true;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Some undocumented extra functions
dpca_plot_fcn = @dpca_plot_default;
assignParams(who,params);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~iscell(conditions), conditions = {conditions}; end
signals = check_signals(trial_data(1),signals);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% THESE SHOULD ALL BECOME INPUT-BASED IN SOME INTUITIVE WAY
% ------------------------------------------------------------------------
% 2. Define parameters
% 1 - condition
% 2 - decision_var
% 3 - time
% [1 3] - condition/time interaction
% [2 3] - decision_var/time interaction
% [1 2] - condition/decision_var interaction
% [1 2 3] - rest
combined_params = { {1,[1 3]}, {2,[2,3]}, {3}, {[1 2],[1 2 3]} };
marg_names      = {'task','target','time','task/target interaction'};
marg_colors     = [23 100 171; 187 20 25; 150 150 150; 114 97 171]/256; % blue, red, grey, purple

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
u = unique([trial_data.(decision_var)]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get the size of the massive input matrix
%   N is the number of signals (e.g. neurons)
%   S is the number of conditions --> tasks in our case
%   D is the number of decisions --> targets in our case
%   T is the number of time points --each trial should have the same duration in time !!!
N = sum(cellfun(@(x) length(x),signals(:,2)));
S = length(conditions);
D = length(u);
T = size(trial_data(1).(signals{1,1}),1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% max number of repetitions
temp = zeros(S,D);
for s = 1:S
    for d = 1:D
        temp(s,d) = length(getTDidx(trial_data(conditions{s}),decision_var,u(d)));
    end
end
max_trial_num = max(max(temp));

% firing_rates: N x S x D x T x max_trial_num -- populated with our
% single_trial_data
firing_rates        = nan(N,S,D,T,max_trial_num);
% Get firing rates
for s = 1:S
    for d = 1:D
        temp = [];
        for n = 1:size(signals,1)
            [~,td] = getTDidx(trial_data(conditions{s}),decision_var,u(d));
            temp = cat(2,temp,cat(3,td.(signals{n,1})));
        end
        for n = 1:size(temp,2)
            firing_rates(n,s,d,:,1:size(temp,3)) = squeeze(temp(:,n,:));
        end
    end
end

% firing_rates_average: N x S x D x T -- these are PSTHs
firing_rates_avg    = nanmean(firing_rates, 5);

% ------------------------------------------------------------------------
% 2. Do dPCA without regularization

[W, V, which_marg]  = dpca( firing_rates_avg, num_dims, 'combinedParams', combined_params );

expl_var            = dpca_explainedVariance(firing_rates_avg, W, V, ...
    'combinedParams', combined_params);

time_events     = 1;
time            = (1:T)*trial_data(1).bin_size;
if do_plot
    dpca_plot(firing_rates_avg, W, V, dpca_plot_fcn, ...
        'explainedVar', expl_var, ...
        'marginalizationNames', marg_names, ...
        'marginalizationColours', marg_colors, ...
        'whichMarg', which_marg,                 ...
        'time', time,                        ...
        'timeEvents', time_events,               ...
        'timeMarginalization', 3, ...
        'legendSubplot', num_dims);
end


% ------------------------------------------------------------------------
% 2. Do dPCA in each marginalization separately

% dpca_perMarginalization(firing_rates_avg, @dpca_plot_default, ...
%    'combinedParams', combined_params);


% ------------------------------------------------------------------------
% project data onto dPC axes
[lat_vars, lat_vars_st] = get_lat_vars_dPCA( firing_rates_avg, firing_rates, W, ...
    'explainedVar', expl_var, ...
    'marginalizationNames', marg_names, ...
    'marginalizationColours', marg_colors, ...
    'whichMarg', which_marg,                 ...
    'time', time,                        ...
    'timeEvents', time_events,               ...
    'timeMarginalization', 3, ...
    'legendSubplot', num_dims);


% ------------------------------------------------------------------------
% 3. Return dPCA results

dPCA_info.W          = W;
dPCA_info.V          = V;
dPCA_info.lat_vars_mn = lat_vars;
dPCA_info.lat_vars   = lat_vars_st;
dPCA_info.which_marg = which_marg;
dPCA_info.marg_names = marg_names;
dPCA_info.expl_var   = expl_var;
dPCA_info.num_comps  = num_dims;
dPCA_info.combined_params = combined_params;


