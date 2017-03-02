%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [trial_data, pca_info] = getPotentSpace(trial_data, params)
%
%   This will compute the geometry of the potent and null spaces for two
% sets of signals. Calls getPCA for neural data processing.
%
% TO DO:
%   - allow for multiple input signals (e.g. M1+PMd arrays onto vel)
%
% INPUTS:
%   trial_data : the struct
%   params    : struct containing parameters
%     .in_signals  : (string) name of input signal
%     .out_signals : (string) name of output signal
%                       See getPCA for formatting. Can include column idx.
%     .in_dims     : dimensionality of input space
%     .out_dims    : dimensionality of output space
%     .use_trials  : vector list of trial indices to use for spaces
%         NOTE: if adds scores to trial_data, will add scores for all, not only use_trials
%     (Also include any parameters desired for getPCA, as params gets
%     passed to that function call as well)
%
% OUTPUTS:
%   trial_data : the struct with potent/null fields added
%                NOTE: processed by getPCA too (e.g. trial averaging, etc)
%   pca_results : struct with PCA output
%       .V_potent : potent space basis vectors
%       .V_null   : null space basis vectors
%       .w_in     : weights for PC space of input signals
%       .w_out    : weights for PC space of output signals
%       .params   : parameter struct
%
% Written by Matt Perich. Updated Feb 2017.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [trial_data,pca_info] = getPotentSpace(trial_data,params)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEFAULT PARAMETERS
in_signals   =  [];
out_signals  =  [];
in_dims      =  [];
out_dims     =  [];
use_trials   =  1:length(trial_data);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% other undocumented PCA parameters
pca_centered   =  true;   % whether to center data
pca_algorithm  =  'svd';  % which PCA algorithm
add_proj_to_td =  true;   % add projections to trial data
assignParams(who,params); % overwrite parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isempty(in_signals), error('Need to specify input signals'); end
if isempty(out_signals), error('Need to specify output signals'); end
if isempty(in_dims), error('Need to specify input dimensionality'); end
if isempty(out_dims), error('Need to specify output dimensionality'); end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
in_signals = check_signals(trial_data(1),in_signals);
out_signals = check_signals(trial_data(1),out_signals);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pca_params = params;
pca_params.pca_centered = pca_centered;
pca_params.pca_algorithm = pca_algorithm;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get input PC space
pca_params.signals = in_signals;
[~,pca_info_in] = getPCA(trial_data(use_trials),pca_params);
w_in = pca_info_in.w;
score_in = pca_info_in.scores;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get output PC space
pca_params.signals = out_signals;
[~,pca_info_out] = getPCA(trial_data(use_trials),pca_params);
w_out = pca_info_out.w;
score_out = pca_info_out.scores;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find null and potent spaces
sig_name_out = cellfun(@(x) strrep(x,'_spikes',''),out_signals(:,1),'uni',0);
sig_name_in = cellfun(@(x) strrep(x,'_spikes',''),in_signals(:,1),'uni',0);
if ~strcmpi([sig_name_in{:}],[sig_name_out{:}])
    y = score_out(:,1:out_dims);
    x = score_in(:,1:in_dims);
    % find the model
    W = zeros( size(y,2), size(x,2) );
    for i = 1:size(y,2)
        %[b_pc, ~, ~, ~, stats_this] = regress(y(:,i),x);
        b_pc = x\y(:,i);
        % fill MIMO matrix W
        W(i,:) = b_pc';
    end
    
    % do SVD of weights to get potent/null spaces
    [U, S, V]                   = svd( W );
    % The output potent spaces is defined by the first m columns of V', where m
    % is the number of dimensions of the output
    V_potent                    = V(1:size(y,2),:)';
    V_null                      = V(size(y,2)+1:end,:)';
else
    disp('Input and output signals are the same. Skipping null/potent.');
    [V_potent,V_null] = deal([]);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% package up PCA weights, etc.
out_params = struct( ...
    'in_signals', {in_signals}, ...
    'out_signals', {out_signals}, ...
    'in_dims', in_dims, ...
    'out_dims', out_dims, ...
    'use_trials', use_trials);
pca_info = struct(        ...
    'V_potent', V_potent, ...
    'V_null',   V_null,   ...
    'w_in',     w_in,     ...
    'w_out',    w_out,    ...
    'params',   out_params);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% populate trial_data struct with PCA scores and null/potent projections
if add_proj_to_td
    if pca_centered
        mu_in  = mean(get_vars(trial_data,in_signals),1);
        mu_out = mean(get_vars(trial_data,out_signals),1);
    else
        n_signals_in  = cellfun(@(x) length(x),in_signals(:,2));
        n_signals_out = cellfun(@(x) length(x),out_signals(:,2));
        mu_in  = zeros(1,sum(n_signals_in));
        mu_out = zeros(1,sum(n_signals_out));
    end
    for trial = 1:length(trial_data)
        % get signals to recreate input PCA
        data = get_vars(trial_data(trial),out_signals);
        trial_data(trial).([sig_name_out{:} '_pca']) = (data - repmat(mu_out,size(data,1),1)) * w_out;
        
        if ~strcmpi(sig_name_in{:},sig_name_out{:}) % if they aren't the same
            % now add output PCA
            data = get_vars(trial_data(trial),in_signals);
            trial_data(trial).([sig_name_in{:} '_pca']) = (data - repmat(mu_in,size(data,1),1)) * w_in;
            
            % now do null/potent
            data = trial_data(trial).([sig_name_in{:} '_pca']);
            data = data(:,1:in_dims);
            trial_data(trial).([sig_name_in{:} sig_name_out{:} '_potent']) = data * V_potent;
            trial_data(trial).([sig_name_in{:} sig_name_out{:} '_null'  ]) = data * V_null;
        end
    end
end

