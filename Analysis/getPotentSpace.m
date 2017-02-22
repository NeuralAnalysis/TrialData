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
%     .in_dims     : dimensionality of input space
%     .out_dims    : dimensionality of output space
%     .in_idx      : list of columns to include for in_signal (default: all)
%     .out_idx     : list of columns to include for out_signal (default: all)
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
%       .mu_in    : means of input signals
%       .mu_out   : means of output signals
%                NOTE: as with PCA, the mu are essential for V/w later
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
in_idx       =  [];
out_idx      =  [];
use_trials   =  1:length(trial_data);
assignParams(who,params); % overwrite parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isempty(in_signals), error('Need to specify input signals name'); end
if isempty(out_signals), error('Need to specify output signals name'); end
if isempty(in_dims), error('Need to specify input dimensionality'); end
if isempty(out_dims), error('Need to specify output dimensionality'); end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~iscell(in_signals), in_signals = {in_signals}; end
if ~iscell(out_signals), out_signals = {out_signals}; end
if isempty(in_idx)
    in_idx = cell(1,length(in_signals));
    for i = 1:length(in_signals), in_idx{i} = 1:size(trial_data(1).(in_signals{i}),2); end
end
if isempty(out_idx)
    out_idx = cell(1,length(out_signals));
    for i = 1:length(out_signals), out_idx{i} = 1:size(trial_data(1).(out_signals{i}),2); end
end
if ~iscell(in_idx), in_idx = {in_idx}; end
if ~iscell(out_idx), out_idx = {out_idx}; end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pca_params = params;

% get output PC space
pca_params.signals = out_signals;
pca_params.signal_idx = out_idx;
[~,pca_info] = getPCA(trial_data(use_trials),pca_params);
w_out = pca_info.w;
mu_out = pca_info.mu;
score_out = pca_info.scores;

% get input PC space
pca_params.signals = in_signals;
pca_params.signal_idx = in_idx;
[~,pca_info] = getPCA(trial_data(use_trials),pca_params);
w_in = pca_info.w;
mu_in = pca_info.mu;
score_in = pca_info.scores;
pca_params = pca_info.params;

if ~strcmpi(in_signals,out_signals)
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

% package up PCA weights, etc.
pca_params = rmfield(pca_params,{'signals','signal_idx'});
pca_info = struct(        ...
    'V_potent', V_potent, ...
    'V_null',   V_null,   ...
    'w_in',     w_in,     ...
    'w_out',    w_out,    ...
    'mu_in',    mu_in,    ...
    'mu_out',   mu_out,   ...
    'params',   pca_params);

% populate trial_data struct with PCA scores and null/potent projections
for trial = 1:length(trial_data)
    % get signals to recreate PCA
    temp = get_data(trial_data(trial),out_signals,out_idx);
    temp_sig_out = cellfun(@(x) strrep(x,'_spikes',''),out_signals,'uni',0);
    trial_data(trial).([temp_sig_out{:} '_pca']) = temp*w_out;
    
    if ~strcmpi(in_signals,out_signals)
        temp = get_data(trial_data(trial),in_signals,in_idx);
        temp_sig_in = cellfun(@(x) strrep(x,'_spikes',''),in_signals,'uni',0);
        trial_data(trial).([temp_sig_in{:} '_pca']) = temp*w_in;
        
        temp = trial_data(trial).([temp_sig_in{:} '_pca']);
        temp = temp(:,1:in_dims);
        trial_data(trial).([temp_sig_in{:} temp_sig_out{:} '_potent']) = temp*V_potent;
        trial_data(trial).([temp_sig_in{:} temp_sig_out{:} '_null']) = temp*V_null;
    end
end

end

function data = get_data(td,signals,signal_idx)
data = [];
for i = 1:length(signals)
    temp_data = cat(1,td.(signals{i}));
    data = [data, temp_data(:,signal_idx{i})];
end

end


